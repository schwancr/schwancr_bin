#!/usr/bin/perl -w
use strict;

#####
# 9/8/06 capNewFF.pl
# input: uncapped pdb file, complete pdb (contains at least the residues that are one past the c-terminal residues)
# output 1: capped pdb file (RNA nucleotides become terminal nucleotides, and amino acids get n-term ACE, c-term NH2
# output 2: list of chains and their break points
# fix: doesn't handle nucleotides/aa's by themselves (i.e. nuc/aa is both nterm and cterm)	
# atom numbers not correct because original file not correct, gromacs will handle
# fix: put in a flag for missing atoms: because this file caps at numeric chain breaks, you have to manually
# check to see if any breaks occur in missing atom regions 
# note: must run rename.pl if get error with sub lookupCapName
#
# Q.  how is this different than cap.pl?
#	1.  fixed bug: the additional atom added in the c-cap had the name of the original atom
#	instead of the new atom
#	2.  atom/res names reflect new amber ff
#		a.  res name (e.g. CALA instead of ALC)
#		b.  atom names (e.g. OC1 and OC2 for c-term)
# todo: fix HIS to HIP (protonated form of HIS), and LYS to LIP (protonated form of LYS)
# KNOWN BUG: if there are two lones in a row, the second one gets mislabeled as a c-term
# KNOWN BUG: if your "complete" pdb contains a concatanation of two or more pdb files (e.g. large and small subunits
# of the ribosome) and the individual files have overlapping chain names, there may be problems with cterm amino acids: 
# the additional atom will come from the +1 residue of the first chain that matches the chain name, even if it's from the
# wrong individual file 
###

if (scalar @ARGV != 2) {die "usage:\nperl cutoutPoints2pdb.pl <cutout.pdb> <complete.pdb>, where complete.pdb contains at least the residues that extend each cap by one\n";};

my $cutoutFilename = $ARGV[0];
my $outputFilename = $cutoutFilename.".capped";
my $capPlusOneFile = $ARGV[1];
my $capLocOutFile = $cutoutFilename.".capLocations";

open (my $inputFH, "<$cutoutFilename") or die "can't open file $cutoutFilename\n";
open (my $outputFH, ">$outputFilename") or die "can't open file $outputFilename\n";
open (my $capLocOutFH, ">$capLocOutFile") or die "can't open file $outputFilename\n";

my $atomCounter = 0;
my $previousResidueNum = 0;
my $previousChain = "";
my $currentResidueNum = 0;
my $currentChain = "";
my $atomsInResidueCounter = 0;
my @currentResidue;
my @previousResidue;
my $previousResidueStatus = "mid";
my $n_cap = 0;
my $previousResidueName = "";
my $currentResidueName = "";
my @outputLinesNonCaps;
my %linesCaps; #keys will be ch#resNum## where # is the chain and ## is the residue number

while (my $line = <$inputFH>){
	if ($line =~ m/\AATOM/ || $line =~ m/\AHETATM/) {
		## initialize variables		
		my $atomNum = substr($line, 6, 5);
		$atomNum =~ s/\s+//;
		my $residueNum = substr($line, 22, 4);
		$residueNum =~ s/\s+//;
		$currentResidueNum = $residueNum;
		my $chain = substr($line, 21, 1);
		$chain =~ s/\s+//;
		$previousChain = $currentChain;
		$currentChain = $chain;
		my $residueName = substr($line, 17, 3);
		if ($residueName eq "HOH" || $residueName eq "MG " || $residueName eq "CL " || $residueName eq "NA " || $residueName eq "K  " || $residueName eq "CD ")
		{	
			next;
		}
		$residueName =~ s/\s+//;
		$previousResidueName = $currentResidueName;
		$currentResidueName = $residueName;
		chomp($line);
		## if this atom is within the current residue, update @currentResidue
		## elsif this atom is the first one of the next residue, 
		##	1.  copy @currentResidue to @previousResidue
		##	2.  clear @currentResidue
		##	3.  add $line to $currentResidue[0]
		##
		##	4.  if ($n_cap){$previousResidueStatus="n-cap"}else{$previousResidueStatus="mid"}
		##	5.  update @previousResidue, given $previousResidueStatus & $previousResidueName
		##	6.  print @previousResidue
		##
		##	7.  $n_cap = 0;
		##	8.  $previousResidueStatus = "";
		##	9.  $atomsInResidueCounter = 0;
		
		## elsif this atom is the first one of a residue that is not "next" in seq order
		##	1.  copy @currentResidue to @previousResidue	
		##	2.  clear @currentResidue
		##	3.  add $line to $currentResidue[0]
		##
		##	4.  algorithm check: residue is either n-term or chain break, n-cap on both conditions
		##		a.  if new chain: n-term
		##		b.  if old chain, not next in seq order: chain break
		##			i.  make sure it's due to end of box, not missing atoms
		##		c.  die otherwise
		##	
		##	5.  $previousResidueStatus = "c-cap"; 
		##		(@previousResidue is c-term or chain break, c-cap on both conditions)
		##	6.  update @previousResidue, given $previousResidueStatus & $previousResidueName
		##	7.  print @previousResidue
		##
		##	8.  $n_cap = 1; 	
		##	9.  $previousResidueStatus = "";	
		##	10.  $atomsInResidueCounter = 0;
	
		##  updating @currentResidue
		if ( ($currentResidueNum == $previousResidueNum) and
		     ($currentChain eq $previousChain)){
			$currentResidue[$atomsInResidueCounter] = $line;
		}
		##  this atom is the start of the next residue in the same chain
		elsif ( ($currentResidueNum == (1 + $previousResidueNum)) and
			($currentChain eq $previousChain)){
			startNewResidue($line);
			
			if ($n_cap) {$previousResidueStatus = "n-cap"}else{$previousResidueStatus = "mid"};
			my @updatedRes = updateResidue($previousResidueStatus, $previousResidueName, $previousResidueNum, $previousChain);
			putResidueInArray(@updatedRes);
	
			$n_cap = 0;
			#$previousResidueStatus = "";
			$atomsInResidueCounter = 0;
		##  this atom is the start of a residue that is not the next residue: CAP
		}elsif ( 
			 (($currentResidueNum !=  $previousResidueNum     ) and ($currentResidueNum != ($previousResidueNum + 1)))
			 or
			 ($currentChain ne $previousChain) 
			){
			startNewResidue($line);

			## new chain? 
			if ( $currentChain ne $previousChain ){ # we're at the n-term of a new chain 
				
			}elsif ( $currentChain eq $previousChain ){ #we're in the same chain, n-term of break
				
			}
			if (($previousResidueStatus eq "c-cap")||($previousResidueStatus eq "lone")){$previousResidueStatus = "lone"}
			else {$previousResidueStatus = "c-cap"};# CRS EDIT: "||($previousResidueStatus eq "lone")"

			my @updatedRes = updateResidue($previousResidueStatus, $previousResidueName, $previousResidueNum,$previousChain);
			putResidueInArray(@updatedRes);

			$n_cap = 1;
			$atomsInResidueCounter = 0;
			
		}else {die "GOT HERE";}

		## update counter variables
		$atomCounter++;
		$previousResidueNum = $currentResidueNum;
		$atomsInResidueCounter++;
	}
	if (eof($inputFH)){#when reach eof, print out the last residue
		if ( ($previousResidueStatus eq "mid") || ($previousResidueStatus eq "n-cap") ){$previousResidueStatus = "c-cap"}
		elsif ( ($previousResidueStatus eq "c-cap") || ($previousResidueStatus eq "lone") ){$previousResidueStatus = "lone"}
		else {die "error: at the end of the file, previousResidueStatus $previousResidueStatus must be mid, n-cap, or c-cap\n";}
		@previousResidue = @currentResidue;
		my @updatedRes = updateResidue($previousResidueStatus, $previousResidueName, $previousResidueNum,$previousChain);
		putResidueInArray(@updatedRes);
	}
}

capAndPrint();

close($inputFH);
close($outputFH);
close($capLocOutFile);

sub capAndPrint{
	my $numItems = scalar @outputLinesNonCaps;
	for (my $i = 0; $i < $numItems; $i++){
		if ($outputLinesNonCaps[$i] =~ m/n-cap/){
			#n-cap
			my @cappedResidue = nCap($outputLinesNonCaps[$i]);
			my $lengthRes = scalar @cappedResidue;
			print $capLocOutFH "n-cap goes here: $outputLinesNonCaps[$i]\n";
			for (my $i=0; $i<$lengthRes; $i++){
				print $outputFH "$cappedResidue[$i]\n";
			}
		}elsif ($outputLinesNonCaps[$i] =~ m/c-cap/){
			#c-cap
			my @cappedResidue = cCap($outputLinesNonCaps[$i]);
			my $lengthRes = scalar @cappedResidue;
			print $capLocOutFH "c-cap goes here: $outputLinesNonCaps[$i]\n";
			for (my $i=0; $i<$lengthRes; $i++){
				print $outputFH "$cappedResidue[$i]\n";
			}
		}elsif ($outputLinesNonCaps[$i] =~ m/lone/){
			#print "hit a lone for ";die;
			my @cappedResidue = loneCap($outputLinesNonCaps[$i]);
			my $lengthRes = scalar @cappedResidue;
			print $capLocOutFH "lone res goes here: $outputLinesNonCaps[$i]\n";
			for (my $i=0; $i<$lengthRes; $i++){
				print $outputFH "$cappedResidue[$i]\n";
			}
			
		}
		else {
			print $outputFH "$outputLinesNonCaps[$i]\n";
		}
	}
}

sub updateResidue{ #uses global variable @previousResidue
	my ($resStatus, $resName, $resNum, $chain) = @_;
	my @updatedResidue;
	my $resLength = scalar @previousResidue;
	my $linesCapsKey = "ch".$chain."resNum".$resNum;
	if ($resStatus eq "mid"){ #no cap needed
		@updatedResidue = @previousResidue;
	}
	elsif($resStatus eq "n-cap"){	
		#1.  look for $resName's n-cap equivalent
		#2.  change the resName to new name on each index of @previousResidue
		#3.  add additional atoms as needed
		for (my $i = 0; $i < $resLength; $i++){
			$linesCaps{$linesCapsKey}[$i]="$previousResidue[$i]";
		}
		$updatedResidue[0] = "n-cap $resNum $chain";
	}
	elsif($resStatus eq "c-cap"){
		if ( (scalar @previousResidue) > 1){# if it's not greater than one, we're at the first atom in file
			for (my $i = 0; $i < $resLength; $i++){
				$linesCaps{$linesCapsKey}[$i]="$previousResidue[$i]";
			}
			$updatedResidue[0] = "c-cap $resNum $chain";
		} #else {$updatedResidue[0] = "we're at first atom"}; #first atom in cutout
	}
	elsif($resStatus eq "lone"){#atom is by itself
		#print "LONE: resName $resName resNum $resNum chain $chain resStatus $resStatus\n";
		for (my $i = 0; $i < $resLength; $i++){
			$linesCaps{$linesCapsKey}[$i]="$previousResidue[$i]";
		}
		$updatedResidue[0] = "lone $resNum $chain";
	}
	else {die "error in sub updateResidue.  residueStatus must be mid/n-cap/c-cap/lone\n";} 
	return @updatedResidue;
}

sub putResidueInArray{
	my @residue = @_;
	my $residueLength = scalar @residue;
	for (my $i=0; $i < $residueLength; $i++){
		#print $outputFH "$residue[$i]\n";
		push @outputLinesNonCaps, $residue[$i];
	}
}

sub startNewResidue {#uses global variables @previousResidue and @currentResidue
	my ($line) = @_;
	@previousResidue = @currentResidue;
	@currentResidue = ();
	$currentResidue[0] = $line;
}
sub nCap{ #n-caps/5' ends new atoms are hydrogens (no heavy atoms), therefore rename only (no atoms location info needed) for aa, nucleotides: new hydrogen replaces three phosphate atoms (therefore delete three atoms)
	my ($line) = @_;
	my @lineItems = split " ", $line;
	my $resNum = $lineItems[1];
	my $chain = $lineItems[2];
	my $linesCapsKey = "ch".$chain."resNum".$resNum;
	my $firstAtom = $linesCaps{$linesCapsKey}[0];
	my $resName = substr($firstAtom, 17, 3);
	$resName =~ s/\s+//;
	my $cappedResName = lookupCapName($resName, "n");
	my @cappedRes = @{ $linesCaps{$linesCapsKey}  };
	my $numAtoms = scalar @cappedRes;
	
	if ( (length $resName) == 3){ #aa
		for (my $i=0; $i < $numAtoms; $i++){
			#substr($cappedRes[$i], 17, 3, $cappedResName);
			substr($cappedRes[$i], 17, 4, $cappedResName);
		}
	}elsif ( (length $resName) == 2) {#nuc
		## CRS: Find out the very first atom name: (if it's a P, then delete, else, don't)
		my $atomName = substr($cappedRes[0],13,1);
		
		if ( ($atomName eq 'P') ) {		
			for (my $i=3; $i < $numAtoms; $i++){
				substr($cappedRes[$i], 17, 3, $cappedResName);
				my $newIndex = $i - 3;
				$cappedRes[$newIndex] = $cappedRes[$i]; 
			}
		$#cappedRes = $numAtoms - 4; # $#array is the number of the last valid index number
		}else {
			for (my $i=0; $i < $numAtoms; $i++){
				substr($cappedRes[$i], 17, 3, $cappedResName);
			}
		$#cappedRes = $numAtoms - 1; # $#array is the number of the last valid index number
		
		}
		
	}else { 
		die "error in sub nCap.  resName $resName must be two (aa) or three (nuc) in length\n";
	}
	return @cappedRes;
}
sub cCap{ #c-caps/3' ends: new atoms are hydrogens (no heavy atoms) except for c-cap on aa, which has a new O2 atom which will go on th next atom's "N" position
	my ($line) = @_;
	my @lineItems = split " ", $line;
	my $resNum = $lineItems[1];
	my $chain = $lineItems[2];
	my $linesCapsKey = "ch".$chain."resNum".$resNum;
	my $firstAtom = $linesCaps{$linesCapsKey}[0];
	my $resName = substr($firstAtom, 17, 3);
	$resName =~ s/\s+//;
	my $cappedResName = lookupCapName($resName, "c");
	my @cappedRes = @{ $linesCaps{$linesCapsKey}  };
	my $numAtoms = scalar @cappedRes;
	if ((length $resName) == 2){#if nucleotide, new atoms are just hydrogen, all we need to do is rename
		for (my $i=0; $i < $numAtoms; $i++){
			substr($cappedRes[$i], 17, 3, $cappedResName);
		}
	}elsif ((length $resName) == 3){
		for (my $i=0; $i < $numAtoms; $i++){
			substr($cappedRes[$i], 17, 4, $cappedResName);
		}
	}
	
	elsif ((length $resName) == 3){#if aa, must add atom and get info about location
		my $capPositionResNum = $resNum + 1;
		my $resNumSize4 = $resNum;
		open (my $capLocationsFH, "<$capPlusOneFile") or die "can't open file $capPlusOneFile\n"; #this file contains atoms that are +/- 1 residue from the cap.  Atoms from these residues will be used as locations for additional atoms in caps
		LINE_LOOP: while (my $capPositionLine = <$capLocationsFH>){
			if ($capPositionLine =~ m/\AATOM/){
				my $currentLineChain = substr($capPositionLine, 21, 1);
				$currentLineChain =~ s/\s+//g; 
				my $currentLineResNum = substr($capPositionLine, 22, 4);
				$currentLineResNum =~ s/\s+//g;
				my $currentLineAtomType= substr($capPositionLine, 13, 3);
				$currentLineAtomType =~ s/\s+//g;
				chomp($capPositionLine);
				if ( ($currentLineChain eq $chain) && ($currentLineResNum == $capPositionResNum) && ($currentLineAtomType eq "N") ){
					my $lengthResNum = length $resNum;
					if ($resNum != 4){
						for (my $i=$lengthResNum; $i < 4; $i++){
							$resNumSize4 = " ".$resNumSize4;
						} 
					}
					substr($capPositionLine, 22, 4, $resNumSize4);
					substr($capPositionLine, 13, 3, "OC2");
					substr($capPositionLine, 76, 2, " O");
					substr($capPositionLine, 17, 4, $cappedResName);
					for (my $i=0; $i < $numAtoms; $i++){
						#substr($cappedRes[$i], 17, 3, $cappedResName);
						substr($cappedRes[$i], 17, 4, $cappedResName);
						if ($cappedRes[$i] =~ m/[0-9]  O /){
							substr($cappedRes[$i], 13, 3, "OC1");
						}
					}
					push @cappedRes, $capPositionLine;	
					close ($capLocationsFH);
					last LINE_LOOP;
				}else {
					next LINE_LOOP;
				}
			}	
		}
	}
	return @cappedRes;
}
sub loneCap{ #lone residues: aa's: need to get eric to make a new res type, in the meantime leave it the same as a "mid" aa; nucleotides: use RXN type
	my ($line) = @_;
	my @lineItems = split " ", $line;
	my $resNum = $lineItems[1];
	my $chain = $lineItems[2];
	my $linesCapsKey = "ch".$chain."resNum".$resNum;
	my $firstAtom = $linesCaps{$linesCapsKey}[0];
	my $resName = substr($firstAtom, 17, 3);
	$resName =~ s/\s+//;
	my $cappedResName = lookupCapName($resName, "lone");
	my @cappedRes = @{ $linesCaps{$linesCapsKey}  };
	my $numAtoms = scalar @cappedRes;
	
	if ( (length $resName) == 3){ #aa
		for (my $i=0; $i < $numAtoms; $i++){
			substr($cappedRes[$i], 17, 3, $cappedResName);
		}
	}elsif ( (length $resName) == 2) {#nuc
		for (my $i=3; $i < $numAtoms; $i++){
			substr($cappedRes[$i], 17, 3, $cappedResName);
			my $newIndex = $i - 3;
			$cappedRes[$newIndex] = $cappedRes[$i]; 
		}
		$#cappedRes = $numAtoms - 4; # $#array is the number of the last valid index number
	}else { 
		die "error in sub loneCap.  resName $resName must be two (aa) or three (nuc) in length\n";
	}
	return @cappedRes;
}


sub lookupCapName{
	my ($resName, $NorCorLone) = @_;
	my $switchKey = $NorCorLone."capResName".$resName;
	my %switch = (
		'ncapResNameRA'	=>	'RA5',
		'ccapResNameRA'	=>	'RA3',
		'lonecapResNameRA'	=>	'RAN',
		'ncapResNameRC'	=>	'RC5',
		'ccapResNameRC'	=>	'RC3',
		'lonecapResNameRC'	=>	'RCN',
		'ncapResNameRG'	=>	'RG5',
		'ccapResNameRG'	=>	'RG3',
		'lonecapResNameRG'	=>	'RGN',
		'ncapResNameRU'	=>	'RU5',
		'ccapResNameRU'	=>	'RU3',
		'lonecapResNameRU'	=>	'RUN',
		'ncapResNameALA'	=>	'NALA',
		'ccapResNameALA'	=>	'CALA',
		'lonecapResNameALA'	=>	'ALA',
		'ncapResNameGLY'	=>	'NGLY',
		'ccapResNameGLY'	=>	'CGLY',
		'lonecapResNameGLY'	=>	'GLY',
		'ncapResNameSER'	=>	'NSER',
		'ccapResNameSER'	=>	'CSER',
		'lonecapResNameSER'	=>	'SER',
		'ncapResNameTHR'	=>	'NTHR',
		'ccapResNameTHR'	=>	'CTHR',
		'lonecapResNameTHR'	=>	'THR',
		'ncapResNameLEU'	=>	'NLEU',
		'ccapResNameLEU'	=>	'CLEU',
		'lonecapResNameLEU'	=>	'LEU',
		'ncapResNameILE'	=>	'NILE',
		'ccapResNameILE'	=>	'CILE',
		'lonecapResNameILE'	=>	'ILE',
		'ncapResNameVAL'	=>	'NVAL',
		'ccapResNameVAL'	=>	'CVAL',
		'lonecapResNameVAL'	=>	'VAL',
		'ncapResNameASN'	=>	'NASN',
		'ccapResNameASN'	=>	'CASN',
		'lonecapResNameASN'	=>	'ASN',
		'ncapResNameGLN'	=>	'NGLN',
		'ccapResNameGLN'	=>	'CGLN',
		'lonecapResNameGLN'	=>	'GLN',
		'ncapResNameARG'	=>	'NARG',
		'ccapResNameARG'	=>	'CARG',
		'lonecapResNameARG'	=>	'ARG',
		'ncapResNameHIS'	=>	'NHIP', #fix: HIS: which one?
		'ccapResNameHIS'	=>	'CHIP',
		'lonecapResNameHIS'	=>	'HIP',
		'ncapResNameTRP'	=>	'NTRP',
		'ccapResNameTRP'	=>	'CTRP',
		'lonecapResNameTRP'	=>	'TRP',
		'ncapResNamePHE'	=>	'NPHE',
		'ccapResNamePHE'	=>	'CPHE',
		'lonecapResNamePHE'	=>	'PHE',
		'ncapResNameTYR'	=>	'NTYR',
		'ccapResNameTYR'	=>	'CTYR',
		'lonecapResNameTYR'	=>	'TYR',
		'ncapResNameGLU'	=>	'NGLU',
		'ccapResNameGLU'	=>	'CGLU',
		'lonecapResNameGLU'	=>	'GLU',
		'ncapResNameASP'	=>	'NASP',
		'ccapResNameASP'	=>	'CASP',
		'lonecapResNameASP'	=>	'ASP',
		'ncapResNameLYS'	=>	'NLYS',
		'ccapResNameLYS'	=>	'CLYS',
		'lonecapResNameLYS'	=>	'LYS',
		'ncapResNamePRO'	=>	'NPRO',
		'ccapResNamePRO'	=>	'CPRO',
		'lonecapResNamePRO'	=>	'PRO',
		'ncapResNameCYS'	=>	'NCYS',
		'ccapResNameCYS'	=>	'CCYS',
		'lonecapResNameCYS'	=>	'CYS',
		'ncapResNameMET'	=>	'NMET',
		'ccapResNameMET'	=>	'CMET',
		'lonecapResNameMET'	=>	'MET',
	);
	if (exists ($switch{$switchKey})){
		return $switch{$switchKey};
	}
	else {
		die "residue $switchKey not in sub lookupCapName\n";
	}	
}
