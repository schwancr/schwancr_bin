#!/usr/bin/perl -w
use strict;

# Input a .pdb file format and output a .txt file formatted for input
# into GAMESS using the COORD=UNIQUE format for atom input

my $pdbIn = $ARGV[0];
my $txtOut = $ARGV[1];
my %atmChg = 
    (
    H => 1.0,
    C => 6.0,
    N => 7.0,
    O => 8.0,
    P => 15.0,
    );

open (my $input, "<$pdbIn") or die "cannot open file $pdbIn\n";
open (my $output, ">$txtOut") or die "cannot open file $txtOut\n";

# Files are open, now search for lines in input that begin with ATOM or HETATM

while (my $line = <$input>)
{
    if ($line =~ m/\AATOM/ || $line =~ m/\AHETATM/)
    {
	 # Get atom name
	 my $atmName = substr($line,13,1);
	 my $tmpCoord = substr($line,30,24);

	 
	 print $output "$atmName $atmChg{$atmName} $tmpCoord\n";
    }
}

close ($input);
close ($output);
