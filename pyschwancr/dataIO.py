#!/usr/bin/env python

from numpy import *
from msmbuilder import Trajectory, io
from matplotlib.backends import backend_pdf
from matplotlib.pyplot import *
import os, re, sys

def readData(FN):
	"""
This function will read data from a filename based on it's extension.
Inputs:
	1) FN: filename to find data
Outputs:
	2) data: numpy array containing data read in.

The function tries load and loadtxt from numpy, and throws an error if it cannot read the file.
Additionally if the filename is lh5 it will load the data as a Trajectory.Trajectory Object but return the array of coordinates.

	"""
	if FN.split('.')[-1] == 'npy':
		data = load(FN)
	elif FN.split('.')[-1] == 'lh5': # Assume this is a trajectory in msmbuilder
		data = Trajectory.load_from_lhdf( FN )['XYZList'] # Only return the coordinates.
	elif FN.split('.')[-1] == 'h5':
		data = io.loadh( FN )
        try: data = data['arr_0']
        except: data = data['Data']
		data = data[ where( data != -1 ) ]
	else:
		try:
			data = loadtxt(FN)
		except:
			print "\n\n dataIO.readData: Cannot read %s. Use numpy.save or numpy.savetxt. Exiting" % FN
			exit()
	if data.shape == (): # If there is only one data point, then turn it into a list
		data = array( [ data ] )
	return data

def writeData(inputs, data, txt=True, dir='./'):
	"""
This function will write a data using the input file names in the list inputs.
Inputs:
	1) inputs: List of things to place in the output Filename
	2) data: (numpy array)
	3) (optional) txt: if True (default) this function will use savetxt, otherwise it will use np.save
	4) (optional) dir: The directory to write output file to. Default = './'
Output:
	1) outName: Filename saved with data

The format of the output will be:
	input[0]_input[1]_...input[-1].ext
	"""
	# Remove the directories from any files:
	temp = [ thing.split('/')[-1] for thing in inputs ]
	# Remove extensions (if they're there)
	temp = [ '.'.join( thing.split('.')[:-1]) if len( thing.split('.') ) > 1 and len(thing.split('.')[-1]) == 3 else thing for thing in temp ]
	# Next join the list:	
	outName = '_'.join( temp )
	# Add the directory:
	outName = os.path.join(dir,outName)
	
	if txt:
		outName += '.dat'
		savetxt(outName, data)
	else:
		outName += '.npy'
		save(outName, data)
	
	return outName

def plotData( xDat, yDat, x_lbl, y_lbl, LinePlot = False, LabelsAreFiles = False, outName = None ):

	if LinePlot:
	   plot( xDat, yDat )
	else:
		plot( xDat, yDat, '.' )

	if LabelsAreFiles:
		x_lbl = '.'.join( x_lbl.split('/')[-1].split('.')[:-1] )
		y_lbl = '.'.join( y_lbl.split('/')[-1].split('.')[:-1] )

	title("%s vs. %s" % ( y_lbl, x_lbl) )
	ylabel(y_lbl)
	xlabel(x_lbl)

	xRange = xDat.max() - xDat.min()
	x0 = xDat.min() - 0.01 * xRange
	x1 = xDat.max() + 0.01 * xRange
	xlim([ x0, x1 ])

	yRange = yDat.max() - yDat.min()
	y0 = yDat.min() - 0.01 * yRange
	y1 = yDat.max() + 0.01 * yRange
	ylim([ y0, y1 ])
	
	if outName == None:
		savefig( "%s_vs_%s.pdf" % ( x_lbl, y_lbl ) )
	else:
		savefig( outName )
	return

def getTrajList( dir, BeginInd=3, EndInd=-4, RegEx=r'^trj\d+\.lh5' ):
	"""
	This function creates a list of data files indexed by integers, and named in some regular way. T
	The files should be named in such a way that contains an integer. 
	
	Inputs:
	1) dir - Directory containing the files
	2) BeginInd [ 3 ] - Index which is the first integer digit in the filename
	3) EndInd [ -4 ] - Index which is one past the last integer digit in the filename
	4) RegEx [ r\'^trj\\d+\\.lh5\' ] - Regular expression to use to filter things to place in the list

	Outputs:
	1) trajList - List of filenames INCLUDING THE DIRECTORY (i.e. os.path.join( dir, fn ) )

	EG)
	
	Trajectories:	
	trj0.lh5  trj1.lh5  trj2.lh5  trj3.lh5 ...

	getTrajList( \'Trajectories\', BeginInd=3, EndInd=-4, RegEx=r\'trj\\d+\\.lh5\' )

	[ Trajectories/trj0.lh5, Trajectories/trj1.lh5, ... ]

	Since fn[3] = 0, fn[-4] = \'.\'

	"""
	trajList = [ ( int( fn[ BeginInd : EndInd ] ), fn ) for fn in os.listdir(dir) if re.search( RegEx, fn ) ]
	trajList.sort()
	trajList = [ os.path.join( dir, fn ) for ( i, fn ) in trajList ]

	return trajList
