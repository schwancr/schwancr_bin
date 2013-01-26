#!/usr/bin/env python -u

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-f','--f-committors',dest='fcFN',default='./FCommittors.dat',help='Forward committors or Pfolds of all the states in your MSM')
parser.add_option('-e','--epsilon',dest='eps',type=float,default=0.01,help='The TSE will be defined by  0.5 - eps <= Pfold <= 0.5 + eps')
parser.add_option('-d','--data',dest='dataFNs',action='append',help='Data to calculate the average and stdev for the TSE. You can input multiple files')
parser.add_option('--dc','--data-column',dest='dataCol',default=1,type=int,help='Column to read from the data. This will be the column for ALL input data')

options,args = parser.parse_args()

from pyschwancr.dataIO import readData as rd
from pyschwancr.dataIO import writeData as wd
from numpy import *
import matplotlib
matplotlib.use('Pdf')
from matplotlib.pyplot import *


def main():
	# First read the data:

	fc = rd( options.fcFN )
	
	dataList = []
	for FN in options.dataFNs:
		temp = rd(FN)
		if len(temp.shape) !=1:
			temp = temp[:,options.dataCol]
		dataList.append( temp )

	TSind = where( (0.5-options.eps<=fc)*(fc<=0.5+options.eps) )[0]

	wd( ['TSE',str(options.eps) ], TSind )

	avgList = []
	stdList = []
#	if dataList:
#		for item in dataList:
#			avgList.append( item[TSind].mean() )
#			stdList.append( item[TSind].std() )

#	avgList = array(avgList)
#	stdList = array(stdList)

	nameFN = wd( ['TSEaverages','Names'], array(dataList) )
#	avgFN = wd( ['TSEaverages','Averages','eps-%.2f'%options.eps], avgList )
#	stdFN = wd( ['TSEaverages','StdDevs','eps-%.2f'%options.eps], stdList )

	figure()
	width = 0.5
	#bar( arange(len(avgList))+width/2., avgList,width, yerr=stdList, color='red',alpha=0.5,ecolor='black',edgecolor='black')
	boxplot( [ item[TSind] for item in dataList ] )
	barNames = [ '.'.join( thing.split('_')[-1].split('.')[:-2] ) for thing in options.dataFNs ]
	title('TSE State Averages Pfold=0.5%s %.2f' % ( u"\u00B1", options.eps) )
	ylabel('Fraction Native Contacts')
	xticks( arange(1,len(options.dataFNs)+1), barNames, rotation=60 )
	ylim([0,1])
	savefig('TSEaverages_box_eps-%.2f.pdf'%options.eps)	 


if __name__ == '__main__':
	main()

