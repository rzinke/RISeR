#!/usr/bin/env python3
"""
	** MCMC Incremental Slip Rate Calculator **
	Quickly plot displacement and age data before calculating 
	 slip rates by MC sampling.

	Rob Zinke 2019, 2020
"""

### IMPORT MODULES ---
import os
import numpy as np 
import matplotlib.pyplot as plt 
from slipRateObjects import *
from calcSlipRates import plotRawData, plotMCresults


### PARSER ---
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Quickly plot displacement and age data before calculating slip rates by MC sampling.')
	# Required
	parser.add_argument('-a','--age-list', dest='ageListFile', type=str, required=True, help='Text file with one age file per line, list in order from youngest (top line) to oldest (bottom).')
	parser.add_argument('-d','--dsp-list', dest='dspListFile', type=str, required=True, help='Text file with one displacement file per line, list in order from smallest (top line) to largest (bottom).')
	parser.add_argument('-o','--out-name', dest='outName', type=str, default='Out', help='Head name for outputs (no extension)')
	# Recommended
	parser.add_argument('-pt','--plot-type', dest='plotType', type=str, default='whisker', help='Plot marker type [\'wkisker\',\'rectangle\']')
	parser.add_argument('-t','--title', dest='title', type=str, default=None, help='Plot title')
	parser.add_argument('-l','--labels', dest='labels', action='store_true', default=False, help='Label features')
	parser.add_argument('-verb','--verbose', dest='verbose', action='store_true', default=False, help='Verbose?')
	parser.add_argument('--plot-inputs',dest='plotInputs', action='store_true', help='Plot inputs')
	parser.add_argument('--plot-outputs',dest='plotOutputs', action='store_true', help='Plot outputs')
	return parser 

def cmdParser(inpt_args=None):
	parser = createParser()
	return parser.parse_args(inpt_args)



### MAIN ---
if __name__ == '__main__':
	inpt=cmdParser()

	## Load files
	# Read age file
	agesIn=open(inpt.ageListFile,'r')
	ageLines=agesIn.readlines() # read lines and save to variable
	nAges=len(ageLines) # how many age measurements
	agesIn.close()

	# Load ages into dictionary
	Ages={}; ageList=[]; xmaxGlobal=0.
	for ageLine in ageLines:
		"""
			Record one dictionary entry per age file. Each entry in an 
			 instance of an "ageDatum" object. Append age name (path and 
			 suffix removed) to list for later.
		"""
		# Isolate basename
		ageLine=ageLine.strip('\n') # remove extraneous newline
		ageName=os.path.basename(ageLine)
		ageName=ageName.split('.')[0] # remove file extension
		ageList.append(ageName) # append to list

		# Spawn age object and record to dictionary
		Ages[ageName]=ageDatum(ageName)
		Ages[ageName].readFromFile(ageLine) # load data to object

		# Format for slip rate analysis
		#  builds inverse interpolation function
		Ages[ageName].format(verbose=inpt.verbose,plot=inpt.plotInputs)

		# Set global limit for plotting
		if Ages[ageName].upperLimit>xmaxGlobal:
			xmaxGlobal=Ages[ageName].upperLimit


	# Read disp file
	dspsIn=open(inpt.dspListFile,'r')
	dspLines=dspsIn.readlines() # read lines and save to variable
	nDsps=len(dspLines) # how many displacement measurements
	dspsIn.close()

	# Load displacements into dictionary
	Dsps={}; dspList=[]; ymaxGlobal=0.
	for dspLine in dspLines:
		"""
			Record one dictionary entry per displacement file. Each entry in an instance of 
			an "dspDatum" object. Append displacement name (path and suffix removed) to list for later.
		"""
		# Isolate basename
		dspLine=dspLine.strip('\n') # remove extraneous newline
		dspName=os.path.basename(dspLine)
		dspName=dspName.split('.')[0] # remove file extension
		dspList.append(dspName) # append to list

		# Spawn age object and record to dictionary
		Dsps[dspName]=dspDatum(dspName)
		Dsps[dspName].readFromFile(dspLine) # load data to object

		# Format for slip rate analysis
		#	builds inverse interpolation function
		Dsps[dspName].format(verbose=inpt.verbose,plot=inpt.plotInputs)

		# Set global limit for plotting
		if Dsps[dspName].upperLimit>ymaxGlobal:
			ymaxGlobal=Dsps[dspName].upperLimit

	# Check input files have same number of lines
	assert nAges == nDsps, 'MUST HAVE SAME NUMBER OF AGE AND DISPLACEMENT MEASUREMENTS'
	m=nAges # assign number of measurements
	if inpt.verbose is True:
		print('Detected m = {} age and displacement measurements'.format(nAges))
		# Confirm pairing
		print('Pairing (youngest at top):')
		for i in range(nAges):
			print('\t{} - {}'.format(ageList[i],dspList[i]))

	
	## Formulate interval names
	#  Intervals are the rates between one measurement and another
	intervalList=[]
	for i in range(m-1):
		intervalName='{}-{}'.format(dspList[i],dspList[i+1])
		intervalList.append(intervalName)


	## Plot raw data (whisker plot)
	if inpt.plotType.lower() in ['whisker','whiskers']:
		# Whisker plot
		Fraw,axRaw=plotRawData(Ages,ageList,Dsps,dspList,
			xmaxGlobal,ymaxGlobal)
	elif inpt.plotType.lower() in ['rectangle','rectangles','box']:
		# Rectangle plot
		Fraw,axRaw=plotMCresults(Ages,ageList,Dsps,dspList,
			AgePicks=-np.ones((1,1)),DspPicks=-np.ones((1,1)),
			xMax=1.1*xmaxGlobal,yMax=1.1*ymaxGlobal,maxPicks=0)

	# Label if desired
	if inpt.labels is True:
		for i in range(m):
			x_mid=Ages[ageList[i]].median
			x_err=np.array([[x_mid-Ages[ageList[i]].lowerLimit],
				[Ages[ageList[i]].upperLimit-x_mid]])
			y_mid=Dsps[dspList[i]].median
			y_err=np.array([[y_mid-Dsps[dspList[i]].lowerLimit],
				[Dsps[dspList[i]].upperLimit-y_mid]])
			feature_name=Dsps[dspList[i]].name
			axRaw.text(1.02*x_mid,1.02*y_mid,feature_name)

	# Format chart
	axRaw.set_xlim([0,1.1*xmaxGlobal]) # x-limits
	axRaw.set_ylim([0,1.1*ymaxGlobal]) # y-limits
	axRaw.set_xlabel('age'); axRaw.set_ylabel('displacement')

	# Title if specified
	if inpt.title:
		axRaw.set_title('{} raw data (95 % limits)'.format(inpt.title))
	else:
		axRaw.set_title('Raw data (95 % limits)')


	## Outputs
	# Save output if specified 
	if inpt.outName:
		if inpt.labels is True:
			Fraw.savefig('{}_RawData-labelled.png'.format(inpt.outName),dpi=600)
		else:
			Fraw.savefig('{}_RawData.png'.format(inpt.outName),dpi=600)

	# Show plots if specified
	if inpt.plotInputs or inpt.plotOutputs:
		plt.show()