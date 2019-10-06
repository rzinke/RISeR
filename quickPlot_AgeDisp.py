#!/usr/bin/env python3
import os
import numpy as np 
import matplotlib.pyplot as plt 
from slipRateObjects import *
from calcSlipRates import plotRawData, plotMCresults

### --- Parser ---
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Quickly plot displacement and age data before calculating slip rates by MC sampling.')
	# Required
	parser.add_argument('-a','--age_list',dest='age_list_file',type=str,required=True,help='Text file with one age file per line, list in order from youngest (top line) to oldest (bottom).')
	parser.add_argument('-d','--dsp_list',dest='dsp_list_file',type=str,required=True,help='Text file with one displacement file per line, list in order from smallest (top line) to largest (bottom).')
	parser.add_argument('-o','--output',dest='outName',type=str,default='Out',help='Head name for outputs (no extension)')
	# Recommended
	parser.add_argument('-plot_type','--plot_type',dest='plot_type',type=str,default='whisker',help='Plot marker type [\'wkisker\',\'rectangle\']')
	parser.add_argument('-t','--title',dest='title',type=str,default=None,help='Plot title')
	parser.add_argument('-l','--labels',dest='labels',action='store_true',default=False,help='Label features')
	parser.add_argument('-verb','--verbose',dest='verbose',action='store_true',default=False,help='Verbose?')
	parser.add_argument('-plot_inputs','--plot_inputs',dest='plot_inputs',action='store_true',help='Plot inputs')
	parser.add_argument('-plot_outputs','--plot_outputs',dest='plot_outputs',action='store_true',help='Plot outputs')
	return parser 


def cmdParser(inpt_args=None):
	parser = createParser()
	return parser.parse_args(inpt_args)



### --- Main function ---
if __name__ == '__main__':
	inpt=cmdParser()

	## Load files
	# Read age file
	agesIn=open(inpt.age_list_file,'r')
	ageLines=agesIn.readlines() # read lines and save to variable
	nAges=len(ageLines) # how many age measurements
	agesIn.close()

	# Load ages into dictionary
	Ages={}; ageList=[]; xmaxGlobal=0.
	for age_line in ageLines:
		'''
		Record one dictionary entry per age file. Each entry in an instance of 
		an "ageDatum" object. Append age name (path and suffix removed) to list for later.
		'''
		# Isolate basename
		age_line=age_line.strip('\n') # remove extraneous newline
		age_name=os.path.basename(age_line)
		age_name=age_name.split('.')[0] # remove file extension
		ageList.append(age_name) # append to list

		# Spawn age object and record to dictionary
		Ages[age_name]=ageDatum(age_name)
		Ages[age_name].readFromFile(age_line) # load data to object

		# Format for slip rate analysis
		#	builds inverse interpolation function
		Ages[age_name].format(verbose=inpt.verbose,plot=inpt.plot_inputs)

		# Set global limit for plotting
		if Ages[age_name].upperLimit>xmaxGlobal:
			xmaxGlobal=Ages[age_name].upperLimit


	# Read disp file
	dspsIn=open(inpt.dsp_list_file,'r')
	dspLines=dspsIn.readlines() # read lines and save to variable
	nDsps=len(dspLines) # how many displacement measurements
	dspsIn.close()

	# Load displacements into dictionary
	Dsps={}; dspList=[]; ymaxGlobal=0.
	for dsp_line in dspLines:
		'''
		Record one dictionary entry per displacement file. Each entry in an instance of 
		an "dspDatum" object. Append displacement name (path and suffix removed) to list for later.
		'''
		# Isolate basename
		dsp_line=dsp_line.strip('\n') # remove extraneous newline
		dsp_name=os.path.basename(dsp_line)
		dsp_name=dsp_name.split('.')[0] # remove file extension
		dspList.append(dsp_name) # append to list

		# Spawn age object and record to dictionary
		Dsps[dsp_name]=dspDatum(dsp_name)
		Dsps[dsp_name].readFromFile(dsp_line) # load data to object

		# Format for slip rate analysis
		#	builds inverse interpolation function
		Dsps[dsp_name].format(verbose=inpt.verbose,plot=inpt.plot_inputs)

		# Set global limit for plotting
		if Dsps[dsp_name].upperLimit>ymaxGlobal:
			ymaxGlobal=Dsps[dsp_name].upperLimit

	# Check input files have same number of lines
	assert nAges == nDsps, 'MUST HAVE SAME NUMBER OF AGE AND DISPLACEMENT MEASUREMENTS'
	m=nAges # assign number of measurements
	if inpt.verbose is True:
		print('Detected m = {} age and displacement measurements'.format(nAges))
		# Confirm pairing
		print('Pairing (youngest at top):')
		for i in range(nAges):
			print('\t{} - {}'.format(ageList[i],dspList[i]))

	# Formulate interval names
	#  Intervals are the rates between one measurement and another
	intervalList=[]
	for i in range(m-1):
		interval_name='{}-{}'.format(dspList[i],dspList[i+1])
		intervalList.append(interval_name)

	## Plot raw data (whisker plot)
	if inpt.plot_type.lower() in ['whisker','whiskers']:
		# Whisker plot
		Fraw,axRaw=plotRawData(Ages,ageList,Dsps,dspList,
			xmaxGlobal,ymaxGlobal)
	elif inpt.plot_type.lower() in ['rectangle','rectangles','box']:
		# Rectangle plot
		Fraw,axRaw=plotMCresults(Ages,ageList,Dsps,dspList,
			AgePicks=-np.ones((1,1)),DspPicks=-np.ones((1,1)),
			xMax=1.1*xmaxGlobal,yMax=1.1*ymaxGlobal,max_picks=0)

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
	axRaw.set_xlabel('age'); axRaw.set_ylabel('offset')
	# Title if specified
	if inpt.title:
		axRaw.set_title('{} raw data (95 % limits)'.format(inpt.title))
	else:
		axRaw.set_title('Raw data (95 % limits)')
	# Save output if specified 
	if inpt.outName:
		if inpt.labels is True:
			Fraw.savefig('{}_Fig1_RawData-labelled.png'.format(inpt.outName),dpi=600)
		else:
			Fraw.savefig('{}_Fig1_RawData.png'.format(inpt.outName),dpi=600)

	# Show plots if specified
	if inpt.plot_inputs or inpt.plot_outputs:
		plt.show()