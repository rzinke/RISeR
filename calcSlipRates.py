#!/usr/bin/env python3
'''
'''
import os
import numpy as np
import matplotlib.pyplot as plt
from SlipRateObjects import *



## Parser
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Main function for calculating incremental slip rates.')
	parser.add_argument('-a','--age_list',dest='age_list_file',type=str,required=True,help='Text file with one age file per line, list in order from youngest (top line) to oldest (bottom).')
	parser.add_argument('-d','--dsp_list',dest='dsp_list_file',type=str,required=True,help='Text file with one displacement file per line, list in order from smallest (top line) to largest (bottom).')
	parser.add_argument('-o','--output',dest='outName',type=bool,default=False,help='Head name for outputs')
	parser.add_argument('-verb','--verbose',dest='verbose',type=bool,default=False,help='Verbose? [True/False]')
	parser.add_argument('-plot_inputs','--plot_inputs',dest='plot_inputs',type=bool,default=False,help='Plot inputs')
	parser.add_argument('-plot_outputs','--plot_outputs',dest='plot_outputs',type=bool,default=False,help='Plot outputs')
	return parser 

def cmdParser(inpt_args=None):
	parser = createParser()
	return parser.parse_args(inpt_args)



## Main
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

	## Plot raw data
	F=plt.figure('RawData')
	ax=F.add_subplot(111)
	for i in range(m):
		x_mid=Ages[ageList[i]].median
		x_err=np.array([[x_mid-Ages[ageList[i]].lowerLimit],
			[Ages[ageList[i]].upperLimit-x_mid]])
		y_mid=Dsps[dspList[i]].median
		y_err=np.array([[y_mid-Dsps[dspList[i]].lowerLimit],
			[Dsps[dspList[i]].upperLimit-y_mid]])
		ax.errorbar(x_mid,y_mid,xerr=x_err,yerr=y_err,
			color=(0.3,0.3,0.6),marker='o')
	ax.set_xlim([0,1.1*xmaxGlobal]) # x-limits
	ax.set_ylim([0,1.1*ymaxGlobal]) # y-limits
	ax.set_xlabel('age'); ax.set_ylabel('offset')
	ax.set_title('Raw data')
	F.savefig('{}_RawData.png'.format(inpt.outName),dpi=300)


	if inpt.plot_inputs==True or inpt.plot_outputs==True:
		plt.show()