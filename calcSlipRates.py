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
	parser.add_argument('-verb','--verbose',dest='verbose',type=bool,default=False,help='Verbose? [True/False]')
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
	Ages={}; ageList=[]
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

	# Read disp file
	dspsIn=open(inpt.dsp_list_file,'r')
	dspLines=dspsIn.readlines() # read lines and save to variable
	nDsps=len(dspLines) # how many displacement measurements
	dspsIn.close()

	# Check input files have same number of lines
	assert nAges == nDsps, 'MUST HAVE SAME NUMBER OF AGE AND DISPLACEMENT MEASUREMENTS'
	if inpt.verbose is True:
		print('Ages: {}'.format(ageList))
		#print('Displacements: {}'.format(dspList))
		print('Detected n = {} age and displacement measurements'.format(nAges))