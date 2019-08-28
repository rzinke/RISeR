#!/usr/bin/env python3
'''
'''
import numpy as np
import matplotlib.pyplot as plt
from FSRsupport import *



## Parser
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Main function for calculating incremental slip rates.')
	parser.add_argument('-a','--age_list',dest='age_list',type=str,required=True,help='Text file with one age file per line, list in order from youngest (top line) to oldest (bottom).')
	return parser 

def cmdParser(inpt_args=None):
	parser = createParser()
	return parser.parse_args(inpt_args)



## Read in age data



## Main
if __name__ == '__main__':
	inpt=cmdParser()

	## Load files
	# Check input files have same number of lines
	nLines_age_file=ageList.readlines()
	nLines_dsp_file=inpt.dsp_list