#!/usr/bin/env python3
import numpy as np 
import matplotlib.pyplot as plt 

'''
Covert an array of numbers into a probability density function (PDF).
'''

## Histogram method
# Convert to PDF using histogram (more stable over fast intervals)
def arrayHist(V,smoothing=None,verbose=False):
	'''
	INPUTS: 
		V is an array of values
	OUTPUTS:
		x is the values
		px is the probability of occurrence
	'''
	if verbose is True:
		print('Converting array to histogram')

	# Compute raw statistics

	# Normalize area to 1.0

	#return x, px


## KDE method
# Convert to PDF using kernel density estimation (inherently smoother)
def arrayKDE(V,smoothing=None,verbose=False):
	'''
	INPUTS: 
		V is an array of values
	OUTPUTS:
		x is the values
		px is the probability of occurrence
	'''
	if verbose is True:
		print('Converting array to kernel density plot')

	# Compute raw statistics

	# Normalize area to 1.0

	# return x, px