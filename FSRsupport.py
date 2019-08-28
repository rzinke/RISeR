#!/usr/bin/env python3
'''
Objects and support functions for MCMC slip rate calculator
'''

import numpy as np 

class ageDatum:
	'''
	PDF representing an age measurement. Units are thousands of years.
	'''
	def __init__(self,name):
		# Basic parameters
		self.name=name

class dispDatum:
	'''
	PDF representing a displacement measurement. Units are in meters.
	'''
	def __init__(self,name):
		# Basic parameters
		self.name=name



##########################
### Plotting functions ###
##########################


#############################
### Statistical functions ###
#############################

# Highest posterior density

# Quantile
