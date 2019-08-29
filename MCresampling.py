#!/usr/bin/env python3
'''
Markov chain Monte Carlo sampling schema for Bayesian analysis of fault slip rate data.
'''
import numpy as np 

def MCMC_resample(Nsamples,Ages,ageList,Dsps,dspList,condition='standard',maxRate=None,seed_value=0,outName=None):	
	'''
	This method uses the inverse transform sampling method (see Zinke et al., 2017; 2019) to 
	 randomly sample the PDFs of age and displacement measurements provided. A Bayesian 
	 condition is applied to the outputs to enforce the constraint of no negative slip rates
	 (e.g., Gold and Cowgill, 2011; Zinke et al., 2017; 2019). This is herein called the 
	 "standard" condition. Optionally, different conditions may be specified. A maximum
	 slip rate may be specified as well to ensure statistically meaningful values.
	Typically, this function is called by the calcSlipRates wrapper.

	INPUTS
		Nsamples is the number of picks or successful sample runs to complete
		Ages is a dictionary of ageDatum objects (see calcSlipRates and slipRateObjects)
		ageList is a list of age measurement names ordered youngest-oldest (see calcSlipRates)
		Dsps is a dictionary like Ages but for displacement measurements
		dspList is a list like ageList with displacement names ordered smallest-largest
	OUTPUTS
		AgePicks is an (m x n) matrix of valid age sample values
		DspPicks is an (m x n) matrix of valid displacement sample values
	'''

	# Setup
	np.random.seed(seed_value) # seed random number generator for consistency

	m=len(ageList) # number of measurements
	r_ages=np.zeros(m)
	r_dsps=np.zeros(m)
	AgePicks=np.zeros((m,Nsamples))
	DspPicks=np.zeros((m,Nsamples))

	# Monte Carlo runs
	i=0 # initialize counter
	while i<Nsamples:
		# Pick random numbers from uniform distribution
		a_rand=np.random.uniform(0,1,(m)) # random uniform numbers for ages
		d_rand=np.random.uniform(0,1,(m)) # random uniform numbers for disps

		# Interpolate CDF to get random age or disp
		for j in range(m):
			r_ages[j]=Ages[ageList[j]].InvCDF(a_rand[j])
			r_dsps[j]=Dsps[dspList[j]].InvCDF(d_rand[j])

		# Condition met
		AgePicks[:,i]=r_ages
		DspPicks[:,i]=r_dsps

		i+=1


	# Save to file
	if outName:
		pass

	return AgePicks, DspPicks#, RatePicks