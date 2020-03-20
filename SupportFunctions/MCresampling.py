"""
	** MCMC Incremental Slip Rate Calculator **
	Markov chain Monte Carlo sampling schema for Bayesian analysis of fault slip rate data.

	Rob Zinke 2019, 2020
"""

import numpy as np 

def MCMC_resample(Nsamples,Ages,ageList,Dsps,dspList,condition='standard',maxRate=None,seed_value=0,verbose=False,outName=None):	
	"""
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
	"""

	## Setup
	if verbose is True:
		print('Initializing Monte Carlo sampling')

	# Bayesian condition
	if verbose is True:
		print('Condition: {}'.format(condition))

	if condition=='standard':
		condition=lambda ageDiffs,dspDiffs: ageDiffs.min()>=0 and dspDiffs.min()>=0

	if maxRate is None:
		maxRate=np.Inf

	# Random number generator
	np.random.seed(seed_value) # seed random number generator for consistency

	# Arrays to fill in
	m=len(ageList) # number of measurements
	randAges=np.zeros(m)
	randDsps=np.zeros(m)
	AgePicks=np.zeros((m,Nsamples))
	DspPicks=np.zeros((m,Nsamples))
	RatePicks=np.zeros((m-1,Nsamples))


	## Monte Carlo sampling
	tossed=0 # 
	successes=0 # success counter
	while successes<Nsamples:
		# Pick random numbers from uniform distribution
		a_rand=np.random.uniform(0,1,(m)) # random uniform numbers for ages
		d_rand=np.random.uniform(0,1,(m)) # random uniform numbers for disps

		# Interpolate CDF to get random age or disp
		for j in range(m):
			randAges[j]=Ages[ageList[j]].InvCDF(a_rand[j])
			randDsps[j]=Dsps[dspList[j]].InvCDF(d_rand[j])

		# Differences
		ageDiffs=np.diff(randAges)
		dspDiffs=np.diff(randDsps)

		# Check against condition
		if condition(ageDiffs,dspDiffs)==False:
			# If condition not met, try again
			tossed+=1
		else:
			# Check max rate
			rates=dspDiffs/ageDiffs
			if rates.max()>maxRate:
				tossed+=1
			else:
				# If condition is met, record values and advance counter
				AgePicks[:,successes]=randAges
				DspPicks[:,successes]=randDsps
				RatePicks[:,successes]=rates
				successes+=1


	## Finishing
	if verbose is True:
		print('Finished\n\tN successes: {}\n\tN tossed: {}'.format(successes,tossed))

	# Save to file
	if outName:
		# Save picks
		np.save('{}_AgePicks.npy'.format(outName),AgePicks)
		np.save('{}_DspPicks.npy'.format(outName),DspPicks)
		np.save('{}_RatePicks.npy'.format(outName),RatePicks)

	return AgePicks, DspPicks, RatePicks