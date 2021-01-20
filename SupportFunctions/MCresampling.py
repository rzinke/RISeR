'''
** MCMC Incremental Slip Rate Calculator **
Markov chain Monte Carlo sampling schema for Bayesian analysis of fault slip rate data.

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import numpy as np


### RESAMPLING FUNCTION ---
def MCMCresample(DspAgeData, Nsamples, condition='standard', maxRate=None, seed_value=0, verbose=False, outName=None):
    '''
    This method uses the inverse transform sampling method (see Zinke et al., 2017; 2019) to
     randomly sample the PDFs of age and displacement measurements provided. A Bayesian
     condition is applied to the outputs to enforce the constraint of no negative slip rates
     (e.g., Gold and Cowgill, 2011; Zinke et al., 2017; 2019). This is herein called the
     "standard" condition. Optionally, different conditions may be specified. A maximum
     slip rate may be specified as well to ensure statistically meaningful values.
    Typically, this function is called by the calcSlipRates wrapper.

    INPUTS
        DspAgeData is a dictionary with one entry per displacement-age datum.
         Each entry has an "Age" entry (ageDatum) and "Dsp" entry (dspDatum).
        Nsamples is the number of picks or successful sample runs to complete
    OUTPUTS
        AgePicks is an (m x n) matrix of valid age sample values
        DspPicks is an (m x n) matrix of valid displacement sample values
    '''

    ## Setup
    if verbose == True:
        print('*******************************')
        print('Initializing Monte Carlo sampling')

    # Bayesian condition
    if verbose == True:
        print('Condition: {}'.format(condition))

    if condition == 'standard':
        condition=lambda ageDiffs, dspDiffs: ageDiffs.min()>=0 and dspDiffs.min()>=0

    if maxRate == None:
        maxRate=np.Inf

    # Random number generator
    np.random.seed(seed_value)  # seed random number generator for consistency

    # Arrays to fill in
    m = len(DspAgeData.keys())  # number of measurements
    randAges = np.zeros(m)
    randDsps = np.zeros(m)
    AgePicks = np.zeros((m, Nsamples))
    DspPicks = np.zeros((m, Nsamples))
    RatePicks = np.zeros((m-1, Nsamples))


    ## Monte Carlo sampling
    tossed=0  # toss counter
    successes=0  # success counter
    while successes < Nsamples:
        # Pick random numbers from uniform distribution
        a_rand=np.random.uniform(0, 1, (m)) # random uniform numbers for ages
        d_rand=np.random.uniform(0, 1, (m)) # random uniform numbers for disps

        # Interpolate CDF to get random age or disp
        for j, datumName in enumerate(DspAgeData.keys()):
            # Age and displacement objects for each datum
            Age = DspAgeData[datumName]['Age']
            Dsp = DspAgeData[datumName]['Dsp']
            # Random samples
            randAges[j] = Age.InvCDF(a_rand[j])
            randDsps[j] = Dsp.InvCDF(d_rand[j])

        # Differences
        ageDiffs = np.diff(randAges)
        dspDiffs = np.diff(randDsps)

        # Check against condition
        if condition(ageDiffs, dspDiffs) == False:
            # If condition not met, try again
            tossed += 1
        else:
            # Check max rate
            rates = dspDiffs/ageDiffs
            if rates.max() > maxRate:
                tossed += 1
            else:
                # If condition is met, record values and advance counter
                AgePicks[:,successes] = randAges
                DspPicks[:,successes] = randDsps
                RatePicks[:,successes] = rates
                successes += 1


    ## Finishing
    if verbose == True:
        print('Finished\n\tN successes: {}\n\tN tossed: {}'.format(successes, tossed))

    # Save to file
    if outName:
        # Save picks
        savename = '{}_Picks'.format(outName)
        np.savez(savename,
            AgePicks=AgePicks,
            DspPicks=DspPicks,
            RatePicks=RatePicks)

    return AgePicks, DspPicks, RatePicks