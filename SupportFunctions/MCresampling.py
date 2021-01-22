'''
** RISeR Incremental Slip Rate Calculator **
Markov chain Monte Carlo sampling schema for Bayesian analysis of fault slip rate data.

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import numpy as np
from slipRateObjects import incrSlipRate


### RESAMPLING FUNCTION ---
def MCMCresample(DspAgeData, Nsamples, condition='standard', maxRate=None, bound=None, seedValue=0,
    verbose=False, outName=None):
    '''
    This method uses the inverse transform sampling method (see Zinke et al., 2017; 2019) to
     randomly sample the PDFs of age and displacement measurements provided. A Bayesian
     condition is applied to the outputs to enforce the constraint of no negative slip rates
     (e.g., Gold and Cowgill, 2011; Zinke et al., 2017; 2019). This is herein called the
     'standard' condition. Optionally, different conditions may be specified. A maximum
     slip rate may be specified as well to ensure statistically meaningful values.
    Typically, this function is called by the calcSlipRates wrapper.

    INPUTS
        DspAgeData is a dictionary with one entry per displacement-age datum.
         Each entry has an 'Age' entry (ageDatum) and 'Dsp' entry (dspDatum).
        Nsamples is the number of picks or successful sample runs to complete
        maxRate is the maximum slip rate to be considered. If any picks exceed
         this rate, they will be discarded.
        bound is the upper limit of sampling. Once this limit is exceeded,
         the loop will break no matter what. If unspecified, the default value
         is set to four times the desired number of samples.
    OUTPUTS
        AgePicks is an (m x n) matrix of valid age sample values
        DspPicks is an (m x n) matrix of valid displacement sample values
    '''
    ## Setup
    if verbose == True:
        print('*'*32)
        print('Initializing Monte Carlo sampling')

    Nsamples = int(Nsamples)  # ensure integer value

    # Bayesian condition
    if verbose == True: print('Condition: {}'.format(condition))

    if condition == 'standard':
        condition = lambda ageDiffs, dspDiffs: ageDiffs.min()>=0 and dspDiffs.min()>=0

    # Other conditions
    if maxRate == None:
        maxRate = np.Inf  # max possible rate to be considered

    if bound == None:
        bound = 4*Nsamples  # default is four times the desired number of samples

    # Random number generator
    np.random.seed(seedValue)  # seed random number generator for consistency

    # Arrays to fill in
    m = len(DspAgeData.keys())  # number of measurements
    randAges = np.zeros(m)
    randDsps = np.zeros(m)
    AgePicks = np.zeros((m, Nsamples))
    DspPicks = np.zeros((m, Nsamples))
    RatePicks = np.zeros((m-1, Nsamples))


    ## Monte Carlo sampling
    # Initialize
    tossed = 0  # toss counter
    successes = 0  # success counter
    if verbose == True: print('Progress:')

    # Loop through runs
    while successes < Nsamples:
        # Pick random numbers from uniform distribution
        a_rand=np.random.uniform(0, 1, (m))  # random uniform numbers for ages
        d_rand=np.random.uniform(0, 1, (m))  # random uniform numbers for disps

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

        # Report progress
        if verbose == True:
            if 100*successes/Nsamples % 10 == 0:
                print('{:.0f} %'.format(100*successes/Nsamples))

        # Check against bound
        if (successes+tossed) > bound:
            print('WARNING! Maximum sample limit exceeded ({:d}).'.format(bound))
            break


    ## Finishing
    if verbose == True:
        print('Finished')
        print('\tN successes: {:d}'.format(successes))
        print('\tN tossed: {:d}'.format(tossed))

    # Save picks to file
    if outName:
        savename = '{:s}_Picks'.format(outName)
        np.savez(savename,
            AgePicks=AgePicks,
            DspPicks=DspPicks,
            RatePicks=RatePicks)

    return AgePicks, DspPicks, RatePicks



### CONVERT PICKS TO PDF ---
def picks2PDF(DspAgeData, RatePicks, method, stepSize, smoothingKernel=None, kernelWidth=2, verbose=False):
    '''
    Convert rate picks to a probability function using the specified
     method.
    '''
    # Parameters
    dataNames = list(DspAgeData.keys())
    m = RatePicks.shape[0]

    if verbose == True:
        print('*'*32)
        print('Converting slip rate picks to PDFs using {:s} method'.format(method))

    intervals = []
    Rates = {}  # dictionary of object instances
    for i in range(m):
        # Formulate interval name
        intvl = '{:s}-{:s}'.format(dataNames[i], dataNames[i+1])
        intervals.append(intvl)

        # Create slip rate object
        Rates[intvl] = incrSlipRate(name=intvl)

        # Convert picks to PDF
        Rates[intvl].picks2PDF(RatePicks[i,:], method, stepSize, smoothingKernel, kernelWidth)

    return Rates



### STATISTICAL FUNCTIONS ---
def rawPercentiles(DspAgeData, RatePicks, txtFile, confidence=68.27, verbose=False):
    '''
    Compute the percentiles of slip rate picks.
    '''
    # Parameters
    dataNames = list(DspAgeData.keys())
    m, Nsamples = RatePicks.shape  # nb displacement markers
    percentiles=[50-confidence/2, 50, 50+confidence/2]

    # Report if requested
    if verbose == True:
        print('Raw slip rate picks: (median values and {:.2f}% confidence)'.format(confidence))

    # Write to text file
    txtFile.append('\nIncremental slip rates based on percentiles of {:d} slip rate picks ({:.2f}% confidence):\n'.\
        format(Nsamples, confidence))

    # Compute 
    for i in range(m):
        # Formulate interval name
        intvl = '{:s}-{:s}'.format(dataNames[i], dataNames[i+1])

        # Calculate percentile
        pct = np.percentile(RatePicks[i,:], percentiles)
        median = pct[1]
        high_err = pct[2]-pct[1]
        low_err = pct[1]-pct[0]
        rawStats = '{0:s}: {1:.2f} +{2:.2f} -{3:.2f}'.format(intvl, median, high_err, low_err)
        
        # Report if requested
        if verbose == True: print(rawStats)

        # Update text file
        txtFile.append(rawStats+'\n')