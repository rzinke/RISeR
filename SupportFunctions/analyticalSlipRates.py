#!/usr/bin/env python3
'''
** RISeR Incremental Slip Rate Calculator **
This function will compute the incremental slip rates of a data
 set consisting of offsets and age markers described as
 probability density functions (PDFs).

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import numpy as np
from slipRateObjects import incrSlipRate
from differencePDFs import PDFdiff
from dividePDFs import PDFquotient



### SLIP RATE FUNCTIONS ---
def analyticalSlipRates(DspAgeData, stepSize=None, maxRate=None, verbose=False, printDetails=False):
    '''
    For each pair of dated displacement markers, compute the incremental slip
     rate.
    First, calculate the difference between pairs of ages and displacements.
    Then, compute the quotient of each pair. Limit the quotient computations
     only to values greater than zero to ensure no negative slip rates.

    INPUTS
        DspAgeData is the dictionary of dated displacement markers loaded using
         the loadDspAgeInputs function in the dataLoading module
        stepSize is the sample size of the quotient axis
        maxRate is the maximum rate to be considered
    '''
    ## SETUP
    # Parameters
    markerNames = list(DspAgeData.keys())
    m = len(markerNames)

    # Report if requested
    if verbose == True:
        print('*'*32)
        print('Computing slip rates using the analytical formulation')


    ## COMPUTATIONS
    Rates = {}  # empty dictionary for storing incremental slip rates

    # Loop through each pair of markers
    for i in range(m-1):
        # Data names
        youngerName = markerNames[i]
        olderName = markerNames[i+1]
        intvl = '{:s}-{:s}'.format(markerNames[i], markerNames[i+1])

        if verbose == True:
            print('*'*32)
            print('Interval: {:s} '.format(intvl))

        # Compute age difference
        if verbose == True: print('Computing age difference')

        youngerAge = DspAgeData[youngerName]['Age']  # younger PDF
        olderAge = DspAgeData[olderName]['Age']  # older PDF

        deltaAge = PDFdiff(olderAge.ages, olderAge.probs, youngerAge.ages, youngerAge.probs, verbose=printDetails)

        # Compute displacement difference
        smallerDsp = DspAgeData[youngerName]['Dsp']  # smaller PDF
        largerDsp = DspAgeData[olderName]['Dsp']  # larger PDF

        deltaDsp = PDFdiff(largerDsp.dsps, largerDsp.probs, smallerDsp.dsps, smallerDsp.probs, verbose=printDetails)

        # Remove points with zero probability
        ageNdx = (deltaAge.D > 0)  # indices with probability greater than zero
        deltaAge.D = deltaAge.D[ageNdx]
        deltaAge.pD = deltaAge.pD[ageNdx]

        dspNdx = (deltaDsp.D > 0)  # indices with probability greater than zero
        deltaDsp.D = deltaDsp.D[dspNdx]
        deltaDsp.pD = deltaDsp.pD[dspNdx]

        # Compute incremental slip rates
        slipRate = PDFquotient(deltaDsp.D, deltaDsp.pD, deltaAge.D, deltaAge.pD, stepSize=stepSize, Qmax=maxRate,
            verbose=printDetails)

        # Record to dictionary
        Rates[intvl] = incrSlipRate(intvl)
        Rates[intvl].rates = slipRate.Q
        Rates[intvl].probs = slipRate.pQ

    return Rates