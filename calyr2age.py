#!/usr/bin/env python3
'''
** RISeR Incremental Slip Rate Calculator **
Convert calendar years to years before present

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from PDFanalysis import smoothPDF
from resultSaving import confirmOutputDir


### PARSER ---
Description = '''Convert a date PDF in CE/BCE to years before present (ybp). This is designed to
work with OxCal posterior distributions in which the left column is age in CE/BCE and the right
column is relative probability.'''

Examples = '''EXAMPLES

# Convert date to age relative to 1950 C.E.
calyr2age.py Date1_posterior.txt -o Age1

# Convert date to age relative to 2020 C.E.
calyr2age.py Date1_posterior.txt -o Age1 -r 2020

# Convert date to age, smooth age and convert to thousands years B.P.
calyr2age.py Date1_posterior.txt -o Age1 -s 3 -f 1000
'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    parser.add_argument(dest='datefile',
        help='File containing the date in OxCal posterior format.')
    parser.add_argument('-o','--outName', dest='outName', required = True,
        type=str,help='Output file name.')
    parser.add_argument('-r','--reference-date', dest='refDate', default=1950, type=float,
        help='Reference calendar date (C.E.). Default = 1950.')
    parser.add_argument('-f','--age-factor', dest='ageFactor', default=1, type=float,
        help='Factor to divide age by (e.g., 1000). Default = 1.')
    parser.add_argument('-s','--smoothing', dest='smoothing', default=None, type=int,
        help='Width of smoothing kernel.')
    parser.add_argument('-k','--smoothing-kernel', dest='smoothingKernel', default='gaussian', type=str,
        help='Smoothing kernel type [boxcar/Gaussian]')
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode.')
    parser.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot outputs.')
    return parser

def cmdParser(inpt_args = None):
    parser = createParser()
    return parser.parse_args(inpt_args)



### ANCILLARY FUNCTIONS ---
def loadDates(datefile, verbose=False):
    '''
    Load data from two column file in the form of OxCal posterior.
    '''
    datePDF = np.loadtxt(datefile)
    xDate = datePDF[:,0]  # dates
    pxDate = datePDF[:,1]  # probabilities

    if verbose == True: print('Loaded: {:s}'.format(datefile))

    return xDate, pxDate


def date2age(xDate, pxDate, refDate=1950, ageFactor=1, verbose=False):
    '''
    Convert calendar date to age before specified date. Scale by age factor.
    '''

    # Shift relative to age
    xAge = refDate - xDate

    # Flip arrays for chronological order
    xAge = np.flip(xAge)
    pxAge = np.flip(pxDate)

    # Scale by age factor
    xAge /= ageFactor

    # Normalize area to 1.0
    P = np.trapz(pxAge, xAge)
    pxAge /= P

    if verbose == True:
        print('Converted to age relative to reference date: {}'.format(refDate))
        if ageFactor != 1: print('Scaled by age factor: {}'.format(ageFactor))
        print('Normalized area to unit mass')

    return xAge, pxAge


def saveOutputs(x, px, outName, verbose=False):
    '''
    Save outputs to text file.
    '''

    fname = outName+'.txt'
    with open(fname, 'w') as outFile:
        outFile.write('# Value,\tProbability\n')
        for i in range(len(x)):
            outFile.write('{0:f}\t{1:f}\n'.format(x[i], px[i]))
        outFile.close()

    if verbose == True:
        print('Saved data to {}'.format(fname))


def plotPDF(xDate, pxDate, xAge, pxAge, refDate, ageFactor):
    '''
    Plot date and equivalent age PDFs.
    '''

    # Calculate CDF
    PxAge=cumtrapz(pxAge, xAge, initial=0)

    # Set up figure
    fig = plt.figure()

    # Plot date
    axDate = fig.add_subplot(211)
    axDate.plot(xDate, pxDate,'k')
    axDate.invert_xaxis()
    axDate.set_yticks([])
    axDate.set_ylabel('Cal date CE/BCE')

    # Plot age
    axAge = fig.add_subplot(212)
    axAge.plot(xAge, pxAge/pxAge.max(), 'k', linewidth=2, zorder=2, label='PDF')
    axAge.plot(xAge, PxAge, 'b', linewidth=2, zorder=1, label='CDF')
    axAge.set_yticks([])
    axAge.set_ylabel('Age ({:1.1e} yr)\nbefore {:.0f}'.format(ageFactor, refDate))
    axAge.legend()



### MAIN ---
if __name__ == '__main__':
    # Gather inputs
    inps = cmdParser()

    # Confirm output directory exists
    confirmOutputDir(inps.outName)

    # Load input data
    xDate, pxDate = loadDates(inps.datefile, verbose=inps.verbose)

    # Convert date to age
    xAge, pxAge = date2age(xDate, pxDate, refDate=inps.refDate, ageFactor=inps.ageFactor, verbose=inps.verbose)

    # Smooth if requested
    if inps.smoothing:
        pxAge = smoothPDF(xAge, pxAge, ktype=inps.smoothingKernel, kwidth=inps.smoothing)
        if inps.verbose == True: print('Smoothed function')

    # Save to text file
    saveOutputs(xAge, pxAge, inps.outName, verbose=inps.verbose)

    # Plot if requested
    if inps.plot:
        plotPDF(xDate, pxDate, xAge, pxAge, inps.refDate, inps.ageFactor)
        plt.show()