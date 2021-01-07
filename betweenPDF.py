#!/usr/bin/env python3
'''
** MCMC Incremental Slip Rate Calculator **
Use this to create a probability density function (PDF) representing the
 likeliness of a value between two PDFs.

Rob Zinke 2019, 2020
'''

### IMPORT MODULES ---
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as intrp


### PARSER ---
Description = '''Calculate a PDF representing the values between two PDFs.
The first PDF should be smaller than the second PDF.
'''

Examples='''EXAMPLES

# Find the "between" PDF
betweenPDF.py YoungerAge.txt OlderAge.txt -o BetweenAge.txt
'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    parser.add_argument(dest='pdf1', type=str, help='Smaller PDF')
    parser.add_argument(dest='pdf2', type=str, help='Larger PDF')
    parser.add_argument('-o', '--output', dest='outName', type=str, required=True,
        help='Output file path/name')
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Print outputs to command line')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true',
        help='Show plot of output')
    return parser

def cmdParser(inpt_args=None):
    parser = createParser()
    return parser.parse_args(args=inpt_args)



### ANCILLARY FUNCTIONS ---
def loadPDF(pdfName, verbose=False):
    '''
    Load PDF from file. File should have two columns:
     value, probability
    '''
    if verbose == True: print('Loading: {:s}'.format(pdfName))

    data = np.loadtxt(pdfName)
    x = data[:,0]
    px = data[:,1]

    return x, px


def resampleCommon(x1, px1, x2, px2, verbose=False):
    '''
    Resample two PDFs onto a common axis.
    INPUTS
        x1, px1 are the values and probabilities of the smaller PDF
        x2, px2 are the values and probabilities of the larger PDF
    OUTPUTS
        x, px1, px2 are the PDFs resampled onto the same axis
    '''
    # Find resampling values
    xmin = x1.min()
    xmax = x2.max()

    dx1 = np.diff(x1).mean()
    dx2 = np.diff(x2).mean()
    dx = np.min([dx1, dx2])

    # Formulate common axis
    x = np.arange(xmin, xmax+dx, dx)

    # Resample onto common axis
    I1 = intrp.interp1d(x1, px1, bounds_error=False, fill_value=0)
    px1 = I1(x)

    I2 = intrp.interp1d(x2, px2, bounds_error=False, fill_value=0)
    px2 = I2(x)

    # Report if requested
    if verbose == True:
        print('Resampling to common axis')
        print('\tmin: {:f}; max: {:f}; dx: {:f}'.format(xmin, xmax, dx))
        print('\tlen: {:d}'.format(len(x)))

    return x, px1, px2


def plotBetween(x1, px1, x2, px2, x, pxB):
    '''
    Plot inputs and between PDF.
    '''
    # Spawn figure
    Fig, [axInputs, axBetween] = plt.subplots(nrows=2)

    # Plot inputs
    axInputs.plot(x1, px1, color='b', label='PDF 1')
    axInputs.plot(x2, px2, color='g', label='PDF 2')

    # Plot between PDF
    axBetween.plot(x, pxB, color='k', label='Between PDF')

    # Format plot
    axInputs.legend()
    axInputs.set_ylabel('Inputs')

    axBetween.legend()
    axBetween.set_ylabel('Between')


def saveOutputs(x, px, outName, verbose=False):
    '''
    Save between PDF to file, with first column representing values, 
     second column representing probabilities.
    '''
    fname = outName+'.txt'
    with open(fname,'w') as outFile:
        outFile.write('# Value,\tProbability\n')
        for i in range(len(x)):
            outFile.write('{0:f}\t{1:f}\n'.format(x[i], px[i]))
        outFile.close()

    if verbose == True:
        print('Saved data to {:s}'.format(fname))



### BETWEEN PDF ---
def betweenPDF(x, px1, px2, verbose=False):
    '''
    Given two pdfs sampled on the same axis, compute the "between" probability.
    '''
    # Parameters
    nx = len(x)
    dx = np.diff(x).mean()

    # Sum to CDFs
    Px1 = np.cumsum(px1)*dx
    Px1 /= Px1.max()

    Px2 = np.cumsum(px2)*dx
    Px2 /= Px2.max()

    # Compute between PDF
    pxB = Px1*(1-Px2)

    # Normalize to unit mass
    pxB = pxB/np.trapz(pxB,x)

    return x, pxB



### MAIN ---
if __name__ == '__main__':
    inps = cmdParser()

    # Load input data
    x1, px1 = loadPDF(inps.pdf1, verbose=inps.verbose)
    x2, px2 = loadPDF(inps.pdf2, verbose=inps.verbose)

    # Resample onto common axis
    x, px1r, px2r = resampleCommon(x1, px1, x2, px2, verbose=inps.verbose)

    # Compute "between" PDF
    x, pxB = betweenPDF(x, px1r, px2r, verbose=inps.verbose)

    # Save outputs
    saveOutputs(x, pxB, inps.outName, verbose=inps.verbose)

    # Plot if requested
    if inps.plot == True:
        plotBetween(x1, px1, x2, px2, x, pxB)
        plt.show()