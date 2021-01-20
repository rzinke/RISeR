#!/usr/bin/env python3
'''
** MCMC Incremental Slip Rate Calculator **
This function applies a form of convolution to analytically find the
 difference between two PDFs.

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from dataLoading import confirmOutputDir


### PARSER ---
Description = '''Find the "delta" difference between two PDFs.'''

Examples='''EXAMPLES
# Subtract Younger_age from Older_age
differencePDFs.py Older_age.txt Younger_age.txt -o Diff_older-younger -v -p
'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    parser.add_argument(dest='PDF1name', type=str,
        help='PDF from which PDF2 is subtracted')
    parser.add_argument(dest='PDF2name', type=str,
        help='PDF to subtract from PDF1')
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



### PDF DIFFERENCE CLASS ---
class PDFdiff:
    '''
    Class for differencing two PDFs.
    '''
    def __init__(self, X1, pX1, X2, pX2, verbose=False):
        # Record data
        self.verbose = verbose

        # Establish difference axis
        self.__estbDiffAxis__(X1, X2)

        # Resample PDFs at same frequency
        self.__resamplePDFs__(X1, pX1, X2, pX2)

        # Compute difference
        self.__diffPDF__()

    def __estbDiffAxis__(self, X1, X2):
        '''
        Establish a set of "difference" values on which to map the difference probabilities.
        '''
        # Establish sampling rate
        dX1 = np.abs(np.diff(X1)).min()
        dX2 = np.abs(np.diff(X2)).min()
        self.dX = np.min([dX1, dX2])

        # Determine min/max differences
        Dmin = X1.min() - X2.max()
        Dmax = X1.max() - X2.min()

        # Establish axis
        self.D = np.arange(Dmin, Dmax+self.dX, self.dX)  # difference values
        self.nD = len(self.D)  # number of values

        # Report if requested
        if self.verbose == True:
            print('Difference parameters:')
            print('\tMin difference: {:f}'.format(self.Dmin))
            print('\tMax difference: {:f}'.format(self.Dmax))
            print('\tStep: {:f}'.format(self.dX))

    def __resamplePDFs__(self, X1, pX1, X2, pX2):
        '''
        Resample PDF onto specified x-axis.
        '''
        # Resample PDFs
        self.X1, self.pX1 = self.__resamplePDF__(X1, pX1, self.dX)
        self.X2, self.pX2 = self.__resamplePDF__(X2, pX2, self.dX)

        # Report if requested
        if self.verbose == True: print('Resampled PDFs')

    def __resamplePDF__(self, X, pX, dx):
        '''
        Resample a PDF with the given resolution dx.
        '''
        # Build interpolation function
        Intp = interp1d(X, pX, kind='linear', bounds_error=False, fill_value=0)

        # Resample values
        Xintp = np.arange(X.min(), X.max()+dx, dx)

        # Resample probabilities
        pXintp = Intp(Xintp)

        return Xintp, pXintp

    def __diffPDF__(self):
        '''
        Compute probability of differenes bewteen PDFs using convolution.
        '''
        # Establish difference probabilities
        self.pD = np.zeros(self.nD)

        # Array lengths
        nX1 = len(self.X1)
        nX2 = len(self.X2)

        # Loop through X1's
        for i in range(nX1):
            # Set/Reset difference arrays
            D = []  # difference values
            pD = []  # difference probabilities

            # Loop through X2's
            for j in range(nX2):
                # Difference
                D.append(self.X1[i] - self.X2[j])  # difference value
                pD.append(self.pX1[i] * self.pX2[j])  # difference probability

            # Add to PDF
            interpDiff = interp1d(D, pD, kind='linear', bounds_error=False, fill_value=0)
            self.pD += interpDiff(self.D)

        # Normalize area to 1.0
        self.pD = self.pD/np.trapz(self.pD, self.D)

    def plot(self, title=None):
        '''
        Plot raw data and difference PDF.
        '''
        # Establish figure
        self.fig = plt.figure()

        # Plot raw data
        self.axRaw = self.fig.add_subplot(211)
        self.axRaw.plot(self.X1, self.pX1, color='b', linewidth=2, label='PDF1')
        self.axRaw.plot(self.X2, self.pX2, color='g', linewidth=2, label='PDF2')
        self.axRaw.legend()

        # Plot difference
        self.axDiff = self.fig.add_subplot(212)
        self.axDiff.plot(self.D, self.pD, color = 'k', linewidth=3, label = 'difference')

        # Format plot
        self.axRaw.set_title('Difference (PDF1 - PDF2)')
        self.axRaw.set_yticks([])
        self.axRaw.set_ylabel('Input functions')

        self.axDiff.set_yticks([])
        self.axDiff.set_ylabel('Difference\n(PDF1 - PDF2)')

        if title: self.fig.suptitle(title)
        self.fig.tight_layout()

    def savePDF(self, outName):
        '''
        Save to file as n x 2 array.
        '''
        outName += '.txt'
        with open(outName, 'w') as outFile:
            outFile.write('# Value,\tProbability\n')
            for i in range(self.nD):
                outFile.write('{0:f}\t{1:f}\n'.format(self.D[i], self.pD[i]))
            outFile.close()


### MAIN ---
if __name__ == '__main__':
    # Gather inputs
    inps = cmdParser()

    # Confirm output directory exists
    confirmOutputDir(inps.outName)

    # Load PDFs from file
    PDF1 = np.loadtxt(inps.PDF1name)
    PDF2 = np.loadtxt(inps.PDF2name)

    X1 = PDF1[:,0]; pX1 = PDF1[:,1]
    X2 = PDF2[:,0]; pX2 = PDF2[:,1]

    # Difference PDFs
    diff = PDFdiff(X1, pX1, X2, pX2, verbose=inps.verbose)

    # Save to file
    if inps.outName: diff.savePDF(inps.outName)

    # Plot if requested
    if inps.plot == True: diff.plot()

    plt.show()
