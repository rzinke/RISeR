#!/usr/bin/env python3
'''
** MCMC Incremental Slip Rate Calculator **
This function applies a form convolution to analytically find the
 difference between two PDFs.

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


### PARSER ---
Description = '''Find the "delta" difference between two PDFs.'''

Examples='''EXAMPLES
# Subtract Sample1_age from Sampl2_age
differencePDFs.py Sample2_age.txt Sample1_age.txt -o Sample12_ageDiff -v -p
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

        # Create common axis and resample PDFs
        self.__establishCommonAxis__(X1, X2)
        self.__resamplePDFs__(X1, pX1, X2, pX2)

        # Compute difference
        self.__diffPDF__()

    def __establishCommonAxis__(self, X1, X2):
        '''
        Establish a common x-axis.
        '''
        # Function parameters
        self.minX = np.min([X1.min(), X2.min()])
        self.maxX = np.max([X1.max(), X2.max()])
        self.Xstep = np.min([np.diff(X1).min(), np.diff(X2).min()])

        # Establish common axis
        self.x = np.arange(self.minX, self. maxX+self.Xstep, self.Xstep)
        self.N = len(self.x)

    def __resamplePDFs__(self, X1, pX1, X2, pX2):
        '''
        Resample PDF onto specified x-axis.
        '''
        I1 = interp1d(X1, pX1, kind='linear', bounds_error=False, fill_value=0)
        I2 = interp1d(X2, pX2, kind='linear', bounds_error=False, fill_value=0)

        self.pX1 = I1(self.x)
        self.pX2 = I2(self.x)

    def __diffPDF__(self):
        '''
        Compute probability of differenes bewteen PDFs.
        '''
        # Array of difference values
        maxD = self.x.max()-self.x.min()
        self.D = np.arange(0, maxD+self.Xstep, self.Xstep)
        self.pD = np.zeros(self.D.shape)

        # Convolution
        for i in range(2, self.N):
            # Array of difference values X1 - X2
            D = self.x[i]-self.x[:i]

            # Joint probability of pX1[i] and all pX2 values less than X1
            pD = self.pX1[i]*self.pX2[:i]

            # Interpolate function at points in difference array
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
        self.axRaw.plot(self.x, self.pX1, color='b', linewidth=2, label='PDF1')
        self.axRaw.plot(self.x, self.pX2, color='g', linewidth=2, label='PDF2')
        self.axRaw.legend()

        # Plot difference
        self.axDiff = self.fig.add_subplot(212)
        self.axDiff.plot(self.D, self.pD, color = 'k', linewidth=3, label = 'difference')

        # Format plot
        self.axRaw.set_title('Difference (blue - green)')
        self.axRaw.set_ylabel('Input functions')
        self.axDiff.set_ylabel('Difference')
        if title: self.fig.suptitle(title)

    def savePDF(self, outName):
        '''
        Save to file as n x 2 array.
        '''
        outName += '.txt'
        with open(outName, 'w') as outFile:
            outFile.write('# Value,\tProbability\n')
            for i in range(self.N):
                outFile.write('{0:f}\t{1:f}\n'.format(self.D[i], self.pD[i]))
            outFile.close()


### MAIN ---
if __name__ == '__main__':
    inps = cmdParser()  # gather inputs

    # Load PDFs from file
    PDF1 = np.loadtxt(inps.PDF1name)
    PDF2 = np.loadtxt(inps.PDF2name)

    X1 = PDF1[:,0]; pX1 = PDF1[:,1]
    X2 = PDF2[:,0]; pX2 = PDF2[:,1]

    # Difference PDFs
    diff = PDFdiff(X1, pX1, X2, pX2, verbose=inps.verbose)

    # Save to file
    if inps.outName:
        diff.savePDF(inps.outName)

    # Plot if requested
    if inps.plot == True:
        diff.plot()

    plt.show()
