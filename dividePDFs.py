#!/usr/bin/env python3
'''
    ** MCMC Incremental Slip Rate Calculator **
    This function applies a form convolution to analytically find the
     quotient two PDFs.

    Rob Zinke 2019, 2020
'''

### IMPORT MODULES ---
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


### PARSER ---
Description = '''Find the quotient of two PDFs.'''

Examples='''EXAMPLES'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    parser.add_argument(dest='numerPDFname', type=str,
        help='PDF to be used as the numerator')
    parser.add_argument(dest='denomPDFname', type=str,
        help='PDF to be used as the denominator')
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
class PDFquotient:
    '''
        Find analytically the quotient of two quantities described by two PDFs.
    '''
    def __init__(self,Xnumer,pXnumer,Xdenom,pXdenom,verbose=False):
        # Record data
        self.verbose = verbose

        # Format data
        self.Xnumer, self.pXnumer = self.__formatPDF__(Xnumer,pXnumer)
        self.Xdenom, self.pXdenom = self.__formatPDF__(Xdenom,pXdenom)

        # Compute quotient
        self.__dividePDFs__()

    def __formatPDF__(self,X,pX):
        '''
            Format an n x 2 array into separate parts and ensure unit mass.
        '''
        # Normalize area to 1.0
        pX = pX/np.trapz(pX,X)

        return X, pX

    def __dividePDFs__(self):
        '''
            Compute quotient Q, as probability function pQ.
        '''
        # Establish quotient axis, Q
        minQ = self.Xnumer.min()/self.Xdenom.max()
        maxQ = self.Xnumer.max()/self.Xdenom.min()
        lenQ = np.min([len(self.Xnumer),len(self.Xdenom)])
        self.Q = np.linspace(minQ, maxQ, lenQ)

        # Establish interpolation function for numerator
        Inumer = interp1d(self.Xnumer,self.pXnumer,
            kind='linear',bounds_error=False,fill_value=0)

        # Compute convolution\
        self.pQ = []
        for q in self.Q:
            # Equivalent numerator at each Q * denominator
            Pnumer = Inumer(q*self.Xdenom)
            pQ = np.sum(self.pXdenom*Pnumer*self.Xdenom)
            self.pQ.append(pQ)
        self.pQ = np.array(self.pQ)

        # Normalize area to unit mass
        self.pQ = self.pQ/np.trapz(self.pQ,self.Q)

    def plot(self,title=None):
        '''
            Plot raw data and difference PDF.
        '''
        # Establish figure
        self.Fig = plt.figure()

        # Plot raw data
        self.axNumer = self.Fig.add_subplot(311)
        self.axNumer.plot(self.Xnumer,self.pXnumer,
            color='b', linewidth=2, label='Numerator PDF')
        self.axNumer.legend()

        self.axDenom = self.Fig.add_subplot(312)
        self.axDenom.plot(self.Xdenom,self.pXdenom,
            color='g', linewidth=2, label='Denominator PDF')
        self.axDenom.legend()

        # Plot quotient
        self.axQuot = self.Fig.add_subplot(313)
        self.axQuot.plot(self.Q,self.pQ,
            color = 'k', linewidth=3, label = 'Quotient')
        self.axQuot.legend()

        # Format plot
        if title:
            self.Fig.suptitle(title)

    def savePDF(self,outName):
        '''
            Save to file as n x 2 array.
        '''
        outName += '.txt'
        with open(outName,'w') as outFile:
            outFile.write('# Value,\tProbability\n')
            for i in range(len(self.Q)):
                outFile.write('{0:f}\t{1:f}\n'.format(self.Q[i],self.pQ[i]))
            outFile.close()



### MAIN ---
if __name__ == '__main__':
    inps = cmdParser() # gather inputs

    # Load PDFs from file
    numerPDF = np.loadtxt(inps.numerPDFname)
    denomPDF = np.loadtxt(inps.denomPDFname)

    Xnumer = numerPDF[:,0]; pXnumer = numerPDF[:,1]
    Xdenom = denomPDF[:,0]; pXdenom = denomPDF[:,1]


    # Difference PDFs
    quot = PDFquotient(Xnumer,pXnumer,Xdenom,pXdenom,verbose=inps.verbose)

    # Save to file
    if inps.outName:
        quot.savePDF(inps.outName)

    # Plot if requested
    if inps.plot == True:
        quot.plot()

    plt.show()