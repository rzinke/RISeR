#!/usr/bin/env python3
"""
    ** MCMC Incremental Slip Rate Calculator **
    Combine probability density functions (PDFs) as union (sum)
     for intersection (multiplication).

    Rob Zinke 2019, 2020
"""

### IMPORT MODULES ---
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from PDFanalysis import HPDpdf


### PARSER ---
Description = '''Combine probability density functions (PDFs) by pointwise summing for the statistical union or mulitplication for the statistical intersection.'''

Examples = '''EXAMPLES

# Create and combine two overlapping PDFs based on their union
makePDF.py -d gauss -v 3.0 0.3 -o age1
makePDF.py -d trap -v 3.1 3.2 3.8 4.2 -o age2
combinePDFs.py age1.txt age2.txt -o Age1-2_union -m union -p

# Combine two overlapping PDFs based on their intersection
combinePDFs.py age1.txt age2.txt -o Age1-2_inters -m intersection -p
'''

def createParser():
    import argparse
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    parser.add_argument(dest='pdfList', nargs='+',
        help='List of PDFs to combine: Filenames separated by spaces.')
    parser.add_argument('-o','--outName', dest='outName', type=str, required=True,
        help='Output name.')
    parser.add_argument('-m','--method', dest='method', default='union', type=str,
        help='Method for combining (union or intersection). Default = union.')
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode')
    parser.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot resulting PDF')
    return parser

def cmdParser(inps_args=None):
    parser=createParser()
    return parser.parse_args(inps_args)



### ANCILLARY FUNCTIONS ---
def loadPDFs(pdfList, verbose = False):
    '''
        Load PDFs from text files.
    '''
    PDFs = {}

    # Load data
    for pdf in pdfList:
        if verbose == True: print('Loading: {}'.format(os.path.basename(pdf)))
        PDFs[pdf] = np.loadtxt(pdf) # load data

    return PDFs



### PDFcombo Class ---
class PDFcombo:
    '''
        Class for combining probability density functions.
    '''
    def __init__(self, PDFs, method = 'union', verbose=False):
        '''
            Combine PDFs.
        '''
        # Initialize
        self.PDFs = PDFs
        PDFs = None
        self.verbose = verbose

        # Check for minimum/maximum values
        self.__findSampleParameters__()

        # Combine PDFs
        self.__combinePDFs__(method)


    def __findSampleParameters__(self):
        '''
            Find the min/max values and sampling rate to be used for the combo.
        '''
        xmin = []; xmax = []
        dx = []

        # Find parameters for each distribution
        for pdf in self.PDFs.values():
            # Parse data
            x = pdf[:,0]
            px = pdf[:,1]

            # Store min/max values
            xmin.append(x.min())
            xmax.append(x.max())

            # Sampling rate
            dx.append(np.mean(np.diff(x)))

        # Record values
        self.xmin = min(xmin)
        self.xmax = max(xmax)
        self.dx = min(dx)

        if self.verbose == True:
            print('min: {} / max {} / sampling {:.4f}'.format(self.xmin, self.xmax, self.dx))


    def __combinePDFs__(self, method):
        '''
            Interpolate the PDFs along a common axis and combine using specified scheme.
        '''
        # Build common axis
        self.xCombo = np.arange(self.xmin, self.xmax+self.dx, self.dx)

        # Build initial values
        if method == 'union':
            self.pxCombo=np.zeros(len(self.xCombo))
        elif method == 'intersection':
            self.pxCombo=np.ones(len(self.xCombo))

        # Loop through each PDF
        for pdf in self.PDFs.values():
            # Parse data
            x = pdf[:,0]
            px = pdf[:,1]

            # Ensure area = unit mass
            P = np.trapz(px,x)
            px /= P

            # Interpolate along common axis
            I = interp1d(x, px, bounds_error = False, fill_value=0., kind='linear')

            # Combine
            if method == 'union':
                # Union sums PDFs
                self.pxCombo+=I(self.xCombo) # add to cumulative function
            elif method == 'intersection':
                # Intersection multiplies PDFs
                self.pxCombo*=I(self.xCombo)

        # Normalize final area to unit mass
        P = np.trapz(self.pxCombo,self.xCombo)
        self.pxCombo /= P

        # Report final statistics if requested
        if self.verbose == True:
            HPDpdf(self.xCombo,self.pxCombo,confidence=95.45,verbose=True)


    def saveToFile(self, outName):
        '''
            Save to two-column file.
        '''
        fname = outName+'.txt'
        with open(fname, 'w') as outFile:
            outFile.write('# Value,\tProbability\n')
            for i in range(len(self.xCombo)):
                outFile.write('{0:f}\t{1:f}\n'.format(self.xCombo[i],self.pxCombo[i]))
            outFile.close()

        # Report if requested
        if self.verbose == True:
            print('Saved data to file: {}'.format(fname))


    def plotPDF(self, outName):
        '''
            Plot the inputs and combined PDFs.
        '''

        # Plot
        Fig = plt.figure()
        ax = Fig.add_subplot(111)
        for pdfName in self.PDFs.keys():
        	data = self.PDFs[pdfName]
        	x = data[:,0]
        	px = data[:,1]
        	ax.plot(x, px, color=(0.6,0.6,0.6))
        	ax.text(x[px==px.max()][0],px[px==px.max()][0],pdfName)
        ax.plot(self.xCombo, self.pxCombo, color='k', linewidth=2)
        ax.set_xlabel('value'); ax.set_ylabel('rel prob')
        ax.set_title('Combined PDF')

        # Save to figure
        figName = outName+'.png'
        Fig.savefig(figName)
        if self.verbose == True: print('Saved figure to {}'.format(figName))

        plt.show()



### MAIN ---
if __name__ == '__main__':
    # Gather arguments
    inps=cmdParser()

    # Load PDFs
    PDFs = loadPDFs(pdfList = inps.pdfList, verbose = inps.verbose)

    # Instantiate object
    combo = PDFcombo(PDFs, method = inps.method, verbose = inps.verbose)

    # Save to file
    combo.saveToFile(outName = inps.outName)

    # Plot if requested
    if inps.plot == True:
        combo.plotPDF(outName = inps.outName)