'''
** MCMC Incremental Slip Rate Calculator **
Objects and support functions for MCMC slip rate calculator

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from PDFanalysis import *
from array2pdf import arrayHist, arrayKDE
from PDFanalysis import IQRpdf, HPDpdf



### AGE, DISPLACEMENT, SLIP RATE CLASSES ---
## Age datum class
class ageDatum:
    '''
    PDF representing an age measurement. Units are thousands of years.
    '''
    def __init__(self, name):
        # Basic parameters
        self.name = name

    # Read in data from csv or similar text file
    def readFromFile(self, filepath):
        # Format raw data
        data = np.loadtxt(filepath)
        self.ages = data[:,0]
        self.probs = data[:,1]

    # Format for use in slip rate analysis
    def format(self, verbose=False, plot=False):
        if verbose == True:
            print('Formatting {:s} for slip rate analysis'.format(self.name))

        # Sum to CDF
        P=np.trapz(self.probs, self.ages) # area under curve
        self.probs /= P # normalize area to 1.0
        cdf_all=cumtrapz(self.probs, self.ages, initial=0)

        # Use only unique values
        uniqueVals, uniqueNdx = np.unique(cdf_all, return_index=True)
        self.ages = self.ages[uniqueNdx]  # unique ages
        self.probs = self.probs[uniqueNdx]  # unique age probs
        self.probs = self.probs/np.trapz(self.probs, self.ages)  # re-norm area
        self.cdf = cumtrapz(self.probs, self.ages, initial=0)  # re-calc CDF
        if verbose == True:
            print('\t...limiting to unique values')
            print('\t...final CDF value: {}'.format(self.cdf[-1]))

        # Inverse interpolation function
        #    use cdf as "x" value for inverse interpolation
        #    leave kind as linear to avoid values < 0 or > 1
        self.InvCDF = interp1d(self.cdf, self.ages, kind='linear')
        if verbose == True:
            print('\t...built inverse interpolation function')

        # Basic statistics
        self.lowerLimit, self.median, self.upperLimit = self.InvCDF([0.025, 0.5, 0.975])
        if verbose == True:
            print('\t95.45 % limits: {0:.5f}; {1:.5f}'.format(self.lowerLimit, self.upperLimit))

        # Plot if requested
        if plot == True:
            fig, ax = plt.subplots()
            ax.plot(self.ages, 0.8*self.probs/self.probs.max(), color=(0.3,0.3,0.6), label='PDF')
            ax.plot(self.ages, self.cdf,color='k', linewidth=2, label='CDF')
            ax.plot([self.lowerLimit, self.upperLimit], [0,0], 'rx', label='95% bounds')
            ax.set_ylim([-0.1,1.1])
            ax.set_xlabel('age'); ax.set_ylabel('prob')
            ax.set_title('INPUT: {:s}'.format(self.name))
            ax.legend()


## Displacement datum class
class dspDatum:
    '''
    PDF representing a displacement measurement. Units are meters.
    '''
    def __init__(self, name):
        # Basic parameters
        self.name = name

    # Read in data from csv or similar text file
    def readFromFile(self, filepath):
        # Format raw data
        data = np.loadtxt(filepath)
        self.dsps = data[:,0]
        self.probs = data[:,1]

    # Format for use in slip rate analysis
    def format(self, verbose=False, plot=False):
        if verbose == True:
            print('Formatting {:s} for slip rate analysis'.format(self.name))

        # Sum to CDF
        P = np.trapz(self.probs, self.dsps)  # area under curve
        self.probs /= P  # normalize area to 1.0
        cdf_all = cumtrapz(self.probs, self.dsps, initial=0)

        # Use only unique values
        uniqueVals, uniqueNdx = np.unique(cdf_all, return_index=True)
        self.dsps = self.dsps[uniqueNdx]  # unique displacements
        self.probs = self.probs[uniqueNdx]  # unique displacement probs
        self.probs = self.probs/np.trapz(self.probs, self.dsps)  # re-norm area
        self.cdf = cumtrapz(self.probs, self.dsps, initial=0)  # re-calc CDF
        if verbose == True:
            print('\t...limiting to unique values')
            print('\t...final CDF value: {}'.format(self.cdf[-1]))

        # Inverse interpolation function
        #    use cdf as "x" value for inverse interpolation
        #    leave kind as linear to avoid values < 0 or > 1
        self.InvCDF = interp1d(self.cdf, self.dsps, kind='linear')
        if verbose == True: print('\t...built inverse interpolation function')

        # Basic stats
        self.lowerLimit, self.median, self.upperLimit = self.InvCDF([0.025, 0.5, 0.975])
        if verbose == True:
            print('\t95.45 % limits: {0:.5f}; {1:.5f}'.format(self.lowerLimit, self.upperLimit))

        # Plot if requested
        if plot == True:
            fig, ax = plt.subplots()
            ax.plot(self.dsps, 0.8*self.probs/self.probs.max(), color=(0.3,0.3,0.6), label='PDF')
            ax.plot(self.dsps, self.cdf, color='k', linewidth=2, label='CDF')
            ax.plot([self.lowerLimit, self.upperLimit], [0,0], 'rx', label='95% bounds')
            ax.set_ylim([-0.1, 1.1])
            ax.set_xlabel('displacement'); ax.set_ylabel('prob')
            ax.set_title('INPUT: {:s}'.format(self.name))
            ax.legend()


## Slip rate object
class incrSlipRate:
    '''
    Incremental slip rate object, with options to convert sample picks into
     pseudo-continuous PDFs.
    '''
    def __init__(self, name):
        self.name = name

    # Convert picks to PDF
    def picks2PDF(self, RatePicks, method, stepsize, smoothingKernel=None, kernelWidth=2, verbose=False, plot=False):
        '''
        Convert slip rate picks to PDF using arrayHist or arrayKDE methods.
        '''
        # Use histogram method
        if method.lower() in ['hist', 'histogram']:
            self.rates, self.probs = arrayHist(RatePicks, stepsize,
                smoothingKernel, kernelWidth, verbose, plot)
        # Use kernel density method
        elif method.lower() in ['kde', 'kernel']:
            self.rates, self.probs = arrayKDE(RatePicks, stepsize,
                smoothingKernel, kernelWidth, verbose,plot)
        else:
            print('Choose PDF conversion method: \'histogram\'/\'kde\'')
            exit()

    # Analyze PDF
    def analyzePDF(self, method='IQR', confidence=68.27):
        '''
        Analyze the slip rate PDF to retrieve values using interquantile
         range (IQR) or highest posterior density (HPD).
        '''
        self.analysisMethod = method

        if method.upper() in ['IQR']:
            PDF = IQRpdf(self.rates, self.probs, confidence=confidence)

            # Translate results
            self.median = PDF.median
            self.xIQR = PDF.xIQR
            self.pxIQR = PDF.pxIQR
            self.lowerValue = PDF.lowerValue  # record lower value
            self.upperValue = PDF.upperValue  # record upper value
        elif method.upper() in ['HPD']:
            PDF = HPDpdf(self.rates, self.probs, confidence=confidence)

            # Translate results
            self.mode = PDF.mode  # peak value
            self.x_clusters = PDF.x_clusters  # record cluster values
            self.px_clusters = PDF.px_clusters  # record cluster probs
            self.Nclusters = len(PDF.x_clusters)  # number of clusters
            self.lowerValue = PDF.lowestValue  # record lower value
            self.upperValue = PDF.highestValue  # record upper value
        else:
            print('Choose PDF analysis method: \'IQR\'/\'HPD\'')
            exit()

    # Save to file
    def save2txt(self, outName):
        '''
        Save to text file.
        '''
        # Check outname formatting
        if outName[-4:] != '.txt': outName += '.txt'
        
        # Save to same format as input data (rates, probabilities)
        with open(outName, 'w') as outFile:
            outFile.write('# Rate\tProbability')  # header
            np.savetxt(outFile, np.column_stack([self.rates, self.probs]))  # save numpy array