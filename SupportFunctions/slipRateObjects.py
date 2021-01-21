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
    def __init__(self, name):
        '''
        PDF representing an age measurement.
        '''
        # Basic parameters
        self.name = name

    def readFromFile(self, filepath):
        '''
        Read data from 2-column file (txt, csv, or similar).
        First col = age; Second col = probability.
        '''
        # Format raw data
        data = np.loadtxt(filepath)
        self.ages = data[:,0]
        self.probs = data[:,1]

    # Format for use in slip rate analysis
    def format(self, verbose=False):
        '''
        Format for use in slip rate analysis:
            Sum probabilities to CDF
            Sort for unique CDF values
            Ensure unit mass
            Build inverse interpolation function
            Compute basic statistics
        '''
        if verbose == True: print('Formatted {:s} for slip rate analysis'.format(self.name))

        # Build initial CDF
        if verbose == True: print('\t...building initial CDF')
        self.__buildCDF__()

        # Limit to unique values
        if verbose == True: print('\t...limiting to unique values')
        self.__uniqueCDF__()

        # Re-normalize mass
        self.__normalizeMass__()

        # Recompute final CDF
        self.__buildCDF__()
        if verbose == True: print('\t...final CDF value: {:f}'.format(self.cdf[-1]))

        # Build inverse interpolation function
        if verbose == True: print('\t...builing inverse interpolation function')
        self.__buildPIT__()

        # Compute basic statistics
        self.__computeStats__()
        if verbose == True:
            print('\t95.45 % limits: {0:.5f}; {1:.5f}'.format(self.lowerLimit, self.upperLimit))

    def __normalizeMass__(self):
        '''
        Normalize a PDF to unit mass.
        '''
        # Compute area
        A = np.trapz(self.probs, self.ages)

        # Normalize probability
        self.probs = self.probs/A

    def __buildCDF__(self, verbose=False):
        '''
        Sum values to cumulative distribution function (CDF).
        '''
        # Sum to CDF
        self.cdf = cumtrapz(self.probs, self.ages, initial=0)

    def __uniqueCDF__(self):
        '''
        Sort for only unique values of the CDF such that it is monotonic and not flat.
        '''
        # Fine unique values
        uniqueVals, uniqueNdx = np.unique(self.cdf, return_index=True)

        # Keep only unique values
        self.ages = self.ages[uniqueNdx]  # unique ages
        self.probs = self.probs[uniqueNdx]  # unique age probs

    def __buildPIT__(self):
        '''
        Build probability inverse transform (PIT) function.
        '''
        # Inverse interpolation function
        #    use cdf as "x" value for inverse interpolation
        #    leave kind as linear to avoid values < 0 or > 1
        self.InvCDF = interp1d(self.cdf, self.ages, kind='linear')

    def __computeStats__(self):
        '''
        Compute basic statistics.
        '''
        # Basic statistics
        self.lowerLimit, self.median, self.upperLimit = self.InvCDF([0.025, 0.5, 0.975])


    def plot(self):
        '''
        Plot for visual inspection.
        '''
        # Establish axis
        fig, ax = plt.subplots()

        # Plot PDF
        ages = np.pad(self.ages, (1, 1), 'edge')
        probs = np.pad(self.probs, (1, 1), 'constant')
        ax.fill(ages, probs/probs.max(), color=(0.6, 0.6, 0.6), label='PDF')

        # Plot CDF
        ax.plot(self.ages, self.cdf, color=(0.0, 0.3, 0.7), linewidth=2, label='CDF')

        # Plot 95% limits
        ax.axvline(self.lowerLimit, color='r', linestyle='--')
        ax.axvline(self.upperLimit, color='r', linestyle='--')
        ax.plot(0, 0, color='r', linestyle='--', label='95%')

        # Format plot
        ax.set_xlim([self.ages.min(), self.ages.max()])
        ax.set_ylim([-0.1, 1.1])
        ax.set_xlabel('age')
        ax.set_ylabel('cum. prob')
        ax.set_title('INPUT: {:s}'.format(self.name))
        ax.legend()


## Displacement datum class
class dspDatum:
    def __init__(self, name):
        '''
        PDF representing an displacement measurement.
        '''
        # Basic parameters
        self.name = name

    def readFromFile(self, filepath):
        '''
        Read data from 2-column file (txt, csv, or similar).
        First col = displacement; Second col = probability.
        '''
        # Format raw data
        data = np.loadtxt(filepath)
        self.dsps = data[:,0]
        self.probs = data[:,1]

    # Format for use in slip rate analysis
    def format(self, verbose=False):
        '''
        Format for use in slip rate analysis:
            Sum probabilities to CDF
            Sort for unique CDF values
            Ensure unit mass
            Build inverse interpolation function
            Compute basic statistics
        '''
        if verbose == True: print('Formatted {:s} for slip rate analysis'.format(self.name))

        # Build initial CDF
        if verbose == True: print('\t...building initial CDF')
        self.__buildCDF__()

        # Limit to unique values
        if verbose == True: print('\t...limiting to unique values')
        self.__uniqueCDF__()

        # Re-normalize mass
        self.__normalizeMass__()

        # Recompute final CDF
        self.__buildCDF__()
        if verbose == True: print('\t...final CDF value: {:f}'.format(self.cdf[-1]))

        # Build inverse interpolation function
        if verbose == True: print('\t...builing inverse interpolation function')
        self.__buildPIT__()

        # Compute basic statistics
        self.__computeStats__()
        if verbose == True:
            print('\t95.45 % limits: {0:.5f}; {1:.5f}'.format(self.lowerLimit, self.upperLimit))

    def __normalizeMass__(self):
        '''
        Normalize a PDF to unit mass.
        '''
        # Compute area
        A = np.trapz(self.probs, self.dsps)

        # Normalize probability
        self.probs = self.probs/A

    def __buildCDF__(self, verbose=False):
        '''
        Sum values to cumulative distribution function (CDF).
        '''
        # Sum to CDF
        self.cdf = cumtrapz(self.probs, self.dsps, initial=0)

    def __uniqueCDF__(self):
        '''
        Sort for only unique values of the CDF such that it is monotonic and not flat.
        '''
        # Fine unique values
        uniqueVals, uniqueNdx = np.unique(self.cdf, return_index=True)

        # Keep only unique values
        self.dsps = self.dsps[uniqueNdx]  # unique displacements
        self.probs = self.probs[uniqueNdx]  # unique displacement probs

    def __buildPIT__(self):
        '''
        Build probability inverse transform (PIT) function.
        '''
        # Inverse interpolation function
        #    use cdf as "x" value for inverse interpolation
        #    leave kind as linear to avoid values < 0 or > 1
        self.InvCDF = interp1d(self.cdf, self.dsps, kind='linear')

    def __computeStats__(self):
        '''
        Compute basic statistics.
        '''
        # Basic statistics
        self.lowerLimit, self.median, self.upperLimit = self.InvCDF([0.025, 0.5, 0.975])


    def plot(self):
        '''
        Plot for visual inspection.
        '''
        # Establish axis
        fig, ax = plt.subplots()

        # Plot PDF
        dsps = np.pad(self.dsps, (1, 1), 'edge')
        probs = np.pad(self.probs, (1, 1), 'constant')
        ax.fill(dsps, probs/probs.max(), color=(0.6, 0.6, 0.6), label='PDF')

        # Plot CDF
        ax.plot(self.dsps, self.cdf, color=(0.0, 0.3, 0.7), linewidth=2, label='CDF')

        # Plot 95% limits
        ax.axvline(self.lowerLimit, color='r', linestyle='--')
        ax.axvline(self.upperLimit, color='r', linestyle='--')
        ax.plot(0, 0, color='r', linestyle='--', label='95%')

        # Format plot
        ax.set_xlim([self.dsps.min(), self.dsps.max()])
        ax.set_ylim([-0.1, 1.1])
        ax.set_xlabel('displacement')
        ax.set_ylabel('cum. prob')
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