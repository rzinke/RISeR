#!/usr/bin/env python3
'''
** RISeR Incremental Slip Rate Calculator **
These functions assist with saving data.

Rob Zinke 2019-2021
'''

### LOAD MODULES ---
import os
from datetime import datetime


### DIRECTORY FUNCTIONS ---
## Confirm output directory
def confirmOutputDir(outName, verbose=False):
    '''
    Confirm the existence of the output directory. If it does not exist, create it.
    '''
    # Convert outName to aboslute path
    outName = os.path.abspath(outName)

    # Get directory name
    dirName = os.path.dirname(outName)

    # Create directory if it does not exist
    if not os.path.exists(dirName):
        os.mkdir(dirName)

        # Report if requested
        if verbose == True: print('Created new directory: {:s}'.format(dirName))

    # Return outName as absolute path
    return outName



### TEXTFILE FUNCTIONS ---
class slipRateTxtFile:
    def __init__(self, outName, scheme):
        '''
        Create a text file describing the results of incremental slip rate
         computations.
        '''
        self.outName = outName
        self.scheme = scheme

        self.__startFile__()

    def __startFile__(self):
        '''
        Ascribe the basic parameters to the text file.
        '''
        # Construct filename
        self.txtName = '{:s}_Slip_Rate_Report.txt'.format(self.outName)

        # Output string
        openingReport = 'Incremental slip rate based on {:s} ({:s})\n\n'

        if self.scheme in ['analytical']:
            schemeReport = 'analytical formulation'
        elif self.scheme in ['mcmc']:
            schemeReport = 'MCMC sampling'

        # Timestamp with today's date
        now = datetime.now().strftime('%Y %m %d, %H:%M:%S')

        # Establish file
        with open(self.txtName, 'w') as TXTout:
            TXTout.write(openingReport.format(schemeReport, now))

    def append(self, txtstr):
        '''
        Append the text string txtstr to the file.
        '''
        with open(self.txtName, 'a+') as TXTout:
            TXTout.write(txtstr)


def printIncSlipRates(txtFile, Rates, analysisMethod, confidence=68.27):
    '''
    Append incremental slip rate statistics to the text file.
    '''
    intervalNames = list(Rates.keys())
    txtFile.append('\nIncremental slip rates based on PDF analysis ')
    txtFile.append('({:s}, {:.2f}% confidence)\n'.format(analysisMethod, confidence))

    # Text string
    txtStr = '{0:s}: {1:.2f} +{2:.2f} -{3:.2f}\n'

    # Loop through intervals
    for intvl in intervalNames:
        Rate = Rates[intvl]

        # IQR method
        if analysisMethod == 'IQR':
            # Record confidence range
            txtFile.append(txtStr.format(intvl, Rate.median,
                Rate.upperValue-Rate.median, Rate.median-Rate.lowerValue))

        # HPD method
        elif analysisMethod == 'HPD':
            # Record overall confidence range
            txtFile.append(txtStr.format(intvl, Rate.mode, Rate.upperValue-Rate.mode, Rate.mode-Rate.lowerValue))

            # Record clusters
            txtFile.append('\tRanges:\n')
            rangeStr = '\t{0:.2f} to {1:.2f}\n'
            for i in range(Rate.Nclusters):
                txtFile.append(rangeStr.format(Rate.x_clusters[i].min(), Rate.x_clusters[i].max()))