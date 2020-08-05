#!/usr/bin/env python3
'''
    ** MCMC Incremental Slip Rate Calculator **
    This function will compute the incremental slip rates of a data
     set consisting of offsets and age markers described as
     probability density functions (PDFs).

    Rob Zinke 2019, 2020
'''

### IMPORT MODULES ---
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from dataLoading import confirmOutputDir, loadDspAgeInputs
from slipRateObjects import incrSlipRate
from plottingFunctions import plotRawData, plotMCresults, plotIncSlipRates
from differencePDFs import PDFdiff
from dividePDFs import PDFquotient


### PARSER ---
Description = '''Calculate incremental slip rates for a dated fault slip history.
Inputs are provided as probability density functions (PDFs) representing the age
and displacement of each marker. An exact computation of the slip rate PDF is
found using specialized convolution functions.

Formula for computing slip rate is eqn 7 from Bird, 2007.

REQUIRED INPUT
This routine requires a YAML file (e.g., data.yaml) which lists each dated displacement marker
in order from youngest and least-offset, to oldest and most-offset. Each entry gives the
marker name, followed by a dictionary-like entry specifying the path to the age PDF file, and
the path to the displacement PDF file.

For example: The file data.yaml might contain
# Commented description of data scheme
T2/T3 riser: {\"ageFile\": \"T2T3age.txt\", \"dspFile\": \"T2T3dsp.txt\"}
T1/T2 riser: {\"ageFile\": \"T1T2age.txt\", \"dspFile\": \"T1T2dsp.txt\"}

Where T2/T3 is younger and less-offset than T1/T2. Change the .txt filenames to the relative
or absolute paths to the PDF files, accordingly.
Note: More than one entry must be present to calculate incremental slip rates.
'''

Examples = '''EXAMPLES
From the Examples/SimpleExample folder

analyticalSlipRates.py DspAgeData.yaml
'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    # Required
    requiredArgs = parser.add_argument_group('ESSENTIAL ARGUMENTS')
    requiredArgs.add_argument(dest='dataFile', type=str,
        help='Data file in YAML format. Each entry provides a unique marker \
name, with ageFile and dspFile specified. Youngest, least offset features at \
the top; oldest, most-offset features at the bottom.')
    requiredArgs.add_argument('-n','--number-points', dest='nPts', type=int, default=1000,
        help='Number of data points in quotient function.')
    requiredArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out',
        help='Head name for outputs (no extension)')

    # Generic arguments
    generalArgs = parser.add_argument_group('GENERIC ARGUMENTS')
    generalArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', default=False,
        help='Verbose mode')
    generalArgs.add_argument('-p','--plot-outputs', dest='plotOutputs', action='store_true',
        help='Plot outputs')
    generalArgs.add_argument('-l','--label-markers', dest='labelMarkers', action='store_true',
        help='Label dated displacement markers on raw data plot')

    # Fine-tuning
    detailAnalysisArgs = parser.add_argument_group('DETAILED SLIP RATE ANALYSIS ARGUMENTS')
    detailAnalysisArgs.add_argument('--pdf-analysis', dest='pdfAnalysis', type=str, default='IQR',
        help='Method for analyzing slip rate PDFs. \'IQR\' for interquantile range; \'HPD\' for highest posterior density. Default = IQR')
    detailAnalysisArgs.add_argument('--rate-confidence', dest='rateConfidence', type=float, default=68.27,
        help='Confidence range for slip rate PDF reporting, [percent, e.g., 68.27, 95.45]. Default = 68.27')

    detailFigureArgs = parser.add_argument_group('DETAILED SLIP RATE PLOT ARGUMENTS')
    detailFigureArgs.add_argument('--max-rate2plot', dest='maxRate2plot', type=float, default=None,
        help='Maximum spreading rate to plot (unlike -max_rate, this will not affect calculations). Default = None')
    detailFigureArgs.add_argument('--plot-inputs', dest='plotInputs', action='store_true',
        help='Plot inputs')
    return parser

def cmdParser(inps_args=None):
    parser = createParser()
    return parser.parse_args(inps_args)



### ANCILLARY FUNCTIONS ---
## Output text file
def startTXTfile(outName):
    '''
        Ascribe the basic parameters to an output text file.
    '''
    # Construct filename
    txtName = '{}_Slip_Rate_Report.txt'.format(outName)

    # Output string
    openingReport='Incremental slip rates based on analytical method\n'

    # Establish file
    with open(txtName,'w') as TXTout:
        TXTout.write(openingReport)

    return txtName


## Print incremental slip rate stats to text file
def printIncSlipRates(Rates,analysisMethod,txtName):
    '''
        Append incremental slip rate statistics to the text file.
    '''
    with open(txtName,'a') as TXTout:
        intervalNames = list(Rates.keys())
        TXTout.write('\nIncremental slip rates based on PDF analysis\n')

        # Text string
        txtStr = '{0}: {1:.2f} +{2:.2f} -{3:.2f}\n'

        # Loop through intervals
        for intvl in intervalNames:
            Rate = Rates[intvl]
            # IQR method
            if analysisMethod == 'IQR':
                # Record confidence range
                TXTout.write(txtStr.format(intvl,Rate.median,
                    Rate.upperValue-Rate.median,Rate.median-Rate.lowerValue))
            # HPD method
            elif analysisMethod == 'HPD':
                # Record overall confidence range
                TXTout.write(txtStr.format(intvl,Rate.mode,
                    Rate.upperValue-Rate.mode,Rate.mode-Rate.lowerValue))
                # Record clusters
                TXTout.write('\tRanges:\n')
                rangeStr = '\t{0:.2f} to {1:.2f}\n'
                for i in range(Rate.Nclusters):
                    TXTout.write(rangeStr.format(Rate.x_clusters[i].min(),
                        Rate.x_clusters[i].max()))



### MAIN ---
if __name__ == '__main__':
    inps = cmdParser()


    ## Load files
    # Check output directory exists
    inps.outName = confirmOutputDir(inps.outName)

    # Start text file for results
    txtName = startTXTfile(inps.outName)

    # Load data from YAML file
    DspAgeData = loadDspAgeInputs(inps.dataFile,
        verbose = inps.verbose,
        plotInputs = inps.plotInputs)

    # Plot raw data
    plotRawData(DspAgeData,
        label = inps.labelMarkers,
        outName = inps.outName)


    ## Differences between data points
    # Compute age differences
    Rates = {}

    # Report if requested
    if inps.verbose == True:
        print('*******************************')
        print('Computing differences bewteen PDFs')

    # Loop through intervals
    m = len(DspAgeData.keys())
    for i in range(m-1):
        youngerName = list(DspAgeData.keys())[i]
        olderName = list(DspAgeData.keys())[i+1]
        intvl = '{}-{}'.format(youngerName,olderName)
        if inps.verbose == True: print(intvl)

        # Difference ages
        youngerPDF = DspAgeData[youngerName]
        olderPDF = DspAgeData[olderName]
        deltaAge = PDFdiff(olderPDF['Age'].ages,olderPDF['Age'].probs,
            youngerPDF['Age'].ages,youngerPDF['Age'].probs,
            verbose=inps.verbose)

        # Remove points with zero probability
        w = (deltaAge.D>0.5) # indices where probability gt zero
        deltaAge.D = deltaAge.D[w]
        deltaAge.pD = deltaAge.pD[w]

        # Plot age differences if requested
        if inps.plotInputs == True:
            ageDiffLabel = intvl+' age diff'
            deltaAge.plot(title=ageDiffLabel)


        # Difference displacements
        deltaDsp = PDFdiff(olderPDF['Dsp'].dsps,olderPDF['Dsp'].probs,
            youngerPDF['Dsp'].dsps,youngerPDF['Dsp'].probs,
            verbose=inps.verbose)

        # Plot displacement differences if requested
        if inps.plotInputs == True:
            dspDiffLabel = intvl+' dsp diff'
            deltaDsp.plot(title=dspDiffLabel)


        # Compute incremental slip rates
        slipRate = PDFquotient(deltaDsp.D,deltaDsp.pD,
            deltaAge.D,deltaAge.pD,
            n=inps.nPts)

        # Record to dictionary
        Rates[intvl] = incrSlipRate(intvl)
        Rates[intvl].rates = slipRate.Q
        Rates[intvl].probs = slipRate.pQ


        # Analyze PDF
        Rates[intvl].analyzePDF(method=inps.pdfAnalysis)


    ## Plot incremental slip rates
    plotIncSlipRates(Rates,inps.pdfAnalysis,
        plotMax=inps.maxRate2plot,
        outName=inps.outName)


    if inps.plotOutputs == True:
        plt.show()