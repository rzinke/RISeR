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
from MCresampling import MCMCresample


### PARSER ---
Description = '''Calculate incremental slip rates for a dated fault slip history. Inputs are
are provided as probability density functions (PDFs) representing the age and displacement of each
marker. From those PDFs, random samples are drawn, and samples resulting in negative slip rates
are discarded. From the n valid random draws, incremental slip rates are calculated.

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

calcSlipRates.py DspAgeData.yaml
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
    requiredArgs.add_argument('-n','--Nsamples', dest='Nsamples', type=int, default=10000,
        help='Number of samples picked in MC run (default 1000; more is often better).')
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
    detailMCargs = parser.add_argument_group('DETAILED MC ARGUMENTS')
    detailMCargs.add_argument('--seed', dest='seed', type=float, default=0,
        help='Seed value for random number generator. Default = 0')
    detailMCargs.add_argument('--max-rate', dest='maxRate', type=float, default=1E4,
        help='Maximum rate considered in MC analysis. Units are <dispalcement units> per <age units>. Default = 10,000')

    detailAnalysisArgs = parser.add_argument_group('DETAILED SLIP RATE ANALYSIS ARGUMENTS')
    detailAnalysisArgs.add_argument('--pdf-method', dest='pdfMethod', type=str, default='kde',
        help='Method used for transforming rate picks into slip rate PDF [\'hist\'/\'kde\']. Default = kde')
    detailAnalysisArgs.add_argument('--rate-step', dest='rateStep', type=float, default=0.01,
        help='Step size for slip rate functions. Default = 0.01.')
    detailAnalysisArgs.add_argument('--smoothing-kernel', dest='smoothingKernel', type=str, default=None,
        help='Smoothing kernel type of slip rate functions. [None/mean/gauss]. Default = None')
    detailAnalysisArgs.add_argument('--kernel-width', dest='kernelWidth', type=int, default=2,
        help='Smoothing kernel width of slip rate functions.')
    detailAnalysisArgs.add_argument('--pdf-analysis', dest='pdfAnalysis', type=str, default='IQR',
        help='Method for analyzing slip rate PDFs. \'IQR\' for interquantile range; \'HPD\' for highest posterior density. Default = IQR')
    detailAnalysisArgs.add_argument('--rate-confidence', dest='rateConfidence', type=float, default=68.27,
        help='Confidence range for slip rate PDF reporting, [percent, e.g., 68.27, 95.45]. Default = 68.27')

    detailFigureArgs = parser.add_argument_group('DETAILED SLIP RATE PLOT ARGUMENTS')
    detailFigureArgs.add_argument('--max-picks', dest='maxPicks', type=int, default=500,
        help='Max number of picks to plot on MC results figure. Default = 500')
    detailFigureArgs.add_argument('--max-rate2plot', dest='maxRate2plot', type=float, default=None,
        help='Maximum spreading rate to plot (unlike -max_rate, this will not affect calculations). Default = None')
    detailFigureArgs.add_argument('--plot-inputs', dest='plotInputs', action='store_true',
        help='Plot inputs')
    return parser

def cmdParser(inps_args=None):
    parser = createParser()
    return parser.parse_args(inps_args)



### FORMATTING FUNCTIONS ---
## Output text file
def startTXTfile(outName,Nsamples):
    '''
        Ascribe the basic parameters to an output text file.
    '''
    # Construct filename
    txtName = '{}_Slip_Rate_Report.txt'.format(outName)

    # Output string
    openingReport='Incremental slip rates based on {} samples\n'.format(Nsamples)

    # Establish file
    with open(txtName,'w') as TXTout:
        TXTout.write(openingReport)

    return txtName



### STATISTICAL FUNCTIONS ---
## Percentiles from raw picks
def rawPercentiles(DspAgeData,RatePicks,txtName=None):
    '''
        Compute the percentiles of slip rate picks.
    '''
    # Parameters
    dataNames = list(DspAgeData.keys())
    m = RatePicks.shape[0] # nb incremental rates
    Nsamples = RatePicks.shape[1]
    percentiles=[50-68.27/2, 50, 50+68.27/2]

    print('Slip rate statistics from {} samples:'.format(Nsamples))
    print('Raw slip rates: (Interval: median values and 68% confidence)')
    if txtName:
        with open(txtName,'a') as TXTout:
            TXTout.write('\nIncremental slip rates based on percentiles of slip \
rate picks:\n')

    for i in range(m):
        # Formulate interval name
        intvl = '{}-{}'.format(dataNames[i],dataNames[i+1])

        # Calculate percentile
        pct = np.percentile(RatePicks[i,:],percentiles)
        median=pct[1]
        high_err=pct[2]-pct[1]
        low_err=pct[1]-pct[0]
        rawStats = '{0}: {1:.2f} +{2:.2f} -{3:.2f}'.format(intvl,
            median,high_err,low_err)
        print(rawStats)

        # Update text file if requested
        if txtName:
            with open(txtName,'a') as TXTout:
                TXTout.write(rawStats+'\n')


## Convert to PDF
def convert2PDF(DspAgeData, RatePicks, method, stepsize,
    smoothingKernel=None, kernelWidth=2,
    verbose=False):
    '''
        Convert rate picks to a probability function using the specified
         method.
    '''
    # Parameters
    dataNames = list(DspAgeData.keys())
    m = RatePicks.shape[0]

    if verbose == True:
        print('*******************************')
        print('Converting slip rate picks to PDFs using {} method'.\
            format(method))

    intervals = []
    Rates = {} # dictionary of object instances
    for i in range(m):
        # Formulate interval name
        intvl = '{}-{}'.format(dataNames[i],dataNames[i+1])
        intervals.append(intvl)

        # Create slip rate object
        Rates[intvl] = incrSlipRate(name=intvl)

        # Convert picks to PDF
        Rates[intvl].picks2PDF(RatePicks[i,:],method,stepsize,
            smoothingKernel,kernelWidth)

    return Rates


## Analyze PDFs
def analyzePDFs(Rates,method,verbose=False):
    '''
        Loop through slip rate PDFs to retrieve values using interquantile range
         (IQR) or highest posterior density (HPD).
    '''
    if verbose == True:
        print('*******************************')
        print('Extracting slip rate values from PDFs using {} method'.\
            format(method))

    for intvl in Rates.keys():
        Rates[intvl].analyzePDF(method=method)


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
    txtName = startTXTfile(inps.outName,inps.Nsamples)

    # Load data from YAML file
    DspAgeData = loadDspAgeInputs(inps.dataFile,
        verbose = inps.verbose,
        plotInputs = inps.plotInputs)

    # Plot raw data
    plotRawData(DspAgeData,
        label = inps.labelMarkers,
        outName = inps.outName)


    ## Monte Carlo resamping
    AgePicks,DspPicks,RatePicks=MCMCresample(DspAgeData,inps.Nsamples,
        condition='standard',maxRate=inps.maxRate,
        seed_value=inps.seed,
        verbose=inps.verbose,
        outName=inps.outName)

    # Plot MC results
    plotMCresults(DspAgeData, AgePicks, DspPicks, maxPicks=inps.maxPicks,
        outName=inps.outName)

    # Compute statistics based on picks and save to file
    rawPercentiles(DspAgeData, RatePicks, txtName=txtName)


    ## Convert MC results to PDFs
    Rates=convert2PDF(DspAgeData,RatePicks,method=inps.pdfMethod,
        stepsize=inps.rateStep,
        smoothingKernel=inps.smoothingKernel,kernelWidth=inps.kernelWidth,
        verbose=inps.verbose)


    ## Analyze PDFs
    # Compute statistics based on PDFs
    analyzePDFs(Rates,method=inps.pdfAnalysis,verbose=inps.verbose)

    # Plot incremental slip rates on same plot
    plotIncSlipRates(Rates,inps.pdfAnalysis,
        plotMax=inps.maxRate2plot,
        outName=inps.outName)

    # Print intervals to file
    printIncSlipRates(Rates,inps.pdfAnalysis,txtName)


    if inps.plotOutputs == True:
        plt.show()