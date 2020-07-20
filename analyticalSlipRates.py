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
from matplotlib.patches import Rectangle
from slipRateObjects import *
from PDFanalysis import *
from differencePDFs import PDFdiff
from dividePDFs import PDFquotient


### PARSER ---
Description = '''Calculate incremental slip rates for a dated fault slip history.
Inputs are provided as probability density functions (PDFs) representing the age
and displacement of each marker. An exact computation of the slip rate PDF is
found using specialized convolution functions.

Formula for computing slip rate is eqn 7 from Bird, 2007.

REQUIRED INPUTS
* age-list is a list of file paths to the PDFs representing marker ages--youngest markers at the
    top of the list
* dsp-list is a list of file paths to the PDFs representing marker displacement measurements--youngest
    markers at the top of the list
'''

Examples = '''EXAMPLES
From the Examples/SimpleExample folder

analyticalSlipRates.py -a AgeList.txt -d DspList.txt
'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    # Required
    requiredArgs = parser.add_argument_group('ESSENTIAL ARGUMENTS')
    requiredArgs.add_argument('-a','--age-list', dest='ageListFile', type=str, required=True,
        help='Text file with one age file per line, list in order from youngest (top line) to oldest (bottom).')
    requiredArgs.add_argument('-d','--dsp-list', dest='dspListFile', type=str, required=True,
        help='Text file with one displacement file per line, list in order from smallest (top line) to largest (bottom).')
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



### DATA LOADING FUNCTIONS ---
## Confirm output directory
def confirmOutputDir(outName,verbose=False):
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

    # Return outName as absolute path
    return outName


## Load age PDFs
def loadAges(ageListFile,verbose=False,plot=False):
    # Load ages into dictionary
    with open(ageListFile,'r') as agesIn:
        ageLines=agesIn.readlines() # read lines and save to variable
        agesIn.close()

    Ages={}; ageList=[]; maxAge=0.
    for ageLine in ageLines:
        '''
            Record one dictionary entry per age file. Each entry in an
             instance of a "ageDatum" object. Append age name (path and
             suffix removed) to list for later.
        '''

        # Isolate basename
        ageLine = ageLine.strip('\n') # remove extraneous newline
        ageName = os.path.basename(ageLine)
        ageName = ageName.split('.')[0] # remove file extension
        ageList.append(ageName) # append to list

        # Spawn age object and record to dictionary
        Ages[ageName] = ageDatum(ageName)
        Ages[ageName].readFromFile(ageLine) # load data to object

        # Format for slip rate analysis
        #    builds inverse interpolation function
        Ages[ageName].format(verbose=verbose,plot=plot)

        # Set global limit for plotting
        if Ages[ageName].upperLimit>maxAge:
            maxAge = Ages[ageName].upperLimit

    return Ages, ageList, maxAge


## Load dsp PDFs
def loadDsps(dspListFile,verbose=False,plot=False):
    # Load displacements into dictionary
    with open(dspListFile,'r') as dspsIn:
        dspLines=dspsIn.readlines() # read lines and save to variable
        dspsIn.close()

    Dsps={}; dspList=[]; maxDsp=0.
    for dspLine in dspLines:
        '''
            Record one dictionary entry per displacement file. Each entry is an
             instance of a "dspDatum" object. Append displacement name (path and
             suffix removed) to list for later.
        '''

        # Isolate basename
        dspLine = dspLine.strip('\n') # remove extraneous newline
        dspName = os.path.basename(dspLine)
        dspName = dspName.split('.')[0] # remove file extension
        dspList.append(dspName) # append to list

        # Spawn age object and record to dictionary
        Dsps[dspName]=  dspDatum(dspName)
        Dsps[dspName].readFromFile(dspLine) # load data to object

        # Format for slip rate analysis
        #    builds inverse interpolation function
        Dsps[dspName].format(verbose=verbose,plot=plot)

        # Set global limit for plotting
        if Dsps[dspName].upperLimit>maxDsp:
            maxDsp = Dsps[dspName].upperLimit

    return Dsps, dspList, maxDsp


## Check inputs
def checkInputs(ageList, dspList, verbose = False):
    '''
        Check that the same number of inputs are provided for age and displacement.
    '''
    # Check
    nAges = len(ageList); nDsps=len(dspList)
    assert nAges == nDsps, 'MUST HAVE SAME NUMBER OF AGE AND DISPLACEMENT MEASUREMENTS'
    m = nAges # assign number of measurements

    # Report if requested
    if verbose == True:
        print('Detected m = {} age and displacement measurements'.format(nAges))
        # Confirm pairing
        print('Pairing (youngest at top):')
        for i in range(nAges):
            print('\t{} - {}'.format(ageList[i],dspList[i]))


## Formulate interval names
def intervalsFromDspList(dspList):
    '''
        List of intervals between one offset marker and another.
    '''
    intervalList=[]
    m = len(dspList)
    for i in range(m-1):
        interval_name='{}-{}'.format(dspList[i],dspList[i+1])
        intervalList.append(interval_name)

    return intervalList



### PLOTTING FUNCTIONS ---
## Plot raw data (whisker plot)
def plotRawData(Ages,ageList,Dsps,dspList,maxAge,maxDsp,label,outName):
    '''
        Plot raw data as whisker plot. Whiskers represent 95 % intervals.
        No MC analyses are done at this point.
    '''
    # Parameters
    m = len(ageList)

    # Establish figure
    Fig = plt.figure()
    ax = Fig.add_subplot(111)

    # Plot data
    for i in range(m):
        x_mid=Ages[ageList[i]].median
        x_err=np.array([[x_mid-Ages[ageList[i]].lowerLimit],
            [Ages[ageList[i]].upperLimit-x_mid]])
        y_mid=Dsps[dspList[i]].median
        y_err=np.array([[y_mid-Dsps[dspList[i]].lowerLimit],
            [Dsps[dspList[i]].upperLimit-y_mid]])
        ax.errorbar(x_mid,y_mid,xerr=x_err,yerr=y_err,
            color=(0.3,0.3,0.6),marker='o')

        # Label if requested
        if label == True:
            ax.text(x_mid+0.01*maxAge, y_mid+0.01*maxDsp, dspList[i])

    # Format figure
    ax.set_xlim([0,1.1*maxAge]) # x-limits
    ax.set_ylim([0,1.1*maxDsp]) # y-limits
    ax.set_xlabel('age'); ax.set_ylabel('displacement')
    ax.set_title('Raw data (95 %% limits)')
    Fig.tight_layout()

    # Save figure
    Fig.savefig('{}_Fig1_RawData.pdf'.format(outName),type='pdf')

    # Return values
    return Fig, ax


## Plot incremental slip rate results
def plotIncSlipRates(Rates,intervalList,analysis_method,outName,plot_max=None):
    # Analysis method is IQR or HPD

    # Setup figure
    Fig=plt.figure('IncrementalRates')
    ax=Fig.add_subplot(111)
    intvl_labels=[]; k=0
    # Loop through each rate
    for i in intervalList:
        # Plot full PDF
        scale_value=Rates[i].probs.max()
        ax.fill(Rates[i].rates,
            0.9*Rates[i].probs/scale_value+k,
            color=(0.4,0.4,0.4),zorder=2)
        # Plot confidence bounds
        if analysis_method=='IQR':
            ax.fill(Rates[i].xIQR,
                0.9*Rates[i].pxIQR/scale_value+k,
                color=(0.3,0.3,0.6),zorder=3)
        elif inps.pdfAnalysis=='HPD':
            for j in range(Rates[i].Nclusters):
                xHPD=Rates[i].x_clusters[j]
                xHPD=np.pad(xHPD,(1,1),'edge')
                pxHPD=Rates[i].px_clusters[j]
                pxHPD=np.pad(pxHPD,(1,1),'constant')
                ax.fill(xHPD,0.9*pxHPD/scale_value+k,
                    color=(0.3,0.3,0.6),zorder=3)
        # Format axis labels
        intvl_name='{}-\n{}'.format(i.split('-')[0],i.split('-')[1])
        intvl_labels.append(intvl_name)
        k+=1

    # Format plot
    ax.set_yticks(np.arange(0.5,m-1))
    ax.set_yticklabels(intvl_labels,rotation='vertical')
    ax.set_ylim([0,m-1])
    if plot_max:
        ax.set_xlim([0,plot_max])
    ax.set_xlabel('slip rate')
    ax.set_title('Incremental slip rates')
    Fig.tight_layout()

    # Save figure
    Fig.savefig('{}_Fig3_Incremental_slip_rates.pdf'.format(outName),type='pdf')



### MAIN ---
if __name__ == '__main__':
    inps = cmdParser()


    ## Load files
    # Check output directory exists
    inps.outName = confirmOutputDir(inps.outName)

    # Read age file
    Ages,ageList,maxAge=loadAges(inps.ageListFile,inps.verbose,inps.plotInputs)

    # Read disp file
    Dsps,dspList,maxDsp=loadDsps(inps.dspListFile,inps.verbose,inps.plotInputs)

    # Check input files have same number of lines
    checkInputs(ageList,dspList,verbose = inps.verbose)
    m = len(ageList) # record number of dated displacement measurements

    # Formulate interval names
    intervalList = intervalsFromDspList(dspList)


    ## Differences between data points
    # Compute age differences
    ageDiffList = []
    deltaAges = {}
    for i in range(m-1):
        youngerAge = ageList[i]
        olderAge = ageList[i+1]
        ageDiffLabel = '{}-{}'.format(olderAge,youngerAge)
        ageDiffList.append(ageDiffLabel)

        # Older minus younger age
        deltaAge=PDFdiff(Ages[olderAge].ages,Ages[olderAge].probs,
            Ages[youngerAge].ages,Ages[youngerAge].probs)

        # Remove points with zero probability
        w = (deltaAge.pD>0) # indices where probability gt zero
        deltaAge.D = deltaAge.D[w]
        deltaAge.pD = deltaAge.pD[w]

        # Plot if requested
        if inps.plotInputs == True:
            deltaAge.plot(title=ageDiffLabel)

        # Record to dictionary
        deltaAges[ageDiffLabel] = deltaAge

    # Compute displacement differences
    dspDiffList = []
    deltaDsps = {}
    for i in range(m-1):
        smallerDsp = dspList[i]
        largerDsp = dspList[i+1]
        dspDiffLabel = '{}-{}'.format(largerDsp,smallerDsp)
        dspDiffList.append(dspDiffLabel)

        # Larger minus smaller displacement
        deltaDsp=PDFdiff(Dsps[largerDsp].dsps,Dsps[largerDsp].probs,
            Dsps[smallerDsp].dsps,Dsps[smallerDsp].probs)

        # Remove points with zero probability
        w = (deltaDsp.pD>0)
        deltaDsp.D = deltaDsp.D[w]
        deltaDsp.pD = deltaDsp.pD[w]

        # Plot if requested
        if inps.plotInputs == True:
            deltaDsps[dspDiffLabel].plot(title=dspDiffLabel)

        # Record to dictionary
        deltaDsps[dspDiffLabel] = deltaDsp


    ## Plot raw data
    plotRawData(Ages,ageList,Dsps,dspList,
        maxAge,maxDsp,
        inps.labelMarkers,
        inps.outName)


    ## Quotients of displacements and ages for slip rate calculation
    Rates = {}
    for i in range(m-1):
        intervalName = intervalList[i]
        if inps.verbose == True:
            print(intervalName)

        # Gather age differences
        t = deltaAges[ageDiffList[i]].D
        pt = deltaAges[ageDiffList[i]].pD

        # Gather displacement differences
        u = deltaDsps[dspDiffList[i]].D
        pu = deltaDsps[dspDiffList[i]].pD

        # Compute incremental slip rates
        slipRate = PDFquotient(u,pu,t,pt)

        # Record to dictionary
        Rates[intervalName] = rateObj(intervalName)
        Rates[intervalName].rates = slipRate.Q
        Rates[intervalName].probs = slipRate.pQ


    ## Analyze PDFs
    TXTout = open('{}_slip_rate_report.txt'.format(inps.outName),'w')
    pdf_report='PDF-based statistics ({0}% confidence) calculated using {1}.'.format(inps.rateConfidence,inps.pdfAnalysis)
    print(pdf_report); TXTout.write('{}'.format(pdf_report))
    for i in range(m-1):

        # Interquantile range
        if inps.pdfAnalysis in ['IQR','quantiles','quantile']:
            '''
                More stable option, but slightly skewed toward larger values
            '''
            inps.pdfAnalysis='IQR' # confirm code for later
            PDF=IQRpdf(Rates[intervalList[i]].rates,
                Rates[intervalList[i]].probs,
                inps.rateConfidence,verbose=False)
            Rates[intervalList[i]].median=PDF.median
            Rates[intervalList[i]].xIQR=PDF.xIQR
            Rates[intervalList[i]].pxIQR=PDF.pxIQR
            Rates[intervalList[i]].lowerValue=PDF.lowerValue # record lower value
            Rates[intervalList[i]].upperValue=PDF.upperValue # record upper value

            # Report relevant stats
            median=PDF.median
            low_err=median-PDF.lowerValue
            high_err=PDF.upperValue-median

            incr_slip_rate_report='{0}: {1:.2f} +{2:.2f} -{3:.2f}'.format(intervalList[i],median,high_err,low_err)
            print(incr_slip_rate_report); TXTout.write('\n{}'.format(incr_slip_rate_report))

        # Highest posterior density
        elif inps.pdfAnalysis in ['HPD','density']:
            '''
                More representative of probable values
            '''
            inps.pdfAnalysis='HPD' # confirm code for later
            PDF=HPDpdf(Rates[intervalList[i]].rates,
                Rates[intervalList[i]].probs,
                inps.rateConfidence,verbose=False)
            Rates[intervalList[i]].mode=PDF.mode # peak value
            Rates[intervalList[i]].x_clusters=PDF.x_clusters # record cluster values
            Rates[intervalList[i]].px_clusters=PDF.px_clusters # record cluster probs
            Rates[intervalList[i]].Nclusters=len(PDF.x_clusters) # number of clusters
            Rates[intervalList[i]].lowerValue=PDF.lowestValue # record lower value
            Rates[intervalList[i]].upperValue=PDF.highestValue # record upper value

            # Report relevant stats
            peak=PDF.mode
            high_err=PDF.highestValue-peak
            low_err=peak-PDF.lowestValue

            incr_slip_rate_report='{0}: {1:.2f} +{2:.2f} -{3:.2f}'.format(intervalList[i],peak,high_err,low_err)
            print(incr_slip_rate_report); TXTout.write('\n{}'.format(incr_slip_rate_report))
            TXTout.write('\n\tRanges:')
            for i in range(PDF.nClusters):
                cluster_low_value=PDF.x_clusters[i].min()
                cluster_hi_value=PDF.x_clusters[i].max()
                cluster_range_report='\t\t({0:.2f} to {1:.2f})'.format(cluster_low_value,cluster_hi_value)
                print(cluster_range_report); TXTout.write('\n{}'.format(cluster_range_report))

    # Plot incremental slip rates on same plot
    plotIncSlipRates(Rates,intervalList,inps.pdfAnalysis,inps.outName,plot_max=inps.maxRate2plot)


    ## Close saved file
    TXTout.close()


    if inps.plotOutputs == True:
        plt.show()