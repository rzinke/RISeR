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
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from slipRateObjects import *
from MCresampling import *
from array2pdf import *
from PDFanalysis import *


### PARSER ---
def createParser():
    import argparse
    parser = argparse.ArgumentParser(description='Main function for calculating incremental slip rates.')
    # Required
    parser.add_argument('-a','--age_list', dest='ageListFile', type=str, required=True,
    	help='Text file with one age file per line, list in order from youngest (top line) to oldest (bottom).')
    parser.add_argument('-d','--dsp_list', dest='dspListFile', type=str, required=True,
    	help='Text file with one displacement file per line, list in order from smallest (top line) to largest (bottom).')
    parser.add_argument('-o','--out-name', dest='outName', type=str, default='Out',
    	help='Head name for outputs (no extension)')
    # Recommended
    parser.add_argument('-n','--Nsamples', dest='Nsamples', type=int, default=10000,
    	help='Number of samples picked in MC run (default 1000; more is often better).')
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true', default=False,
    	help='Verbose?')
    parser.add_argument('--plot-inputs', dest='plotInputs', action='store_true',
    	help='Plot inputs')
    parser.add_argument('--plot-outputs', dest='plotOutputs', action='store_true',
    	help='Plot outputs')
    # Highly optional
    parser.add_argument('--max-rate', dest='maxRate', type=float, default=1E4,
    	help='Maximum rate considered in MC analysis. Units are <dispalcement units> per <age units>. Default = 10,000')
    parser.add_argument('--seed', dest='seed', type=float, default=0, help='Seed value for random number generator. Default = 0')
    parser.add_argument('--max-picks', dest='maxPicks', type=int, default=500,
    	help='Max number of picks to plot on MC results figure. Default = 500')
    parser.add_argument('--pdf-method', dest='pdfMethod', type=str, default='kde',
    	help='Method used for transforming rate picks into slip rate PDF [\'hist\'/\'kde\']. Default = kde')
    parser.add_argument('--rate-step', dest='rateStep', type=float, default=0.01,
    	help='Step size for slip rate functions. Default = 0.01.')
    parser.add_argument('--smoothing-kernel', dest='smoothingKernel', type=str, default=None,
    	help='Smoothing kernel type of slip rate functions. [None/mean/gauss]. Default = None')
    parser.add_argument('--kernel-width', dest='kernelWidth', type=int, default=2,
    	help='Smoothing kernel width of slip rate functions.')
    parser.add_argument('--pdf-analysis', dest='pdfAnalysis', type=str, default='IQR',
    	help='Method for analyzing slip rate PDFs. \'IQR\' for interquantile range; \'HPD\' for highest posterior density. Default = IQR')
    parser.add_argument('--rate-confidence', dest='rateConfidence', type=float, default=68.27,
    	help='Confidence range for slip rate PDF reporting, [percent, e.g., 68.27, 95.45]. Default = 68.27')
    parser.add_argument('--max-rate2plot', dest='maxRate2plot', type=float, default=None,
    	help='Maximum spreading rate to plot (unlike -max_rate, this will not affect calculations). Default = None')
    return parser

def cmdParser(inpt_args=None):
    parser = createParser()
    return parser.parse_args(inpt_args)



### DATA LOADING FUNCTIONS ---
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
        ageLine=ageLine.strip('\n') # remove extraneous newline
        ageName=os.path.basename(ageLine)
        ageName=ageName.split('.')[0] # remove file extension
        ageList.append(ageName) # append to list

        # Spawn age object and record to dictionary
        Ages[ageName]=ageDatum(ageName)
        Ages[ageName].readFromFile(ageLine) # load data to object

        # Format for slip rate analysis
        #    builds inverse interpolation function
        Ages[ageName].format(verbose=verbose,plot=plot)

        # Set global limit for plotting
        if Ages[ageName].upperLimit>maxAge:
            maxAge=Ages[ageName].upperLimit

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
        dspLine=dspLine.strip('\n') # remove extraneous newline
        dspName=os.path.basename(dspLine)
        dspName=dspName.split('.')[0] # remove file extension
        dspList.append(dspName) # append to list

        # Spawn age object and record to dictionary
        Dsps[dspName]=dspDatum(dspName)
        Dsps[dspName].readFromFile(dspLine) # load data to object

        # Format for slip rate analysis
        #    builds inverse interpolation function
        Dsps[dspName].format(verbose=verbose,plot=plot)

        # Set global limit for plotting
        if Dsps[dspName].upperLimit>maxDsp:
            maxDsp=Dsps[dspName].upperLimit

    return Dsps, dspList, maxDsp



### PLOTTING FUNCTIONS ---
## Plot raw data (whisker plot)
def plotRawData(Ages,ageList,Dsps,dspList,maxAge,maxDsp,outName=None):
    # Parameters
    m=len(ageList)

    # Plot
    Fig=plt.figure()
    ax=Fig.add_subplot(111)
    for i in range(m):
        x_mid=Ages[ageList[i]].median
        x_err=np.array([[x_mid-Ages[ageList[i]].lowerLimit],
            [Ages[ageList[i]].upperLimit-x_mid]])
        y_mid=Dsps[dspList[i]].median
        y_err=np.array([[y_mid-Dsps[dspList[i]].lowerLimit],
            [Dsps[dspList[i]].upperLimit-y_mid]])
        ax.errorbar(x_mid,y_mid,xerr=x_err,yerr=y_err,
            color=(0.3,0.3,0.6),marker='o')
    ax.set_xlim([0,1.1*maxAge]) # x-limits
    ax.set_ylim([0,1.1*maxDsp]) # y-limits
    ax.set_xlabel('age'); ax.set_ylabel('displacement')
    ax.set_title('Raw data (95 % limits)')
    if outName:
        Fig.savefig('{}_Fig1_RawData.png'.format(outName),dpi=600)
        Fig.savefig('{}_Fig1_RawData.pdf'.format(outName),type='pdf')
    return Fig, ax


## Plot rectangles
def plotRectangles(Ages,ageList,Dsps,dspList):
    # Parameters
    m=len(ageList) # number of measurements

    # Establish plot
    Fig=plt.figure()
    ax=Fig.add_subplot(111)

    # Plot rectangles
    for i in range(m):
        ageLower=Ages[ageList[i]].lowerLimit # load bottom
        dspLower=Dsps[dspList[i]].lowerLimit # load left
        boxWidth=Ages[ageList[i]].upperLimit-ageLower
        boxHeight=Dsps[dspList[i]].upperLimit-dspLower
        ax.add_patch(Rectangle((ageLower,dspLower), # LL corner
            boxWidth,boxHeight, # dimensions
            edgecolor=(0.3,0.3,0.6),fill=False,zorder=3))
    return Fig, ax


## Plot MC results
def plotMCresults(Ages,ageList,Dsps,dspList,AgePicks,DspPicks, \
    xMax,yMax,maxPicks,outName=None):
    # Parameters
    m=len(ageList) # number of measurements
    n=AgePicks.shape[1] # number of picks

    # Plot rectangles
    Fig,ax=plotRectangles(Ages,ageList,Dsps,dspList)

    # Plot picks
    if n<=maxPicks:
        ax.plot(AgePicks,DspPicks,color=(0,0,0),alpha=0.1,zorder=1)
        ax.plot(AgePicks,DspPicks,color=(0,0,1),marker='.',linewidth=0,alpha=0.5,zorder=2)
    else:
        ax.plot(AgePicks[:,:maxPicks],DspPicks[:,:maxPicks],
            color=(0,0,0),alpha=0.1,zorder=1)
        ax.plot(AgePicks[:,:maxPicks],DspPicks[:,:maxPicks],
            color=(0,0,1),marker='.',linewidth=0,alpha=0.4,zorder=2)
    ax.set_xlim([0,1.1*xMax]); ax.set_ylim([0,1.1*yMax])
    ax.set_xlabel('age'); ax.set_ylabel('displacement')
    ax.set_title('MC Picks (N = {})'.format(n))
    if outName:
        Fig.savefig('{}_Fig2_MCpicks.png'.format(outName),dpi=600)
        Fig.savefig('{}_Fig2_MCpicks.pdf'.format(outName),type='pdf')
    return Fig,ax


## Plot incremental slip rate results
def plotIncSlipRates(Rates,intervalList,analysis_method,plot_max=None,outName=None):
    # Analysis method is IQR or HPD

    # Setup figure
    F=plt.figure('IncrementalRates')
    ax=F.add_subplot(111)
    intvl_labels=[]; k=0
    # Loop through each rate
    for i in intervalList:
        # Plot full PDF
        scale_value=Rates[i].probs.max()
        ax.fill(Rates[i].rate,
            0.9*Rates[i].probs/scale_value+k,
            color=(0.4,0.4,0.4),zorder=2)
        # Plot confidence bounds
        if analysis_method=='IQR':
            ax.fill(Rates[i].xIQR,
                0.9*Rates[i].pxIQR/scale_value+k,
                color=(0.3,0.3,0.6),zorder=3)
        elif inpt.pdfAnalysis=='HPD':
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
    ax.set_yticks(np.arange(0.5,m-1))
    ax.set_yticklabels(intvl_labels,rotation='vertical')
    ax.set_ylim([0,m-1])
    if plot_max:
        ax.set_xlim([0,plot_max])
    ax.set_xlabel('slip rate')
    ax.set_title('Incremental slip rates')
    if outName:
        F.savefig('{}_Fig3_Incremental_slip_rates.png'.format(outName),dpi=600)
        F.savefig('{}_Fig3_Incremental_slip_rates.pdf'.format(outName),type='pdf')



### MAIN ---
if __name__ == '__main__':
    inpt=cmdParser()

    ## Load files
    # Read age file
    Ages,ageList,maxAge=loadAges(inpt.ageListFile,inpt.verbose,inpt.plotInputs)

    # Read disp file
    Dsps,dspList,maxDsp=loadDsps(inpt.dspListFile,inpt.verbose,inpt.plotInputs)

    # Check input files have same number of lines
    nAges=len(ageList); nDsps=len(dspList)
    assert nAges == nDsps, 'MUST HAVE SAME NUMBER OF AGE AND DISPLACEMENT MEASUREMENTS'
    m=nAges # assign number of measurements
    if inpt.verbose is True:
        print('Detected m = {} age and displacement measurements'.format(nAges))
        # Confirm pairing
        print('Pairing (youngest at top):')
        for i in range(nAges):
            print('\t{} - {}'.format(ageList[i],dspList[i]))


    ## Formulate interval names
    #  Intervals are the rates between one measurement and another
    intervalList=[]
    for i in range(m-1):
        interval_name='{}-{}'.format(dspList[i],dspList[i+1])
        intervalList.append(interval_name)


    ## Plot raw data
    plotRawData(Ages,ageList,Dsps,dspList,
        maxAge,maxDsp,
        outName=inpt.outName)


    ## Monte Carlo resamping
    AgePicks,DspPicks,RatePicks=MCMC_resample(inpt.Nsamples,
        Ages,ageList,Dsps,dspList,
        condition='standard',maxRate=inpt.maxRate,
        seed_value=inpt.seed,
        verbose=inpt.verbose,outName=None)

    # Save picks to npz file
    if inpt.outName:
        savename='{}_Picks'.format(inpt.outName)
        np.savez(savename,AgePicks=AgePicks,DspPicks=DspPicks,RatePicks=RatePicks)


    ## Plot MC results
    plotMCresults(Ages,ageList,Dsps,dspList,
        AgePicks, DspPicks,
        maxAge,maxDsp,inpt.maxPicks,
        outName=inpt.outName)

    # Save to file as script progresses
    TXTout=open('{}_slip_rate_report.txt'.format(inpt.outName),'w')

    # Raw statistics
    percentiles=[50-68.27/2, 50, 50+68.27/2]
    pct_report='''Slip rate statistics from {} samples:
Raw slip rates: (Interval: median values and 68% confidence)'''.format(inpt.Nsamples)
    print(pct_report)
    TXTout.write(pct_report)

    for i in range(m-1):
        # Calculate percentile
        pct=np.percentile(RatePicks[i,:],percentiles)
        median=pct[1]
        high_err=pct[2]-pct[1]
        low_err=pct[1]-pct[0]
        raw_stats_report='{0}: {1:.2f} +{2:.2f} -{3:.2f}'.format(intervalList[i],median,high_err,low_err)
        print(raw_stats_report); TXTout.write('\n{}'.format(raw_stats_report))


    ## Convert MC results to PDFs
    Rates={} # dictionary of object instances similar to Ages, Dsps
    for i in range(m-1):
        # Convert points to temporary (px2) array
        if inpt.pdfMethod.lower() in ['hist','histogram']:
            R=arrayHist(RatePicks[i,:],inpt.rateStep,
                smoothingKernel=inpt.smoothingKernel,
                kernelWidth=inpt.kernelWidth,
                verbose=inpt.verbose,plot=False)
        elif inpt.pdfMethod.lower() in ['kde','kernel']:
            R=arrayKDE(RatePicks[i,:],inpt.rateStep,
                smoothingKernel=inpt.smoothingKernel,
                kernelWidth=inpt.kernelWidth,
                verbose=inpt.verbose,plot=False)

        # Ascribe to object
        Rates[intervalList[i]]=rateObj(intervalList[i])
        Rates[intervalList[i]].rate=R[:,0] # rate values
        Rates[intervalList[i]].probs=R[:,1] # corresponding probabilities


    ## Analyze PDFs
    pdf_report='PDF-based statistics ({0}% confidence) calculated using {1}.'.format(inpt.rateConfidence,inpt.pdfAnalysis)
    print(pdf_report); TXTout.write('\n\n{}'.format(pdf_report))
    for i in range(m-1):

        # Interquantile range
        if inpt.pdfAnalysis in ['IQR','quantiles','quantile']:
            '''
                More stable option, but slightly skewed toward larger values
            '''
            inpt.pdfAnalysis='IQR' # confirm code for later
            PDF=IQRpdf(Rates[intervalList[i]].rate,
                Rates[intervalList[i]].probs,
                inpt.rateConfidence,verbose=False)
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
        elif inpt.pdfAnalysis in ['HPD','density']:
            '''
                More representative of probable values
            '''
            inpt.pdfAnalysis='HPD' # confirm code for later
            PDF=HPDpdf(Rates[intervalList[i]].rate,
                Rates[intervalList[i]].probs,
                inpt.rateConfidence,verbose=False)
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
    plotIncSlipRates(Rates,intervalList,inpt.pdfAnalysis,plot_max=inpt.maxRate2plot,outName=inpt.outName)


    ## Close saved file
    TXTout.close()


    if inpt.plotInputs or inpt.plotOutputs:
        plt.show()