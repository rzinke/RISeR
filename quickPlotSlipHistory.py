#!/usr/bin/env python3
'''
    ** MCMC Incremental Slip Rate Calculator **
    Quickly plot displacement and age data before calculating
     slip rates by MC sampling.
    This script relies heavily on functions from calcSlipRates.py

    Rob Zinke 2019, 2020
'''

### IMPORT MODULES ---
import argparse
from calcSlipRates import *


### PARSER ---
Description = '''Quickly plot displacement and age data provided as probability functions (PDFs).
The plot elements represent the 95%% confidence limits of the age and displacement represented by
each data point.'''

Examples = '''EXAMPLES
Try these with the SimpleExample data set.

# Whisker plot of raw data
quickPlotSlipHistory.py -a AgeList.txt -d DspList.txt

# Labeled plot
quickPlotSlipHistory.py -a AgeList.txt -d DspList.txt -l

# Rectangle plot of raw data
quickPlotSlipHistory.py -a AgeList.txt -d DspList.txt -pt rectangle
'''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    # Required
    parser.add_argument('-a','--age-list', dest='ageListFile', type=str, required=True,
        help='Text file with one age file per line, list in order from youngest (top line) to oldest (bottom).')
    parser.add_argument('-d','--dsp-list', dest='dspListFile', type=str, required=True,
        help='Text file with one displacement file per line, list in order from smallest (top line) to largest (bottom).')
    parser.add_argument('-o','--out-name', dest='outName', type=str, default='Out',
        help='Head name for outputs (no extension)')
    # Recommended
    parser.add_argument('-pt','--plot-type', dest='plotType', type=str, default='whisker',
        help='Plot marker type [\'wkisker\',\'rectangle\']')
    parser.add_argument('-t','--title', dest='title', type=str, default='Raw data',
        help='Plot title')
    parser.add_argument('-x','--xlabel', dest='xlabel', type=str, default='age',
        help='x-axis label')
    parser.add_argument('-y','--ylabel', dest='ylabel', type=str, default='displacement',
        help='y-axis label')
    parser.add_argument('-l','--labels', dest='labelFeatures', action='store_true', default=False,
        help='Label displaced markers')
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true', default=False,
        help='Verbose?')
    parser.add_argument('--plot-inputs',dest='plotInputs', action='store_true',
        help='Plot inputs')
    return parser

def cmdParser(inpt_args=None):
    parser = createParser()
    return parser.parse_args(inpt_args)



### PLOTTING FUNCTIONS ---
def whiskerPlot(ax,Ages,ageList,Dsps,dspList,labels=False):
    '''
        Plot uncertainties as error bars.
    '''
    m = len(ageList)
    for i in range(m):
        # Format age data
        x_mid = Ages[ageList[i]].median
        x_err = np.array([[x_mid-Ages[ageList[i]].lowerLimit],
            [Ages[ageList[i]].upperLimit-x_mid]])
        # Format disp data
        y_mid = Dsps[dspList[i]].median
        y_err = np.array([[y_mid-Dsps[dspList[i]].lowerLimit],
            [Dsps[dspList[i]].upperLimit-y_mid]])
        # Plot data
        ax.errorbar(x_mid, y_mid,
            xerr=x_err, yerr=y_err,
            color=(0.3,0.3,0.6), marker='o')
        # Labels
        if labels == True:
            ax.text(x_mid,y_mid,dspList[i])


def rectanglePlot(ax,Ages,ageList,Dsps,dspList,labels=False):
    '''
        Plot uncertainties as rectangles.
    '''
    m = len(ageList)
    for i in range(m):
        # Rectangle parameters
        ageLower = Ages[ageList[i]].lowerLimit # load bottom
        dspLower = Dsps[dspList[i]].lowerLimit # load left
        boxWidth = Ages[ageList[i]].upperLimit-ageLower
        boxHeight = Dsps[dspList[i]].upperLimit-dspLower
        # Plot data
        ax.add_patch(Rectangle((ageLower,dspLower), # LL corner
            boxWidth,boxHeight, # dimensions
            edgecolor=(0.3,0.3,0.6), fill=False, zorder=3))
        # Labels
        if labels == True:
            ax.text(ageLower+boxWidth,dspLower,dspList[i])


### MAIN ---
if __name__ == '__main__':
    inps = cmdParser()

    ## Load files
    # Read age file
    Ages,ageList,maxAge=loadAges(inps.ageListFile,inps.verbose,inps.plotInputs)

    # Read disp file
    Dsps,dspList,maxDsp=loadDsps(inps.dspListFile,inps.verbose,inps.plotInputs)

    # Check input files have same number of lines
    checkInputs(ageList,dspList,verbose = inps.verbose)

    # Formulate interval names
    intervalList = intervalsFromDspList(dspList)


    ## Plot data
    # Establish figure
    Fig = plt.figure()
    ax = Fig.add_subplot(111)

    # Plot data
    if inps.plotType == 'whisker':
        whiskerPlot(ax,Ages,ageList,Dsps,dspList,
            labels = inps.labelFeatures)
    elif inps.plotType == 'rectangle':
        rectanglePlot(ax,Ages,ageList,Dsps,dspList,
            labels = inps.labelFeatures)


    ## Finish plot
    # X-axis
    ax.set_xlim([0, 1.1*maxAge])
    ax.set_xlabel(inps.xlabel)

    # Y-axis
    ax.set_ylim([0, 1.1*maxDsp])
    ax.set_ylabel(inps.ylabel)

    # Figure
    ax.set_title(inps.title)
    Fig.tight_layout()


    ## Save if requested
    if inps.outName:
        Fig.savefig(inps.outName)


    plt.show()