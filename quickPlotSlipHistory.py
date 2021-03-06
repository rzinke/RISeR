#!/usr/bin/env python3
'''
    ** RISeR Incremental Slip Rate Calculator **
    Quickly plot displacement and age data before calculating
     slip rates by MC sampling.
    This script relies definitions plottingFunctions.py

    Rob Zinke 2019, 2020
'''

### IMPORT MODULES ---
import argparse
from dataLoading import loadDspAgeInputs
from resultSaving import confirmOutputDir
from plottingFunctions import *


### PARSER ---
Description = '''Quickly plot displacement and age data provided as probability functions (PDFs).
The plot elements represent the 95%% confidence limits of the age and displacement represented by
each data point.'''

Examples = '''EXAMPLES
Try these with the SimpleExample data set.

# Whisker plot of raw data
quickPlotSlipHistory.py DspAgeList.yaml

# Labeled plot
quickPlotSlipHistory.py DspAgeList.yaml -l

# Rectangle plot of raw data
quickPlotSlipHistory.py DspAgeList.yaml -pt rectangle
'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    # Required
    parser.add_argument(dest='dataFile', type=str,
        help='Data file in YAML format. Each entry provides a unique marker \
name, with ageFile and dspFile specified. Youngest, least offset features at \
the top; oldest, most-offset features at the bottom.')
    parser.add_argument('-o','--out-name', dest='outName', type=str, default='Out',
        help='Head name for outputs (no extension)')
    # Recommended
    parser.add_argument('-pt','--plot-type', dest='plotType', type=str, default='whisker',
        help='Plot marker type [\'whisker\', \'rectangle\', \'PDF\']')
    parser.add_argument('-t','--title', dest='title', type=str, default='Raw data',
        help='Plot title')
    parser.add_argument('-x','--xlabel', dest='xlabel', type=str, default='age',
        help='x-axis label')
    parser.add_argument('-y','--ylabel', dest='ylabel', type=str, default='displacement',
        help='y-axis label')
    parser.add_argument('-l','--labels', dest='labelFeatures', action='store_true', default=False,
        help='Label displaced markers')
    parser.add_argument('-c','--cmap', dest='cmap', type=str, default='Greys',
        help='Cmap for headmap plots')
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true', default=False,
        help='Verbose?')
    parser.add_argument('--plot-inputs',dest='plotInputs', action='store_true',
        help='Plot inputs')
    return parser

def cmdParser(inpt_args=None):
    parser = createParser()
    return parser.parse_args(inpt_args)



### MAIN ---
if __name__ == '__main__':
    # Gather arguments
    inps = cmdParser()

    # Check output directory exists
    inps.outName = confirmOutputDir(inps.outName)

    # Load data from YAML file
    DspAgeData = loadDspAgeInputs(inps.dataFile, verbose=inps.verbose, plotInputs=inps.plotInputs)

    # Establish plot maxima
    maxAge, maxDsp = findPlotLimits(DspAgeData)

    ## Plot data
    # Plot data
    if inps.plotType.lower() in ['w', 'whisker']:
        fig, ax = plotWhiskers(DspAgeData, label=inps.labelFeatures)
    elif inps.plotType.lower() in ['r', 'rectangle']:
        fig, ax = plotRectangles(DspAgeData, label=inps.labelFeatures)
    elif inps.plotType.lower() in ['p', 'pdf']:
        fig, ax = plotJointProbs(DspAgeData, cmap=inps.cmap)


    ## Finish plot
    # X-axis
    ax.set_xlim([0, 1.1*maxAge])
    ax.set_xlabel(inps.xlabel)

    # Y-axis
    ax.set_ylim([0, 1.1*maxDsp])
    ax.set_ylabel(inps.ylabel)

    # Figure
    ax.set_title(inps.title)
    fig.tight_layout()


    ## Save if requested
    if inps.outName:
        savename = '{}.pdf'.format(inps.outName)
        fig.savefig(savename, format='pdf')


    plt.show()