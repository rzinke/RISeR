#!/usr/bin/env python3
'''
** RISeR Incremental Slip Rate Calculator **
This script serves as a wrapper for computing incremental fault slip rates using analytical formulations.

It structures the input arguments and passes them to a generic wrapper to ensure consistency in slip rate formatting
 and results.

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import argparse
from slipRateComputation import computeIncrRates


### PARSER ---
Description = '''Calculate incremental slip rates for a dated fault slip history using an analytical forumlation. 
This method is ideal for data sets in which dated markers are independent (i.e., do not overlap within
uncertainty).

Inputs are provided as probability density functions (PDFs) representing the age and displacement of each marker. From 
those PDFs, random samples are drawn, and samples resulting in negative slip rates are discarded. From the n valid 
random draws, incremental slip rates are calculated.

REQUIRED INPUT
This routine requires a YAML file (e.g., data.yaml) which lists each dated displacement marker in order from youngest
and least-offset, to oldest and most-offset. Each entry gives the marker name, followed by a dictionary-like entry
specifying the path to the age PDF file, and the path to the displacement PDF file.

For example: The file data.yaml might contain
# Commented description of data scheme
T2/T3 riser: {\"ageFile\": \"T2T3age.txt\", \"dspFile\": \"T2T3dsp.txt\"}
T1/T2 riser: {\"ageFile\": \"T1T2age.txt\", \"dspFile\": \"T1T2dsp.txt\"}

Where T2/T3 is younger and less-offset than T1/T2. Change the .txt filenames to the relativeor absolute paths to the PDF
files, accordingly.

Note: More than one entry must be present to calculate incremental slip rates.
'''

Examples = '''EXAMPLES
From the Examples/SimpleExample folder

calcSlipRates_Analytical.py DspAgeData.yaml
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
    requiredArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out',
        help='Head name for outputs (no extension). [Default = \'Out\'].')

    # Generic arguments
    generalArgs = parser.add_argument_group('GENERIC ARGUMENTS')
    generalArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', default=False,
        help='Verbose mode.')
    generalArgs.add_argument('-vv','--extra-verbose', dest='xtrVerbose', action='store_true',
        help='Print all possible statements to screen.')
    generalArgs.add_argument('-p','--plot-outputs', dest='plotOutputs', action='store_true',
        help='Plot outputs.')
    generalArgs.add_argument('-l','--label-markers', dest='labelMarkers', action='store_true',
        help='Label dated displacement markers on raw data plot.')

    # Fine-tuning
    detailArgs = parser.add_argument_group('DETAILED MC ARGUMENTS')
    detailArgs.add_argument('--max-rate', dest='maxRate', type=float, default=1E2,
        help='Maximum rate considered in analysis. Units are <dispalcement units> per <age units>. [Default = 100].')
    detailArgs.add_argument('--step-size', dest='stepSize', type=float, default=1E-2,
        help='Step size of quotient axis, in units of <numerator units>/<denominator units>. [Default 0.01].')

    detailAnalysisArgs = parser.add_argument_group('DETAILED SLIP RATE ANALYSIS ARGUMENTS')
    detailAnalysisArgs.add_argument('--pdf-method', dest='pdfMethod', type=str, default='kde',
        help='Method used for transforming rate picks into slip rate PDF (\'hist\', [\'kde\']).')
    detailAnalysisArgs.add_argument('--pdf-analysis', dest='pdfAnalysis', type=str, default='IQR',
        help='Method for analyzing slip rate PDFs. ([\'IQR\'] for interquantile range; \'HPD\' for highest posterior density).')
    detailAnalysisArgs.add_argument('--rate-confidence', dest='rateConfidence', type=float, default=68.27,
        help='Confidence range for slip rate PDF reporting, (percent, e.g., [68.27], 95.45, etc.).')

    detailFigureArgs = parser.add_argument_group('DETAILED SLIP RATE PLOT ARGUMENTS')
    detailFigureArgs.add_argument('--max-rate2plot', dest='maxRate2plot', type=float, default=40,
        help='Maximum spreading rate to plot (unlike -max-rate, this will not affect calculations). [Default = None].')
    detailFigureArgs.add_argument('--plot-inputs', dest='plotInputs', action='store_true',
        help='Plot inputs.')
    return parser

def cmdParser(inps_args=None):
    parser = createParser()
    return parser.parse_args(inps_args)



### MAIN ---
if __name__ == '__main__':
    ## Gather inputs
    inps = cmdParser()


    ## Pass arguments to slip rate wrapper
    computeIncrRates('analytical', inps)