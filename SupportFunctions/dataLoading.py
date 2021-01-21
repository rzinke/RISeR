#!/usr/bin/env python3
'''
** MCMC Incremental Slip Rate Calculator **
These functions are designed to load and parse properly
 formatted data.

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import sys
import os
try:
    import yaml
except:
    print('Please install pyyaml'); exit()
from slipRateObjects import ageDatum, dspDatum


### ANCILLARY FUNCTIONS ---
## Check Python version
def checkVersion():
    '''
    Confirm Python v3.6+ is used due to the necessity of orderered
     dictionary keys.
    '''
    Vs = sys.version_info
    VsMajor = Vs[0]
    VsMinor = Vs[1]
    Vs = float('{}.{}'.format(VsMajor, VsMinor))
    if Vs < 3.6:
        print('Python version must be 3.6+ to use loadInputs function due to \
the necessity of ordered dictionary keys. Please upgrade to v. 3.6 or higher.')
        exit()


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

    # Return outName as absolute path
    return outName



### LOADING FUNCTIONS ---
## Load displacement-age inputs from YAML file for slip rate analysis
def loadDspAgeInputs(fname, verbose=False, plotInputs=False):
    '''
    Load age and displacement data based on YAML inputs.
    Inputs should be specified as one input per line.
     The name of the data point is given, followed by a dictionary-like
     entry specifying the age file and displacement file. E.g.,
      T1/T2 riser: {"ageFile": "T2age.txt", "dspFile": "T1T2dsp.txt"}
    Returns a dictionary of displacement-age markers, where each marker has
     an associated age PDF ("Age") and displacement PDF ("Dsp"). The PDFs
     are loaded as ageDatum and dspDatum objects, respectively.

    NOTE: This function should only be used with Python v.s 3.6 or higher
     due to the necessity of ordered dictionary keys.
    '''
    # Check Python version
    checkVersion()

    # Report if requested
    if verbose == True: print('Loading displacement-age data')

    # Open .yaml file
    with open(fname, 'r') as inputFile:
        # Parse data within file
        DspAgeData = yaml.load(inputFile, Loader=yaml.FullLoader)

        # Check there are two or more markers given
        if len(DspAgeData) < 2:
            print('More than one marker must be specified for incremental slip rate calculation.')
            exit()

        # Loop through each displacement-age datum to load data
        for datumName in DspAgeData.keys():
            # Dictionary entry
            datum = DspAgeData[datumName]

            ageFile = datum['ageFile']
            ageName = os.path.basename(ageFile).split('.')[0]

            dspFile = datum['dspFile']
            dspName = os.path.basename(dspFile).split('.')[0]

            # Report if requested
            if verbose == True:
                print('*******************************')
                print('Datum name: {:s}'.format(datumName))
                print('\tAge file: {:s}'.format(os.path.basename(ageFile)))
                print('\tDisp file: {:s}'.format(os.path.basename(dspFile)))

            # Load age PDF
            datum['Age'] = ageDatum(name = ageName)
            datum['Age'].readFromFile(filepath = ageFile)
            datum['Age'].format(verbose = verbose)
            if plotInputs == True: datum['Age'].plot()

            # Load dsp PDF
            datum['Dsp'] = dspDatum(name = dspName)
            datum['Dsp'].readFromFile(filepath = dspFile)
            datum['Dsp'].format(verbose = verbose)
            if plotInputs == True: datum['Dsp'].plot()

    return DspAgeData