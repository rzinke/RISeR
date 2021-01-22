'''
** RISeR Incremental Slip Rate Calculator **
Generic wrapper to ensure consistency in incremental slip rate formatting and
 results.

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import matplotlib.pyplot as plt
from dataLoading import loadDspAgeInputs
from resultSaving import confirmOutputDir, slipRateTxtFile, printIncSlipRates
from plottingFunctions import plotRawData, plotIncSlipRates



### GENERIC WRAPPER FOR SLIP RATE COMPUTATION ---
def computeIncrRates(scheme, args):
    '''
    This function acts to ensure consistency in slip rate formatting and output
     structure between the two methods for computation: MCMC, and analytical
     formulation.
    INPUTS
        scheme is the type of method to be used (MCMC, analytical)
        args come directly from argparse, and is an object comprising all
         input arguments
    '''
    ## Setup
    # Determine workflow
    scheme = scheme.lower()  # format

    # Check output directory exists
    outName = confirmOutputDir(args.outName, verbose=args.verbose)

    # Start text file
    txtFile = slipRateTxtFile(args.outName, scheme)


    ## Load input data
    # Load data from YAML file
    DspAgeData = loadDspAgeInputs(args.dataFile, 
        verbose=args.verbose, printDetails=args.xtrVerbose, plotInputs=args.plotInputs)

    # Plot raw data
    figNb = 1  # start figure counter
    _, _, figNb = plotRawData(DspAgeData, figNb, label=args.labelMarkers,
        ageUnits=args.ageUnits, dspUnits=args.dspUnits, outName=args.outName)


    ## Compute incremental slip rates
    if scheme in ['analytical']:
        # Compute rates using analytical wrapper function below
        Rates = computeAnalyticalRates(DspAgeData, args, txtFile)

    elif scheme in ['mcmc']:
        # Compute rates using MC wrapper function below
        Rates = computeMCMCrates(DspAgeData, args, txtFile)

        # Update figure counter
        figNb += 1


    ## Analyze PDFs
    # Compute statistics based on PDFs
    analyzePDFs(Rates, method=args.pdfAnalysis, confidence=args.rateConfidence,
        verbose=args.verbose, outName=args.outName)

    # Plot incremental slip rates on same plot
    plotIncSlipRates(Rates, args.pdfAnalysis, figNb, plotMax=args.maxRate2plot, 
        rateUnits=args.rateUnits, outName=args.outName)

    # Print intervals to file
    printIncSlipRates(txtFile, Rates, args.pdfAnalysis, args.rateConfidence, rateUnits=args.rateUnits)

    if args.plotOutputs == True or args.plotInputs == True:
        plt.show()



### ANALYTICAL SLIP RATES ---
def computeAnalyticalRates(DspAgeData, args, txtFile):
    '''
    Wrapper function to compute slip rates using analytical methods.
    '''
    # Import appropriate modules
    from analyticalSlipRates import analyticalSlipRates

    # Use analytical methods
    Rates = analyticalSlipRates(DspAgeData, 
        stepSize=args.stepSize, maxRate=args.maxRate,
        verbose=args.verbose, printDetails=args.xtrVerbose)

    return Rates



### MONTE CARLO SLIP RATES ---
def computeMCMCrates(DspAgeData, args, txtFile):
    '''
    Wrapper function to compute slip rates using Monte Carlo methods.
    '''
    # Import appropriate modules
    from MCresampling import MCMCresample, picks2PDF, rawPercentiles
    from plottingFunctions import plotMCresults

    # Use Markov chain Monte Carlo method
    AgePicks, DspPicks, RatePicks = MCMCresample(DspAgeData,
        Nsamples=args.Nsamples, condition='standard',
        maxRate=args.maxRate, bound=args.MCbound,
        seedValue=args.seed,
        verbose=args.verbose,
        outName=args.outName)

    # Compute raw percentiles
    rawPercentiles(DspAgeData, RatePicks, txtFile, args.rateConfidence, verbose=args.verbose)

    # Plot MC results
    figNb = 2
    plotMCresults(DspAgeData, AgePicks, DspPicks, figNb, maxPicks=args.maxPicks2plot,
        ageUnits=args.ageUnits, dspUnits=args.dspUnits, outName=args.outName)

    # Convert MC results to PDFs
    Rates = picks2PDF(DspAgeData, RatePicks, 
        method=args.pdfMethod,
        stepSize=args.rateStep,
        smoothingKernel=args.smoothingKernel, kernelWidth=args.kernelWidth,
        verbose=args.verbose)

    return Rates



### PDF ANALYSIS ---
def analyzePDFs(Rates, method, confidence, verbose=False, outName=None):
    '''
    Loop through slip rate PDFs to retrieve values using interquantile range
     (IQR) or highest posterior density (HPD).
    '''
    if verbose == True:
        print('*'*32)
        print('Extracting slip rate values from PDFs using the {:s} method'.format(method))

    for intvl in Rates.keys():
        Rates[intvl].analyzePDF(method=method, confidence=confidence)

        # Save rate PDF to file if specified
        if outName:
            pdfName = '{:s}_Incr_slip_rate_{:s}'.format(outName, intvl.replace('/', ''))
            Rates[intvl].save2txt(pdfName)