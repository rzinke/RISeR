#!/bin/bash

##########
# Create an incremental slip rate calculations from scratch
#  using RISeR routines.
# R Zinke 19 Aug 2020
##########


## CREATE INPUT DATA
# Make age PDFs
#  Specify distribution to model ages and give parameters in kyrs
makePDF.py -d gauss -v 4.0 0.8 -o T3T4age

makePDF.py -d gauss -v 6.0 1.0 -o T2T3age

makePDF.py -d gauss -v 10.0 1.0 -o T1T2age


# Make displacement PDFs
#  Specify distribution and give parameters in meters
makePDF.py -d tri -v 4.0 6.0 9.0 -o T3T4offset

makePDF.py -d trap -v 7.0 8.0 10.0 11.0 -o T2T3offset

makePDF.py -d tri -v 12.0 15.0 16.0 -o T1T2offset


## CREATE INPUT FILE
# File for specifying inputs
cat > Inputs.yaml << EOM
T3T4riser: {"ageFile": "T3T4age.txt", "dspFile": "T3T4offset.txt"}
T2T3riser: {"ageFile": "T2T3age.txt", "dspFile": "T2T3offset.txt"}
T1T2riser: {"ageFile": "T1T2age.txt", "dspFile": "T1T2offset.txt"}
EOM


## CALCULATE INCREMENTAL SLIP RATES
Outname=QuickExample  # output head name
Nsamples=100000   # number of valid samples
maxRate=100       # maximum slip rate to be considered (e.g., mm/yr)
PDFmethod=hist    # method to convert sample picks into continuous function
PDFanalysis=HPD   # report PDF range of values
SlipRateStep=0.1  # step between slip rate PDF x-axis values
SmoothType=gauss  # smoothing recommended for histogram of picks
SmoothWidth=3     # how strong is smoothing
Confidence=68.27  # percent confidence to report
PlotMax=10        # maximum of incr. slip rate plot

calcSlipRates.py Inputs.yaml -o $Outname -n $Nsamples \
--pdf-method $PDFmethod --pdf-analysis $PDFanalysis \
--max-rate $maxRate --rate-step $SlipRateStep \
--smoothing-kernel $SmoothType --kernel-width $SmoothWidth \
--rate-confidence $Confidence --max-rate2plot $PlotMax \
--plot-outputs --verbose
