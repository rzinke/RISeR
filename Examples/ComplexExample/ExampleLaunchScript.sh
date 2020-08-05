#!/bin/bash
# This script will compute the incremental slip rates between four dated
#  displacement markers for a hypothetical data set.
#  10,000 valid, random samples will be drawn from the PDFs described in
#  AgeList.txt and DspList.txt.
#
# Unlike the "Simple" example, dated displacement markers in this script overlap
#  in both offset and age. The incremental slip rate computations applied here
#  account for that overlap by rejecting sample combinations that give negative
#  slip rates.
#
# Because of the large uncertainties, I have increased the number of sample
#  draws from the usual 10,000, to 50,000. Hypothetically, the maximum allowable
#  slip rate from these data is infinite, so I manually set the --max-rate to
#  50 mm/yr. I used a histogram method for constructing the slip rate PDFs from
#  the random picks, which I smoothed using a smoothing kernel with a width of
#  5 smaples (i.e., 2-sigma width of ~0.5 mm/yr).

calcSlipRates.py Dsp-AgeList.yaml -o outputs/Practice -n 10 \
    --max-rate 50 --pdf-method hist --rate-step 0.1 --smoothing-kernel gauss \
    --kernel-width 5 --pdf-analysis HPD \
    -v -p -l --plot-inputs
