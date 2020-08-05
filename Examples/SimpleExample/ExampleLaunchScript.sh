#!/bin/bash
# This script will compute the incremental slip rates between three dated
#  displacement markers for a hypothetical data set. 10,000 valid, random
#  samples will be drawn from the PDFs described in AgeList.txt and DspList.txt.
#  Options for verbose mode and live plotting are selected.

calcSlipRates.py DspAgeList.yaml -o outputs/run1 -n 100000 --pdf-analysis HPD -v -p -l