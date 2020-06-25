#!/bin/bash

python3 ../calcSlipRates.py -a AgeList.txt -d DspList.txt -o Practice -n 20000 -max_rate 60 -pdf_method hist -rate_step 0.1 -smoothing_kernel gauss -kernel_width 5 -pdf_analysis HPD --verbose --plot_outputs
