#!/bin/bash

python3 ../calcSlipRates.py -a AgeExs.txt -d DspExs.txt -o Practice -n 10000 -plot_outputs True -max_rate 60 -pdf_method hist -rate_step 0.1 -smoothing_kernel gauss -kernel_width 5 -pdf_analysis HPD
