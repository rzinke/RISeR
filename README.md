# Rejection sampling for Incremental Slip Rate calculation (RISeR)
This code calculates incremental fault slip rates with the assumption that the fault did not slip backwards (inversely to its overall kinematic history) at any point in its history. Inputs are given as probability density functions (PDFs) describing displacement measurements and the ages of offset markers at different epochs in a common fault history. This methodology therefore provides a Bayesian approach to incremental slip rate computation for which prior data are modified by explicit assumptions. This type of estimation is especially important in situations for which uncertainties in age or displacement are large, or measurements overlap within uncertainty.

Conditions are enforced using Markov Chain Monte Carlo (MCMC) sampling. Inputs are sampled via the probability inverse transform method. If any of the samples drawn produces a negative slip rate, that set of samples is thrown out. Sampling continues until the desired number of samples is reached. Slip rates are reported first as percentiles of the allowable picks, and then based on analysis of a continous nonparametric function determined from the picks.

This code is based on an earlier Python2 version developed in Zinke et al., 2017 (GRL) and Zinke et al., 2019 (GRL), and builds on methods described in Gold and Cowgill, 2011 (EPSL). These codes use Python 3.x, and were tested on versions 3.6, 3.7.
If you use these scripts, please cite:
* Zinke, R., Dolan, J.F., Rhodes, E.J., Van Dissen, R., McGuire, C.P. (2017) Highly Variable Latest Pleistocene-Holocene Incremental Slip Rates on the Awatere Fault at Saxton River, South Island, New Zealand, Revealed by Lidar Mapping and Luminescence Dating, Geophyical Research Letters, 44, https://doi.org/10.1002/2017GL075048
* Zinke, R., Dolan, J.F., Rhodes, E.J., Van Dissen, R., McGuire, C.P., Hatem, A.E., Brown, N.D., Langridge, R.M. (2019) Multimillennial Incremental Slip Rate Variability of the Clarence Fault at the Tophouse Road Site, Marlborough Fault System, New Zealand, Geophysical Research Letters, 46, https://doi.org/10.1029/2018GL080688



SETUP:
These scripts have been tested for Linux, Mac, and Windows systems.
To use these scripts, clone the GitHub repository to a location of your choice, which will be referred to as ```MCMC_HOME```. 
Once cloned, add the following paths to your ```~/.bash_rc``` file for a Linux system or ```~/.bash_profile``` for a Mac system:

```
export PATH="${PATH}:${MCMC_HOME}:${MCMC_HOME}/SupportFunctions"
export PYTHONPATH="${PYTHONPATH}:${MCMC_HOME}:${MCMC_HOME}/SupportFunctions"
```

If you are using Windows, it is recommended you set additional environmental variables. 



CONTENTS:
Executable functions:
* makePDF.py - Use this to create a probability density function (PDF) as a text file. A distribution is specified, as well as the "parameters" or specific points in that distribution: For a gaussian distribution, two numbers are required 'mean stddev'; for a triangular distribution, three numbers are required 'min preferred max'; for a trapezoidal function, four numbers are required 'min shoulder1 shoulder2 max'. Supported types currently are Gaussian, tringular, and trapezoidal shapes. Gaussian distributions are printed out to 4-sigma. Uniform (boxcar) sampling distributions are specified using the trapezoidal function, with steeply sloping sides (e.g., closely approximates a boxcar). Users should run this function as a command, e.g., ```makePDF.py -d triangular -v 3 5 6.5 -o Offset_5m```

* calyr2age.py - This function converts calendar ages in C.E./B.C.E. to years before "physics" (i.e., 1950 per radiocarbon convention) or years before some given calendar year (e.g., 2019). This is especially useful outputs from age calibration programs such as OxCal, which report ages in calendar years, whereas the slip rate calculator here requires age data be expressed as functions increasing in age from the present. Optional arguments include the "reference date", which defaults to 1950 C.E., and can be set to the current year, e.g., 2019. An age factor is also available, and used to scale the numbers to the desired age unit of choice (e.g., an age factor of 1000 scales the units from years to kilo-years). One may also wish to smooth their data, especially if the probability function were computed from binned values of a Monte Carlo resampling method (e.g., OxCal). Example calendar year-formatted ages can be found in the folder /Examples/ExampleAges/SampleX_Cal_Yr.txt. They can be converted into the kybp (ka) sample ages using the command ```calyr2age.py Sample1_CalYr.txt -o Sample1_age -r 2019 -f 1000 -s 3 -p``` [note the -p option to plot the result].

* combinePDFs.py - If mulitple age or displacement samples are used to characterize a feature, such as two independent terrace sample ages, the combinePDFs function can be used to combine them. The user is strongly encouraged to explicitly specify the method by which the PDFs are "combined". Options include point-wise addition, and point-wise multiplication. The example sample ages in /ExampleAges/ are similar in age, and were sampled from the same stratigraphic unit, so one might be interested in the probability that a date exists within any one of those ages. To "combine" the ages by summing them, use the command ```combinePDFs.py Sample1_age.txt Sample2_age.txt -o Sample1-2_age_union -m union -p```  Alternatively, one may want to find the product (intersection) of the ages. For this, instead use ```combinePDFs.py Sample1_age.txt Sample2_age.txt -o Sample1-2_age_union -m intersection -p```

* plotAges.py - If the user wants to plot a series of age PDFs on a common plot, the plotAges function may be used. This requires construction of a "list" file that specifies the type of data, file location, and legend entry. An example list (AgeList.txt) can be found in the ExampleAges folder. To generate the plot shown there, use ```plotAges.py AgeList.txt -x 'ages (ka)' -t 'Sample Examples' -r 10 -o AgePlotExample```

* plotDisplacements.py - Similar to ```plotAges.py``` but used for displacement data. The displacement values are shown on the y-axis. For an example, try, use ```plotDisplacements.py DspList.txt -y 'displacements (m)' -t 'Measurement Examples' -r 80 -o DspPlotExample --generic-color b```

* quickPlotSlipHistory - Show the dated displacement history of the fault (offset as a function of age) without computing the incremental slip rates between markers. It is recommended to use this step once the displacement and age PDFs have been established, and before slip rate computation. That way the user can get a "feel" for whether the inputs make sense in the context of this pipeline. For example, from the Examples/SimpleExample folder, try ```quickPlotSlipHistory.py -a AgeList.txt -d DspList.txt -pt rectangle```

* calcSlipRates.py - This is the main function for calculating incremental slip rates based on previously developed inputs. Users should run this function as a command, e.g., calcSlipRates.py -a List-of-ages.txt -d List-of-displacements.txt. This function requires two lists of filenames, as described in the INPUTS section below. This function calculates the incremental slip rates by sampling the input data (priors) and rejecting samples based on the condition that no slip rates should be negative at any point in their history. There are many optional parameters for the calcSlipRates function. Use ```calcSlipRates.py -h``` for help.


Support scripts:
* SlipRateObjects.py - Contains the classes for age and displacement PDFs, as well as an empty class for slip rate PDFs. The age and displacement PDFs carry subroutines for computing basic statistics, and an interpolation function for inverse transform sampling.
* MCresampling.py - Used for sampling the input data and calculating the slip rates by enforcing the no-negative-rates condition. Other conditions can be specified and applied, though this is not recommended. Additionally, a maximum physically reasonable slip rate to be considered can be specified based on the user's judgement to avoid statistically implausible calculations.
* array2pdf.py - Converts an unordered array of sample picks into a continuous PDF. Two options are available. Kernel density estimation (KDE) gives a weight to data points using an automatic bandwidth determination scheme-- this tends to overweight slow slip rates and overweight fast slip rates due to inherently uneven sampling. A most reliable, alternative method is to bin the samples in a histogram using the 'hist' option. If sampling is uneven, this may lead to artificially spiky, poorly conditioned results. To overcome this limitation, smoothing methods are available. All these parameters can be specified in the calcSlipRates call.
* PDFanalysis - The range of possible slip rates can be quite large; it is often useful to report a range representing the most probable slip rates based on the data. To do this, two functions are provided within PDFanalysis: IQR reports the requested inter-quantile range of the PDF. This is more stable than HPD, but can be skewed, especially toward larger values. HPD reports the highest posterior density (most probable values) of a PDF. This method can give more meaningful results than IQR, but can result in anomalous values in spiky, non-smooth functions. For HPD, multiple value ranges are reported, depending on the continuity of probable values in the PDF.


INPUTS:
The calcSlipRates function requires two text files: (1) A list of age filenames, written out from youngest at the top, to oldest at the bottom; and (2) a list of displacement names, similarly written out from yougest (least offset) at the top, to oldest (most offset) at the bottom. Note that these lists do not contain any data; they simply list the files in which the data are contained.

Each "data file" consists of a two-column list describing a single measurement of displacement or age. In each, the first column is the value measured (i.e., age or displacement), and the second column is the relative probability of that value. If a measurement can be described as parametric function, a properly formatted PDF can be generated using makePDF.py, described above. Any arbitrary, pseudo-continuous PDF can be used as an input. 

Two examples are provided. The data files are named, e.g., T1T2age.txt, T3T4dsp.txt, etc.

**Note!** for OxCal outputs or any file in CE/BCE (AD/BC) format must be converted to years before present or years before physics. This can be done using the ```calyr2age.py``` function listed above.


**Disclaimer:**
This is research code. The author(s) assumes no responsibility for any errors in the methods or scripts, or any damages resulting from the code's use.
Do not use this code without citing Zinke et al., 2017; 2019. Please contact Robert Zinke with any questions or comments.
