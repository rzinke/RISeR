# Rejection sampling for Incremental Slip Rate calculation (RISeR)
This code provides a Bayesian framework for computating incremental fault slip rates, where the age and/or displacement recorded by different markers is known only to within some nonzero, finite uncertainty. The slip rate computations herein enforce the assumption that the fault did not slip backwards (inversely to its overall kinematics) at any point in its history. This type of estimation is especially important in situations for which uncertainties in age or displacement are large, or measurements overlap within uncertainty.

The incremental slip rate calcuations are carried out using one of two methods:
* Markov Chain Monte Carlo (MCMC) sampling. Inputs are sampled via the probability inverse transform method. If any of the samples drawn produces a negative slip rate, that set of samples is thrown out. Sampling continues until the desired number of samples is reached. Slip rates are reported first as percentiles of the viable picks, and then based on analysis of a pseudo-continous nonparametric function (PDF) determined from the picks. This method is ideal for data sets in which two or more observations overlap within uncertainty.
* Analytical formulation. Incremental slip rates for a dated fault slip history are computed using an analytical forumlation. Age and displacements differences between marker pairs are assumed to be always positive, and as such some increase in the accuracy and precision of the output slip rates is achieved. For complex, highly uncertain data sets, however, this functionality does not provide the full Bayesian leverage of eliminating all geologically implausible displacement-time paths, as done by the MCMC method. This analytical method is therefore ideal for data sets in which dated markers are independent (i.e., do not overlap within uncertainty).

This code builds on methods described in Gold and Cowgill, 2011 (EPSL), and is based on work developed in Zinke et al., 2017 (GRL) and Zinke et al., 2019 (GRL); and implemented in its current form in Zinke et al., 2021 (G-cubed). These codes use Python 3.6 and above, and have been tested on versions 3.6 - 3.9.
If you use these scripts, please cite:
* Zinke, R., Dolan, J.F., Rhodes, E.J., Van Dissen, R., Hatem, A.E., McGuire, C.P., Brown, N.D., Grenader, J.R. (2021) Latest Pleistocene-Holocene Incremental Slip Rates of the Wairau Fault: Implications for Long-Distance and Long-Term Coordination of Faulting Between North and South Island, New Zealand. Geochemistry, Geophysics, Geosystems, 22, e2021GC009656. https://doi.org/10.1029/2021GC009656
* Zenodo: https://zenodo.org/records/4733235. https://doi.org/10.5281/zenodo.4733235




## SETUP
The tool suite comprises a series of command-line tools, forming an adaptable pipeline for handling and computing probability functions.

These scripts have been tested for Linux/UNIX systems.
To use these scripts, clone the GitHub repository to a location of your choice, which will be referred to as ```RISeR_HOME```.
Once cloned, you must append the filepaths to your $PATH AND $PYTHONPATH variables.

If you are running BASH (```echo $SHELL``` = bash), add the following paths to your ```~/.bashrc``` file for a Linux system or ```~/.bash_profile``` for a Mac system:

```
RISeR_HOME=<path to repo>

export PATH="${PATH}:${RISeR_HOME}"
export PYTHONPATH="${PYTHONPATH}:${RISeR_HOME}/SupportFunctions"
```

If you are running CSH/TCSH (```echo $SHELL``` = (t)csh), add the following paths to your ```~/.cshrc``` or ```~/.tcshrc``` file:

```
set RISeR_HOME=<path to repo>

setenv PATH $PATH\:"$RISeR_HOME"

setenv PYTHONPATH <path to Python3>
setenv PYTHONPATH $PYTHONPATH\:"${RISeR_HOME}/SupportFunctions"
```

If you are using Windows, it is recommended you set additional environmental variables.



## EXAMPLES
Example data sets and launch scripts can be found in the ./Examples folder. An example age data set and a formatted plot are in the ExampleAges folder. An example displacement data set and corresponding plot are in the ExampleDisplacements folder. A simple example for incremental slip rate calculation (i.e., no age or displacement inversions) is found in the SimpleExample folder. A more complex example, with overlapping ages and displacement measurements is found in the ComplexExample folder.



## CONTENTS
### Executable functions
* makePDF.py - Use this to create a probability density function (PDF) as a text file. A distribution is specified, as well as the "parameters" or specific points in that distribution: For a gaussian distribution, two numbers are required 'mean stddev'; for a triangular distribution, three numbers are required 'min preferred max'; for a trapezoidal function, four numbers are required 'min shoulder1 shoulder2 max'. Supported types currently are Gaussian, tringular, and trapezoidal shapes. Gaussian distributions are printed out to 4-sigma. Uniform (boxcar) sampling distributions are specified using the trapezoidal function, with steeply sloping sides (e.g., closely approximates a boxcar). Users should run this function as a command, e.g., ```makePDF.py -d triangular -v 3 5 6.5 -o Offset_5m```

* viewPDF.py - View any PDF stored in the standard text file format, with the first column being the value of interest, and the second column being probability. For example, ```viewPDF.py Offset_5m.txt ```

* calyr2age.py - This function converts calendar ages in C.E./B.C.E. to years before "physics" (i.e., 1950 per radiocarbon convention) or years before some given calendar year (e.g., 2019). This is especially useful outputs from age calibration programs such as OxCal, which report ages in calendar years, whereas the slip rate calculator here requires age data be expressed as functions increasing in age from the present. Optional arguments include the "reference date", which defaults to 1950 C.E., and can be set to the current year, e.g., 2019. An age factor is also available, and used to scale the numbers to the desired age unit of choice (e.g., an age factor of 1000 scales the units from years to kilo-years). One may also wish to smooth their data, especially if the probability function were computed from binned values of a Monte Carlo resampling method (e.g., OxCal). Example calendar year-formatted ages can be found in the folder /Examples/ExampleAges/SampleX_Cal_Yr.txt. They can be converted into the kybp (ka) sample ages using the command ```calyr2age.py Sample1_CalYr.txt -o Sample1_age -r 2019 -f 1000 -s 3 -p```

* combinePDFs.py - If mulitple age or displacement samples are used to characterize a feature, such as two independent terrace sample ages, the combinePDFs function can be used to combine them. The user is strongly encouraged to explicitly specify the method by which the PDFs are "combined". Options include point-wise addition, and point-wise multiplication. The example sample ages in /ExampleAges/ are similar in age, and were sampled from the same stratigraphic unit, so one might be interested in the probability that a date exists within any one of those ages. To "combine" the ages by summing them, use the command ```combinePDFs.py Sample1_age.txt Sample2_age.txt -o Sample1-2_age_union -m union -p```  Alternatively, one may want to find the product (intersection) of the ages. For this, instead use ```combinePDFs.py Sample1_age.txt Sample2_age.txt -o Sample1-2_age_intersection -m intersection -p```

* betweenPDF.py - Find the probability function representing the interval between two PDFs. Suppose there are two bracketing ages each represented by a PDF: One younger than some event, and the other older than some event. This function can be used to find the probability that any given age is the "true" age of an event, i.e., older than the youngest bracketing age and younger than the oldest bracketing age. For example, ```betweenPDF.py youngest-possible_age.txt oldest-possible_age.txt -o between_age -p```

* differencePDFs.py - Compute the "delta" between two PDFs. Note that this is not a pointwise difference between two functions, that would give some similarity between the shapes of the functional forms. Rather, this function uses a modifiied convolution operation to compute the delta X_(i+1) - X_i. For example, the time length between two ages can be computed: ```differencePDFs.py OlderAge.txt YoungerAge.txt -o age_difference -p```

* dividePDFs.py - A weighted convolution function is used to compute the quotient of one PDF and another. See Bird (2007, eqns A7, A8) for derivations. For instance, the slip rate of a single slip rate marker can be computed by: ```dividePDFs.py Offset.txt Age.txt -o slipRate -p```

* plotAges.py - If the user wants to plot a series of age PDFs on a common plot, the plotAges function may be used. This requires construction of a "list" file encoded in YAML format that specifies the datum name, type of data, file location, and other parameters. An example list (AgeList.yaml) can be found in the ExampleAges folder. To generate the plot shown there, use ```plotAges.py AgeList.yaml -x 'ages (ka)' -t 'Sample Examples' -r 10 -o AgePlotExample```

* plotDisplacements.py - Similar to ```plotAges.py``` but used for displacement data. The displacement values are shown on the y-axis. For an example, try, use ```plotDisplacements.py DspList.yaml -x 'displacements (m)' -t 'Measurement Examples' -r 80 -o DspPlotExample --generic-color b```

* quickPlotSlipHistory.py - Show the dated displacement history of the fault (offset as a function of age) without computing the incremental slip rates between markers. It is recommended to use this step once the displacement and age PDFs have been established, and before slip rate computation. That way the user can get a "feel" for whether the inputs make sense in the context of this pipeline. For example, from the Examples/SimpleExample folder, try ```quickPlotSlipHistory.py DspAgeList.yaml -pt rectangle -l```

* calcSlipRates_MCMC.py - This is the main function for calculating incremental slip rates based on previously developed inputs. Users should run this function as a command, e.g., calcSlipRates.py Data.yaml. This function requires a list of dated displacement markers encoded in YAML format, as described in the INPUTS section below. This function calculates the incremental slip rates by sampling the input data (priors) and rejecting samples based on the condition that no slip rates should be negative at any point in their history. There are many optional parameters for the calcSlipRates function. Use ```calcSlipRates_MCMC.py -h``` for help.

* calcSlipRates_Analytical.py - Similar to ```calcSlipRates_MCMC.py```, this function will compute the incremental slip rates and output the results as PDFs. In this case, however, the calculations are performed using analytical formulations rather than bootstrap sampling. NOTE: This function does not yield valid results when measurements overlap within uncertainty! An error will be thrown if uncertainties in age or displacement overlap. Function syntax is similar to that of ```calcSlipRates_MCMC.py```, e.g., ```calcSlipRates_Analytical.py DspAgeData.yaml --pdf-analysis HPD```


### Support routines
* slipRateComputation.py - Provides a wrapper script to ensure consistency and proper formatting of the slip rate calculations, whether using the MCMC or analytical methods.
* slipRateObjects.py - Contains the Python classes for age, displacement, and slip rate PDFs. The age and displacement PDFs carry subroutines for computing basic statistics, and an interpolation function for inverse transform sampling. The slip rate class carries a function for converting sampled slip rate picks to a pseudo-continuous PDF.
* MCresampling.py - Used for sampling the input data and calculating the slip rates by enforcing the no-negative-rates condition. Other conditions can be specified and applied, though this is not recommended. Additionally, a maximum physically reasonable slip rate to be considered can be specified based on the user's judgement to avoid statistically implausible calculations.
* array2pdf.py - Converts an unordered array of sample picks into a continuous PDF. Two options are available. Kernel density estimation (KDE) gives a weight to data points using an automatic bandwidth determination scheme-- this tends to under weight slow slip rates and over weight fast slip rates due to inherently uneven sampling. A most reliable, alternative method is to bin the samples in a histogram using the 'hist' option. If sampling is uneven, this may lead to artificially spiky, poorly conditioned results. To overcome this limitation, smoothing methods are available. All these parameters can be specified in the calcSlipRates call.
* analyticalSlipRates.py - Computes the slip rate PDFs based on the displacement and age input PDFs.
* PDFanalysis.py - The range of possible slip rates can be quite large; it is often useful to report a range representing the most probable slip rates based on the data. To do this, two functions are provided within PDFanalysis: IQR reports the requested inter-quantile range of the PDF. This is more stable than HPD, but can be skewed, especially toward larger values. HPD reports the highest posterior density (most probable values) of a PDF. This method can give more meaningful results than IQR, but can result in anomalous values in spiky, non-smooth functions. For HPD, multiple value ranges are reported, depending on the continuity of probable values in the PDF.
* plottingFunctions - Contains functions for plotting used across modules (e.g., whisker and rectangle plots).
* dataLoading - Contains functions used across modules to handle loading of data from .yaml files, and other bookkeeping functions.
* resultSaving - Contains functions used across modules to handle formatting and saving of data.


## INPUTS
The ```calcSlipRates_XXXX.py``` routines require that the user first define a probability density function (PDF) describing the displacement and age of each marker. Each "data file" consists of a two-column list describing a single measurement of displacement or age. In each, the first column is the value measured (i.e., age or displacement), and the second column is the relative probability of that value. If a measurement can be described as parametric function, a properly formatted PDF can be generated using makePDF.py, described above. Any arbitrary, pseudo-continuous PDF can be used as an input.

Input data are specified to the ```calcSlipRates_XXXX.py``` function using a YAML (.yaml) file, which lists each dated displacement marker in order from youngest and least-offset, to oldest and most-offset. Each entry gives themarker name, followed by a dictionary-like entry specifying the path to the age PDF file, and the path to the displacement PDF file. Note that the .yaml file does not contain any data; it simply lists the feature name and files in which the data are contained.


**Note!** for OxCal outputs or any file in CE/BCE (AD/BC) format must be converted to years before present or years before physics. This can be done using the ```calyr2age.py``` function listed above.


**Disclaimer:**
This is research code. The author(s) assumes no responsibility for any errors in the methods or scripts, or any damages resulting from the code's use.
Do not use this code without citing Zinke et al., 2017; 2019. Please contact Robert Zinke with any questions or comments.