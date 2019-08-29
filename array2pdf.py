#!/usr/bin/env python3
import numpy as np 
import matplotlib.pyplot as plt 
from slipRateObjects import gauss

'''
Covert an array of numbers into a probability density function (PDF).
'''

## Histogram method
# Convert to PDF using histogram (more stable over fast intervals)
def arrayHist(V,stepsize,smoothing_kernel=None,kernel_width=2,verbose=False,plot=False):
	'''
	INPUTS: 
		V is an array of values
	OUTPUTS:
		x is the values
		px is the probability of occurrence
	'''
	if verbose is True:
		print('Converting array to histogram')
		# Compute raw statistics
		pct=np.percentile(V,[50-95.45/2,50-68.27/2,50,50+68.27/2,50+95.45/2])
		print('Raw statistics:')
		print('\tmedian: {0:.3f}'.format(pct[2]))
		print('\t68.27% range: {0:.3f}-{1:.3f}'.format(pct[1],pct[3]))
		print('\t95.45% range: {0:.3f}-{1:.3f}'.format(pct[0],pct[4]))
		print('\tmin: {0:.3f}; max: {1:.3f}'.format(V.min(),V.max()))

	# Bin edges/independent axis
	bins=np.arange(V.min(),V.max()+stepsize,stepsize)

	# Compute histogram
	H,Hedges=np.histogram(V,bins=bins)
	Hcntrs=(Hedges[:-1]+Hedges[1:])/2 # centers of bins

	# Taper histogram edges
	Hcntrs=np.pad(Hcntrs,(1,1),'constant',constant_values=(Hedges[0],Hedges[-1]))
	H=np.pad(H,(1,1),'constant')

	# Smooth if requested
	if smoothing_kernel:
		if verbose is True:
			print('Applying smoothing kernel:\n\ttype: {}\twidth: {}'.format(smoothing_kernel,kernel_width))
		# Moving mean smoothing
		if smoothing_kernel.lower() in ['mean']:
			# Form kernel
			K=np.ones(kernel_width)
			# Apply via convolution
			H=np.convolve(H,K,'same')
			# Set ends back to zero
			H[0]=0; H[-1]=0
		elif smoothing_kernel.lower() in ['gauss','gaussian']:
			# Form kernel
			kernel_width=kernel_width/2
			x=np.linspace(-2*kernel_width,2*kernel_width,6*kernel_width) # width of kernel = +- 2 std dev
			K=gauss(x,0,kernel_width)
			# Apply via convolution
			H=np.convolve(H,K,'same')
			# Set ends back to zero
			H[0]=0; H[-1]=0

	# Normalize area to 1.0

	# Plot if requested
	if plot is True:
		F=plt.figure()
		ax=F.add_subplot(111)
		ax.plot(Hcntrs,H,'k',linewidth=1)

	#return x, px


## KDE method
# Convert to PDF using kernel density estimation (inherently smoother)
def arrayKDE(V,smoothing=None,verbose=False):
	'''
	INPUTS: 
		V is an array of values
	OUTPUTS:
		x is the values
		px is the probability of occurrence
	'''
	if verbose is True:
		print('Converting array to kernel density plot')

	# Compute raw statistics

	# Normalize area to 1.0

	# return x, px