#!/usr/bin/env python3
import os
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
from slipRateObjects import gauss_kernel

'''
This script accepts two probability density functions and computes the
 normalized cross-correlation value of them.
Definitions within this script are intended for generic use, and so 
 some operations might be redundant. 
'''

### --- PARSER --- ###
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='This script computes the normalized cross-correlation for two PDFs')
	parser.add_argument(dest='filenames',nargs='+',help='Functions, stored as files formatted with col 0 = x-values; col 1 = y-values')
	parser.add_argument('-s','--smoothing',dest='smoothing',type=int,default=None,help='Smoothing filter width. Default is None')
	parser.add_argument('-interp_kind','--interp_kind',dest='interp_kind',type=str,default='linear',help='Kind of interpolation. Default is \'linear\'')
	parser.add_argument('-p','--plot',dest='plot',action='store_true',help='Plot inputs')
	return parser

def cmdParser(iargs = None):
	parser = createParser()
	return parser.parse_args(args=iargs)



### --- PROCESSING FUNCTIONS --- ###
# Normalize 1D function
def normFct(y):
	y/=np.sqrt(np.sum(y**2))
	return y 

# Normalized 1D cross-correlation
def normd_cross(y1,y2):
	# Reshape into vectors
	y1=y1.reshape(1,-1)
	y2=y2.reshape(1,-1)

	# Normalize
	Y1=np.sqrt(np.sum(y1**2)) 
	Y2=np.sqrt(np.sum(y2**2))

	# Cross-correlate
	Y=np.dot(y1,y2.T)/(Y1*Y2)
	return Y 

# Print triangular matrix of values
def print_table(X):
	m,n=X.shape # rows and columns
	for i in range(1,m):
		X[i,:i]=None
	print(X)


### --- MAIN --- ###
if __name__=='__main__':
	inpt=cmdParser()

	## Load data
	# Set up parameters
	N=len(inpt.filenames) # number of functions
	Fns=[]; xmin=[]; xmax=[]; dx=[] # empty variables
	for fname in inpt.filenames:
		data=np.loadtxt(fname)
		Fns.append(data)

		# Record function parameters for global function
		xmin.append(data[:,0].min())
		xmax.append(data[:,0].max())
		dx.append(np.mean(np.diff(data[:,0])))

	# Global function parameters
	xmin=np.min(xmin)
	xmax=np.max(xmax)
	dx=np.mean(dx)

	## Format functions for analysis
	# Define common range of values
	xresampled=np.arange(xmin,xmax+dx,dx).reshape(-1,1) 
	nx=len(xresampled) # number of data points in xresampled

	# Format each function for analysis
	ResampledFns=[] # empty list for resampled functions
	for i in range(N):
		x=Fns[i][:,0] # x-values
		y=Fns[i][:,1] # y-values

		# Smooth if requested
		if inpt.smoothing:
			kernel=gauss_kernel(inpt.smoothing)
			y=np.convolve(y,kernel,'same')
		
		# Interpolate
		I=interp1d(x,y,kind=inpt.interp_kind,bounds_error=False,fill_value=0.)

		# Resample on common axis
		yresampled=I(xresampled)

		# Normalize values
		ynorm=np.sqrt(np.sum(yresampled**2))
		yresampled/=ynorm

		resampled_data=np.hstack([xresampled,yresampled])
		ResampledFns.append(resampled_data)


	## Compute cross-product
	# For only 2 functions
	if N==2:
		C=np.dot(ResampledFns[0][:,1],ResampledFns[1][:,1])
		print('Correlation: {:.3f}'.format(C))

	# For more than two functions
	elif N>2:
		# Use covariance matrix
		X=np.zeros((N,nx))
		for i in range(N):
			X[i,:]=ResampledFns[i][:,1].reshape(1,nx)
		C=X.dot(X.T) # covariance matrix

		# Print values in table
		print_table(C)

	else:
		print('Must specify more than one function to correlate')

	## Plot if requested
	if inpt.plot is True:
		F=plt.figure()
		ax=F.add_subplot(111)
		# Plot resampled functions
		for i in range(N):
			resampled_data=ResampledFns[i]
			label=os.path.basename(inpt.filenames[i])
			ax.plot(resampled_data[:,0],resampled_data[:,1],label=label)
		ax.legend(loc='best')
		ax.set_title('Resampled and normalized data')



	if inpt.plot is True:
		plt.show()