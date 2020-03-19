"""
	** MCMC Incremental Slip Rate Calculator **
	Objects and support functions for MCMC slip rate calculator

	Rob Zinke 2019, 2020
"""

### IMPORT MODULES ---
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import cumtrapz
from PDFanalysis import *



### AGE, DISPLACEMENT, SLIP RATE CLASSES ---
## Age datum class
class ageDatum:
	"""
		PDF representing an age measurement. Units are thousands of years.
	"""
	def __init__(self,name):
		# Basic parameters
		self.name=name

	# Read in data from csv or similar text file
	def readFromFile(self,filepath):
		# Format raw data
		data=np.loadtxt(filepath)
		self.ages=data[:,0]
		self.probs=data[:,1]

	# Format for use in slip rate analysis
	def format(self,verbose=False,plot=False):
		if verbose is True:
			print('Formatting {} for slip rate analysis'.format(self.name))

		# Sum to CDF
		P=np.trapz(self.probs,self.ages) # area under curve
		self.probs/=P # normalize area to 1.0
		cdf_all=cumtrapz(self.probs,self.ages,initial=0)

		# Use only unique values
		uniqueVals,uniqueNdx=np.unique(cdf_all,return_index=True)
		self.ages=self.ages[uniqueNdx] # unique ages
		self.probs=self.probs[uniqueNdx] # unique age probs
		self.probs=self.probs/np.trapz(self.probs,self.ages) # re-norm area
		self.cdf=cumtrapz(self.probs,self.ages,initial=0) # re-calc CDF
		if verbose is True:
			print('\t...limiting to unique values')
			print('\t...final CDF value: {}'.format(self.cdf[-1]))

		# Inverse interpolation function
		#	use cdf as "x" value for inverse interpolation
		#	leave kind as linear to avoid values < 0 or > 1
		self.InvCDF=interp1d(self.cdf,self.ages,kind='linear')
		if verbose is True:
			print('\t...built inverse interpolation function')

		# Basic statistics
		self.lowerLimit,self.median,self.upperLimit=self.InvCDF([0.025,0.5,0.975])
		if verbose is True:
			print('\t95.45 % limits: {0:.5f}; {1:.5f}'.format(self.lowerLimit,self.upperLimit))

		# Plot if requested
		if plot is True:
			F=plt.figure()
			ax=F.add_subplot(111)
			ax.plot(self.ages,0.8*self.probs/self.probs.max(),color=(0.3,0.3,0.6),label='PDF')
			ax.plot(self.ages,self.cdf,color='k',linewidth=2,label='CDF')
			ax.plot([self.lowerLimit,self.upperLimit],[0,0],'rx',label='95% bounds')
			ax.set_ylim([-0.1,1.1])
			ax.set_xlabel('age'); ax.set_ylabel('prob')
			ax.set_title('INPUT: {}'.format(self.name))
			ax.legend()


## Displacement datum class
class dspDatum:
	"""
		PDF representing a displacement measurement. Units are meters.
	"""
	def __init__(self,name):
		# Basic parameters
		self.name=name

	# Read in data from csv or similar text file
	def readFromFile(self,filepath):
		# Format raw data
		data=np.loadtxt(filepath)
		self.dsps=data[:,0] 
		self.probs=data[:,1]

	# Format for use in slip rate analysis
	def format(self,verbose=False,plot=False):
		if verbose is True:
			print('Formatting {} for slip rate analysis'.format(self.name))

		# Sum to CDF
		P=np.trapz(self.probs,self.dsps) # area under curve
		self.probs/=P # normalize area to 1.0
		cdf_all=cumtrapz(self.probs,self.dsps,initial=0)

		# Use only unique values
		uniqueVals,uniqueNdx=np.unique(cdf_all,return_index=True)
		self.dsps=self.dsps[uniqueNdx] # unique displacements
		self.probs=self.probs[uniqueNdx] # unique displacement probs
		self.probs=self.probs/np.trapz(self.probs,self.dsps) # re-norm area
		self.cdf=cumtrapz(self.probs,self.dsps,initial=0) # re-calc CDF
		if verbose is True:
			print('\t...limiting to unique values')
			print('\t...final CDF value: {}'.format(self.cdf[-1]))

		# Inverse interpolation function
		#	use cdf as "x" value for inverse interpolation
		#	leave kind as linear to avoid values < 0 or > 1
		self.InvCDF=interp1d(self.cdf,self.dsps,kind='linear')
		if verbose is True:
			print('\t...built inverse interpolation function')

		# Basic stats
		self.lowerLimit,self.median,self.upperLimit=self.InvCDF([0.025,0.5,0.975])
		if verbose is True:
			print('\t95.45 % limits: {0:.5f}; {1:.5f}'.format(self.lowerLimit,self.upperLimit))

		# Plot if requested
		if plot is True:
			F=plt.figure()
			ax=F.add_subplot(111)
			ax.plot(self.dsps,0.8*self.probs/self.probs.max(),color=(0.3,0.3,0.6),label='PDF')
			ax.plot(self.dsps,self.cdf,color='k',linewidth=2,label='CDF')
			ax.plot([self.lowerLimit,self.upperLimit],[0,0],'rx',label='95% bounds')
			ax.set_ylim([-0.1,1.1])
			ax.set_xlabel('displacement'); ax.set_ylabel('prob')
			ax.set_title('INPUT: {}'.format(self.name))
			ax.legend()


## Slip rate object
class rateObj:
	def __init__(self,name):
		self.name=name
		# Add parameters as needed



### MISCELLANEOUS ---
## Definition of gaussian functionn
def gauss(x,mu,sigma):
	y=1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2);
	return y


## Gaussian kernel for filtering
def gauss_kernel(sigma):
	# Return a kernel with Gaussian shape over 2-sigma range with 
	#  width = sigma (samples)
	x=np.arange(-2*sigma,2*sigma+1)
	y=1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*(x/sigma)**2);
	return y