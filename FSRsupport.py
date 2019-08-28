'''
Objects and support functions for MCMC slip rate calculator
'''

import numpy as np 

def gauss(x,mu,sigma):
	y=1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2);
	return y

def erf(x,mu,sigma):
	# Do this the stupid way, haha 
	y=gauss(x,mu,sigma)
	Y=integrate.cumtrapz(y,x,initial=0) 
	Y=Y/Y.max() # normalize max value to 1.0
	return Y 