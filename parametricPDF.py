#!/usr/bin/env python3
'''
Use this to create a probability density function (PDF) as a text file. 
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as intrp 

## Parser
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Create a probability density function (PDF) given reference points and a specified, parametric distribution.')
	parser.add_argument('-d', '--distribution', dest='dstrb', type=str, required=True, help='Distribution type [\'gauss\'/\'triangle\'/\'trapezoid\']]')
	parser.add_argument('-v', '--values', dest='values', type=str, required=True, help='Values to specify distribution.\nTwo values for gaussian ([1] +/- [2])\nThree values for triangle ([min], [pref], [max])\nFour values for trapezoid.')
	parser.add_argument('-o', '--output', dest='outName', type=str, required=True, help='Output file path/name')
	parser.add_argument('-n', '--nDataPts', dest='n', type=int, default=100, help='Number of data points. Default = 100')
	parser.add_argument('-p', '--plot', dest='plot', type=bool, default=False, help='Show plot of output [True/False]')
	parser.add_argument('-verb', '--verbose', dest='verbose', type=bool, default=False, help='Print outputs to command line [True/False]')
	return parser

def cmdParser(inpt_args=None):
	parser = createParser()
	return parser.parse_args(args=inpt_args)


## Create PDF
def makePDF(inpt):
	# Convert values to list
	values=[] # empty list
	[values.append(float(v)) for v in inpt.values.split(' ')]

	# Build distribution
	if inpt.dstrb.lower() in ['gaussian','gauss']:
		# Gaussian
		mu=values[0] # mean 
		sd=values[1] # std deviation
		x=np.linspace(mu-4*sd,mu+4*sd,inpt.n)
		px=gauss(x,mu,sd)
	elif inpt.dstrb.lower() in ['tri','triangle','triangular']:
		# Triangular
		x=np.linspace(values[0],values[2],inpt.n)
		Ix=intrp.interp1d(values,[0,1,0],kind='linear')
		px=Ix(x) # interpolate
	elif inpt.dstrb.lower() in ['trap','trapezoid','trapezoidal']:
		# Trapezoidal
		x=np.linspace(values[0],values[3],inpt.n)
		Ix=intrp.interp1d(values,[0,1,1,0],kind='linear')
		px=Ix(x) # interpolate

	# Normalize to area = 1.0
	P=np.trapz(px,x) # integrated area
	px/=P # normalize
	P=np.trapz(px,x) # normalized area (should be 1.0)

	# Save to text file
	Fout=open(inpt.outName,'w')
	Fout.write('Value, Probability\n')
	for i in range(inpt.n):
		Fout.write('{0:f}\t{1:f}\n'.format(x[i],px[i]))
	Fout.close()

	# Print outputs if specified
	if inpt.verbose is True:
		print('Created distribution: {0:s}'.format(inpt.dstrb))
		print('Final area: {0:f}'.format(P))
		print('Saved to: {0:s}'.format(inpt.outName))

	# Plot if specified
	if inpt.plot is True:
		F=plt.figure()
		ax=F.add_subplot(111)
		ax.plot(x,px,'b',linewidth=2)
		ax.set_title('Probability Density Function')
		ax.set_xlabel('values'); ax.set_ylabel('rel prob')

		plt.show()


# --- User inputs --- 
# Output file name
# outName=''

# Number of data points in output
# nDataPts=100 

# Values
# U=[4.2,0.5]
# if U is two numbers long, it will be a gaussian distribution ([1] +/- [2])
# if U is three numbers long, it will be triangle (min pref max) 
# if U is four numbers long, it will be a trapezoid (min low high max) 



# # --- PROCESSING ---
# # Processing 
# import numpy as np 
# from scipy import interpolate; interp=interpolate.interp1d 
# import matplotlib.pyplot as plt 
# from Gauss import * 

# U=np.array(U) # convert list to array 

# if len(U)==2:
# 	u=np.linspace(U[0]-4*U[1],U[0]+4*U[1],num=n)
# 	p=gauss(u,U[0],U[1]) # normal distribution 
# else:
# 	u=np.linspace(U[0],U[-1],num=n)
# 	if len(U)==3:
# 		pI=interp(U,np.array([0,1,0])) # triangular distribution 
# 		p=pI(u) # interpolate 
# 	elif len(U)==4:
# 		pI=interp(U,np.array([0,1,1,0])) # trapezoidal distribution 
# 		p=pI(u) # interpolate 

# # Normalize area to 1.0
# p=p/np.trapz(p,u)
# print('Area = {}'.format(np.trapz(p,u))) # confirm area = 1.0 

# # Write to text file 
# np.savetxt(outfile,(u,p)) 
# print('Saved to: {}'.format(outfile))

# # Plot figure 
# plt.figure(1) 
# plt.plot(u,p,'b')
# plt.xlabel('offset');plt.ylabel('rel. prob.')
# plt.title('PDF')

# plt.show()

# Definition of gaussian functionn
def gauss(x,mu,sigma):
	y=1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2);
	return y



## Main
if __name__ == '__main__':
	inpt=cmdParser()
	makePDF(inpt)