"""
	** MCMC Incremental Slip Rate Calculator **
	Convert calendar years to years before present

	Rob Zinke 2019, 2020
"""

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from slipRateObjects import gauss_kernel


### PARSER ---
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Convert an OxCal date in CE/BCE to years before present.')
	parser.add_argument(dest='datefile',help='File containing ')
	parser.add_argument(dest='outName',type=str,help='Base name of output file.')
	parser.add_argument('-r','--reference_date',dest='ref_date',default=1950,type=float,help='Reference calendar date (C.E.). Default = 1950.')
	parser.add_argument('-f','--age_factor',dest='age_factor',default=1,type=float,help='Factor to divide age by (e.g., 1000). Default = 1.')
	parser.add_argument('-s','--smoothing',dest='smoothing',default=None,type=int,help='Width of smoothing kernel.')
	parser.add_argument('-p','--plot',dest='plot',action='store_true',help='Plot outputs.')
	return parser

def cmdParser(inpt_args=None):
	parser=createParser()
	return parser.parse_args(inpt_args)



### MAIN ---
if __name__ == '__main__':
	inpt=cmdParser()

	# Load input data
	PDFdate=np.loadtxt(inpt.datefile)
	xDate=PDFdate[:,0] # values
	pxDate=PDFdate[:,1] # probabilities

	# Convert date to age
	print('Relative to {}'.format(inpt.ref_date))
	xAge=inpt.ref_date-xDate
	xAge=np.flip(xAge)
	pxAge=np.flip(pxDate)

	if inpt.age_factor:
		xAge/=inpt.age_factor

	# Smooth if requested
	if inpt.smoothing:
		kernel=gauss_kernel(inpt.smoothing)
		pxAge=np.convolve(pxAge,kernel,'same')

	# Normalize area to 1.0
	P=np.trapz(pxAge,xAge)
	pxAge/=P 

	# Calculate CDF
	PxAge=cumtrapz(pxAge,xAge,initial=0)

	# Save to text file
	if inpt.outName:
		Fout=open(inpt.outName,'w')
		Fout.write('# Age ({0:1.1e}yr)\tProbability\n'.format(inpt.age_factor))
		for i in range(len(xAge)):
			Fout.write('{}\t{}\n'.format(xAge[i],pxAge[i]))
		Fout.close()
		print('Saved to: {}'.format(inpt.outName))

	# Plot if requested
	if inpt.plot:
		F=plt.figure()
		axDate=F.add_subplot(211)
		axDate.plot(xDate,pxDate,'k')
		axDate.invert_xaxis()
		axDate.set_ylabel('Cal date')
		axAge=F.add_subplot(212)
		axAge.plot(xAge,pxAge/pxAge.max(),'k',linewidth=2,zorder=2,label='PDF')
		axAge.plot(xAge,PxAge,'b',linewidth=2,zorder=1,label='CDF')
		axAge.legend()
		axAge.set_ylabel('Age ({0:1.1e} yr)'.format(inpt.age_factor))

		plt.show()