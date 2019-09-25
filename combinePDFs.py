#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from PDFanalysis import HPDpdf


# Command line parser
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Combine probability density functions (PDFs) by pointwise summing or mulitplication.')
	parser.add_argument(dest='pdfList',nargs='+',help='List of PDFs to combine: Filenames separated by spaces.')
	parser.add_argument('-o','--outName',dest='outName',type=str,help='Output name.')
	parser.add_argument('-m','--method',dest='method',default='union',type=str,help='Method for combining (union or intersection). Default = union.')
	parser.add_argument('-v','--verbose',dest='verbose',action='store_true',help='Verbose mode?')
	parser.add_argument('-p','--plot',dest='plot',action='store_true',help='Plot?')
	return parser


def cmdParser(inpt_args=None):
	parser=createParser()
	return parser.parse_args(inpt_args)


if __name__ == '__main__':
	inpt=cmdParser()

	# Loop through PDFs to get basic parameters
	nPDFs=len(inpt.pdfList)
	PDFs={} # store pdf data here
	minx=[]; maxx=[]
	meanSampling=[]

	for pdf in inpt.pdfList:
		print('PDF: {}'.format(pdf))
		# Load data
		PDF=np.loadtxt(pdf)
		PDFs[pdf]=PDF # store PDF to dictionary for later
		x=PDF[:,0]; px=PDF[:,1]

		# Assure probability = 1.0
		P=np.trapz(px,x)
		PDF[:,1]=PDF[:,1]/P # area = 1

		print(np.trapz(PDF[:,1],PDF[:,0]))

		# Stats
		minx.append(x.min())
		maxx.append(x.max())
		meanSampling.append(np.mean(np.diff(x)))

		if inpt.verbose:
			print('min: {} / max {} / sampling {}'.format(\
				minx[-1],maxx[-1],meanSampling[-1]))

	# Interpolate along common axis
	minx=np.min(minx) 
	maxx=np.max(maxx)
	sampling=np.min(meanSampling)
	xCombo=np.arange(minx,maxx+sampling,sampling)

	if inpt.verbose:
		print('Formed common axis:\nmin {} / max {} / sampling {}'.format(\
			minx,maxx,sampling))

	# Interpolate probabilities
	if inpt.method in ['union','add','sum','+']:
		pxCombo=np.zeros(len(xCombo))
	elif inpt.method in ['intersection','multiply','product','x']:
		pxCombo=np.ones(len(xCombo))

	for pdf in inpt.pdfList:
		# Recall values and probs from PDFs dictionary
		x=PDFs[pdf][:,0]
		px=PDFs[pdf][:,1]

		# Interpolate
		I=interp1d(x,px,bounds_error=False,fill_value=0.,kind='linear')

		# Combine using method of choice
		if inpt.method in ['union','add','sum','+']:
			# Union sums PDFs
			pxCombo+=I(xCombo) # add to cumulative function
		elif inpt.method in ['intersection','multiply','product','x']:
			# Intersection of PDFs
			pxCombo*=I(xCombo) # multiply cumulative function

	# Normalize combined area to 1.0
	P=np.trapz(pxCombo,xCombo)
	pxCombo/=P # area = 1.0

	# Report final statistics if requested
	if inpt.verbose:
		HPDpdf(xCombo,pxCombo,confidence=95.45,verbose=True)

	# Save to file
	if inpt.outName:
		Fout=open(inpt.outName,'w')
		Fout.write('# value\tProbability\n')
		for i in range(len(xCombo)):
			Fout.write('{}\t{}\n'.format(xCombo[i],pxCombo[i]))
		Fout.close()
		print('Saved to: {}'.format(inpt.outName))

	# Plot
	if inpt.plot:
		F=plt.figure()
		ax=F.add_subplot(111)
		# Add contributing PDFs
		for pdf in inpt.pdfList:
			pdfLabel=pdf.split('/')[-1]
			x=PDFs[pdf][:,0]; px=PDFs[pdf][:,1]
			ax.plot(x,px,color=(0.6,0.6,0.6))
			ax.text(x[px==px.max()][0],px[px==px.max()][0],pdf)
		ax.plot(xCombo,pxCombo,color='k',linewidth=2)
		ax.set_xlabel('value'); ax.set_ylabel('rel prob')
		ax.set_title('Combined PDF')

		plt.show()