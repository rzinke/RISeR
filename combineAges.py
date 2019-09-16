#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from PDFanalysis import HPDpdf


# Command line parser
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Combine (sum) age probability density functions. Provide a list of ages.')
	parser.add_argument(dest='ageList',nargs='+',help='List of ages to combine, separated by spaces.')
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

	# Loop through ages to get basic parameters
	nAges=len(inpt.ageList)
	Ages={} # store age data here
	minAge=[]; maxAge=[]
	meanSampling=[]

	for age in inpt.ageList:
		print('Age: {}'.format(age))
		# Load data
		PDF=np.loadtxt(age)
		Ages[age]=PDF # store PDF to dictionary for later
		x=PDF[:,0]; px=PDF[:,1]

		# Assure probability = 1.0
		P=np.trapz(px,x)
		PDF[:,1]=PDF[:,1]/P # area = 1

		print(np.trapz(PDF[:,1],PDF[:,0]))

		# Stats
		minAge.append(x.min())
		maxAge.append(x.max())
		meanSampling.append(np.mean(np.diff(x)))

		if inpt.verbose:
			print('min: {} / max {} / sampling {}'.format(\
				minAge[-1],maxAge[-1],meanSampling[-1]))

	# Interpolate along common axis
	minAge=np.min(minAge) 
	maxAge=np.max(maxAge)
	sampling=np.min(meanSampling)
	xCombo=np.arange(minAge,maxAge+sampling,sampling)

	if inpt.verbose:
		print('Formed common axis:\nmin {} / max {} / sampling {}'.format(\
			minAge,maxAge,sampling))

	# Interpolate age probabilities
	if inpt.method in ['union','add','sum','+']:
		pxCombo=np.zeros(len(xCombo))
	elif inpt.method in ['intersection','multiply','product','x']:
		pxCombo=np.ones(len(xCombo))

	for age in inpt.ageList:
		# Recall ages and probs from Ages dictionary
		x=Ages[age][:,0]
		px=Ages[age][:,1]

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
		Fout.write('# Age\tProbability\n')
		for i in range(len(xCombo)):
			Fout.write('{}\t{}\n'.format(xCombo[i],pxCombo[i]))
		Fout.close()
		print('Saved to: {}'.format(inpt.outName))

	# Plot
	if inpt.plot:
		F=plt.figure()
		ax=F.add_subplot(111)
		# Add contributing PDFs
		for age in inpt.ageList:
			ageLabel=age.split('/')[-1]
			x=Ages[age][:,0]; px=Ages[age][:,1]
			ax.plot(x,px,color=(0.6,0.6,0.6))
			ax.text(x[px==px.max()][0],px[px==px.max()][0],age)
		ax.plot(xCombo,pxCombo,color='k',linewidth=2)
		ax.set_xlabel('age'); ax.set_ylabel('rel prob')
		ax.set_title('Combined age PDF')

		plt.show()