#!/usr/bin/env python3
"""
	** MCMC Incremental Slip Rate Calculator **
	View and compute statistics for a probability density function

	Rob Zinke 2019, 2020
"""

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt


### PARSER ---
# Command line parser
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Quickly view and compute statistics for a probability density function (PDF)')
	parser.add_argument(dest='pdfFile' ,help='File with PDF data, first column = x, second column = px.')
	parser.add_argument('-t','--title', dest='title', default=None, type=str, help='x-axis label')
	parser.add_argument('--x-label', dest='xLabel', default=None, type=str, help='x-axis label')
	parser.add_argument('--y-label', dest='yLabel', default='rel. probability', type=str, help='y-axis label')
	parser.add_argument('-o','--out-name', dest='outName',type=str,default=None, help='Save plot to outName.')
	return parser

def cmdParser(inpt_args=None):
	parser=createParser()
	return parser.parse_args(inpt_args)



### MAIN ---
if __name__ == '__main__':
	## Gather arguments
	inpt=cmdParser()


	## Load file
	PDF=np.loadtxt(inpt.pdfFile)
	x=PDF[:,0] # values
	px=PDF[:,1] # probabilities


	## Plot figure
	F=plt.figure()
	ax=F.add_subplot(111)
	ax.plot(x,px,'k',linewidth=2)
	if inpt.title:
		ax.set_title(inpt.title)
	if inpt.xLabel:
		ax.set_xlabel(inpt.xLabel)
	if inpt.yLabel:
		ax.set_ylabel(inpt.yLabel)
	if inpt.outName:
		if inpt.outName[-4:] in ['.png','.PNG']:
			inpt.outName=inpt.outName[:-4] 
		F.savefig('{}.png'.format(inpt.outName),dpi=600)

	plt.show()