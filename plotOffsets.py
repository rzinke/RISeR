#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Parser
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Main function for calculating incremental slip rates.')
	parser.add_argument(dest='offsets_list',type=str,help='File with list of dates.')
	parser.add_argument('-u','--offset_units',dest='offset_units',default='m',type=str,help='Units of offset measurements')
	parser.add_argument('-s','--scale',dest='scale',default=1.0,type=float,help='Vertical scale for PDFs')
	parser.add_argument('-c','--color',dest='color',default='b',type=str,help='Color for generic PDFs')
	parser.add_argument('-m','--smoothing',dest='smoothing',default=None,type=int,help='If selected, specify kernel width')
	parser.add_argument('-r','--label_rotation',dest='label_rotation',default=0,type=float,help='Label rotation')
	parser.add_argument('-p','--plot',dest='plot',action='store_true',help='Plot outcome?')
	parser.add_argument('-t','--title',dest='title',default=None,type=str,help='Base for graph title.')
	parser.add_argument('-o','--outName',dest='outName',default=None,type=str,help='Base for graph title.')
	parser.add_argument('-a','--alpha',dest='alpha',default=1,type=float,help='Opacity value for \"contributing\" PDFs')
	parser.add_argument('--axis_values',dest='axis_values',default=None,type=str,help='Min/Max axis values: \'min max\'')
	return parser

def cmdParser(inpt_args=None):
	parser=createParser()
	return parser.parse_args(inpt_args)

# Plotting functions
def plotContributing(ax,x,px,alpha=1.0):
	ax.fill(x,px,color=(0.6,0.6,0.6),alpha=alpha)

def plotGeneric(ax,x,px,color):
	ax.fill(x,px,color=color)

def plotFullBreak(ax,y):
	ax.axhline(y,color=(0.3,0.35,0.35))

def plotBreak(ax,y):
	ax.axhline(y,color=(0.7,0.75,0.8),linestyle='--')


# Main
if __name__ == '__main__':
	# Parse command line arguments
	inpt=cmdParser()

	# Setup initial figure
	Fmain=plt.figure(figsize=(10,10))
	axMain=Fmain.add_subplot(111)

	label_list=[]
	k=0

	# Read lines of input file
	Fin=open(inpt.offsets_list,'r')
	Lines=Fin.readlines()
	Fin.close()

	# Flip so top line is at the chart top
	Lines=Lines[::-1]

	# Work line by line
	for Line in Lines:
		Line=Line.strip('\n') # remove trailing new line
		print(Line)
		Offsets=Line.split() # parse dates in line
		nOffsets=len(Offsets) # 
		
		# Work date by date
		for Offset in Offsets:
			offsetInfo=Offset.split(':')

			# Breaks
			if offsetInfo[0] in ['FullBreak','Break']:
				if offsetInfo[0] in ['FullBreak']:
					plotFullBreak(axMain,k+0.5)
					offsetname=offsetInfo[1]
				elif offsetInfo[0] in ['Break']:
					plotBreak(axMain,k+0.5)
					offsetname=offsetInfo[1]

			# Offsets
			else:
				offsettype=offsetInfo[0]
				filename=offsetInfo[1]
				offsetname=offsetInfo[2]

				# Load data
				offsetPDF=np.loadtxt(filename)
				x=offsetPDF[:,0]; px=offsetPDF[:,1]

				# Format probability curve
				if inpt.smoothing:
					smoothing_kernel=np.ones(inpt.smoothing)
					px=np.convolve(px,smoothing_kernel,'same')
				px=inpt.scale*px/px.max()+k

				# Plot
				if offsettype=='Contributing':
					plotContributing(axMain,x,px,alpha=inpt.alpha)
				else:
					plotGeneric(axMain,x,px,inpt.color)

		# Update counter for each line
		label_list.append(offsetname)
		k+=1

	# Finishing plot
	if inpt.title:
		axMain.set_title(inpt.title)
	axMain.set_xlabel('offset ({})'.format(inpt.offset_units))
	axMain.set_yticks(np.arange(0.5,k+1.5))
	axMain.set_yticklabels(label_list,rotation=inpt.label_rotation)
	if inpt.axis_values:
		vals=inpt.axis_values.split()
		axMain.set_xlim([float(vals[0]),float(vals[1])])

	if inpt.outName:
		savename='{}.png'.format(inpt.outName)
		Fmain.savefig(savename,dpi=600)
		print('Saved to: {}'.format(savename))


	if inpt.plot:
		plt.show()