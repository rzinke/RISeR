#!/usr/bin/env python3
"""
	** MCMC Incremental Slip Rate Calculator **
	Plot ages on a single figure

	Rob Zinke 2019, 2020
"""

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt


### PARSER ---
def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Plot ages on a single figure.')
	parser.add_argument(dest='datesList', type=str, help='File with list of dates.')
	parser.add_argument('-d','--date-offset', dest='dateOffset', default=None, type=int, help='Reference age for date (e.g., 2019)')
	parser.add_argument('-s','--scale', dest='scale', default=1.0, type=float, help='Vertical scale for PDFs')
	parser.add_argument('-m','--smoothing', dest='smoothing', default=None, type=int, help='If selected, specify kernel width')
	parser.add_argument('-r','--label-rotation', dest='labelRotation', default=0, type=float, help='Label rotation')
	parser.add_argument('-p','--plot', dest='plot', action='store_true', help='Plot outcome?')
	parser.add_argument('-t','--title', dest='title', default=None, type=str, help='Base for graph title.')
	parser.add_argument('-o','--outName', dest='outName', default=None, type=str, help='Base for graph title.')

	parser.add_argument('--generic-color',dest='genericColor',default='k',help='Color for generic plot.')
	parser.add_argument('--generic-alpha',dest='genericAlpha',type=float,default=1.0,help='Opacity for generic plot.')
	return parser

def cmdParser(inpt_args=None):
	parser=createParser()
	return parser.parse_args(inpt_args)



### PLOTTING FUNCTIONS ---
def plotPrior(ax,x,px):
	ax.fill(x,px,color=(0.6,0.6,0.6))

def plotPosterior(ax,x,px):
	ax.fill(x,px,color=(0,0,0))

def plotSilt(ax,x,px):
	ax.fill(x,px,color=(0.843,0.867,0.235))

def plotSandySilt(ax,x,px):
	ax.fill(x,px,color=(0.529,0.831,0.898))

def plotSand(ax,x,px):
	ax.fill(x,px,color=(0.961,0.580,0.192))

def plotGravel(ax,x,px):
	ax.fill(x,px,color=(0.922,0.129,0.184))

def plotBetweenAge(ax,x,px):
	ax.fill(x,px,color='b')

def plotPhaseAge(ax,x,px):
	ax.fill(x,px,color='m')

def plotGeneric(ax,x,px,color,alpha):
	ax.fill(x,px,color=color,alpha=alpha)

def plotFullBreak(ax,y):
	ax.axhline(y,color=(0.3,0.35,0.35))

def plotBreak(ax,y):
	ax.axhline(y,color=(0.7,0.75,0.8),linestyle='--')



### MAIN ---
if __name__ == '__main__':
	## Arguments and setup
	# Parser
	inpt=cmdParser()

	# Setup initial figure
	Fmain=plt.figure(figsize=(10,10))
	axMain=Fmain.add_subplot(111)

	labelList=[]
	k=0


	## Plot data
	# Read lines of input file
	Fin=open(inpt.datesList,'r')
	Lines=Fin.readlines()
	Fin.close()

	# Flip so top line is at the chart top
	Lines=Lines[::-1]

	# Work line by line
	for Line in Lines:
		Line=Line.strip('\n') # remove trailing new line
		print(Line)
		Dates=Line.split() # parse dates in line
		nDates=len(Dates) # 
		
		# Work date by date
		for Date in Dates:
			dateInfo=Date.split(':')

			# Breaks
			if dateInfo[0] in ['FullBreak','Break']:
				if dateInfo[0] in ['FullBreak']:
					plotFullBreak(axMain,k+0.5)
					datename=dateInfo[1]
				elif dateInfo[0] in ['Break']:
					plotBreak(axMain,k+0.5)
					datename=dateInfo[1]

			# Dates
			else:
				datetype=dateInfo[0]
				filename=dateInfo[1]
				datename=dateInfo[2]

				# Load data
				datePDF=np.loadtxt(filename)
				x=datePDF[:,0]; px=datePDF[:,1]

				# Age/Calendar years
				if inpt.dateOffset:
					x=inpt.dateOffset-x

				# Format probability curve
				if inpt.smoothing:
					smoothingKernel=np.ones(inpt.smoothing)
					px=np.convolve(px,smoothingKernel,'same')
				px=inpt.scale*px/px.max()

				x=np.pad(x,(1,1),'edge')
				px=np.pad(px,(1,1),'constant')

				px+=k

				# Plot
				if datetype=='Prior':
					plotPrior(axMain,x,px)
				elif datetype=='Posterior':
					plotPosterior(axMain,x,px)
				elif datetype=='Silt':
					plotSilt(axMain,x,px)
				elif datetype=='SandySilt':
					plotSandySilt(axMain,x,px)
				elif datetype=='Sand':
					plotSand(axMain,x,px)
				elif datetype=='Gravel':
					plotGravel(axMain,x,px)
				elif datetype.lower() in ['date','between_age','betweenage','riser']:
					plotBetweenAge(axMain,x,px)
				else:
					plotGeneric(axMain,x,px,inpt.genericColor,inpt.genericAlpha)

		# Update counter for each line
		labelList.append(datename)
		k+=1


	## Finishing plot
	if inpt.title:
		axMain.set_title(inpt.title)
	if inpt.dateOffset:
		axMain.set_xlabel('age yb{}'.format(inpt.dateOffset))
	else:
		axMain.set_xlabel('age')
	axMain.invert_xaxis()
	axMain.set_yticks(np.arange(0.5,k+1.5))
	axMain.set_yticklabels(labelList,rotation=inpt.labelRotation)
	if inpt.outName:
		savename='{}.png'.format(inpt.outName)
		Fmain.savefig(savename,dpi=600)
		print('Saved to: {}'.format(savename))


	if inpt.plot:
		plt.show()