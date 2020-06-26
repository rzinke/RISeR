#!/usr/bin/env python3
'''
    ** MCMC Incremental Slip Rate Calculator **
    Plot ages on a single figure

    Rob Zinke 2019, 2020
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt


### PARSER ---
Description = '''Plot ages on a single figure. This script is meant for display purposes only and
will not affect the slip rate calculations.

This script accepts a list of ages, encoded with the sample type, sample label, and path to file,
separated by colons. I.e.,
DatumType:[SampleName]:[FilePath]

For example:
Age:Sample1:Sample1_age.txt

shows that the sample is an Age datum, labelled Sample1, and the path to file Sample1_age.txt is
provided. The age sample type can be color coded by one of several statistical or stratigraphic
types:
* Prior
* Posterior
* PhaseAge
* BetweenAge
* Silt
* SandySilt
* Sand
* Gravel

Generic age PDFs can also be specified by the user with the --generic-color and --generic-alpha
options.

Other plot elements are also available, including line separators for ease of visualization:
* Line
* DashedLine
* BoldLine

The line elements do not require file specification. For example:
Line:Strat Boundary
'''

Examples = '''EXAMPLES

# Plot ages from age list
plotAges.py AgeList.txt -x 'ages (ka)' -t 'Sample Examples' -r 10 -o AgePlotExample
'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    parser.add_argument(dest='agesList', type=str,
        help='Text document with list of age files. See above for example.')
    parser.add_argument('-r','--label-rotation', dest='labelRotation', default=0, type=float,
        help='Label rotation')
    parser.add_argument('-t','--title', dest='title', default=None, type=str,
        help='Base for graph title.')
    parser.add_argument('-x','--xlabel', dest='xlabel', default='Age', type=str,
        help='X-axis label')
    parser.add_argument('-o','--outName', dest='outName', default=None, type=str,
        help='Base for graph title')
    parser.add_argument('--pdf-scale', dest='pdfScale', default=1.0, type=float,
        help='Vertical scale for PDFs')
    parser.add_argument('--generic-color', dest='genericColor', default='k',
        help='Color for generic plot')
    parser.add_argument('--generic-alpha', dest='genericAlpha', type=float, default=1.0,
        help='Opacity for generic plot')

    return parser

def cmdParser(inpt_args=None):
    parser=createParser()
    return parser.parse_args(inpt_args)



### PLOT AGES ---
class agePlot:
    '''
        Plot provided ages on a single figure
    '''
    def __init__(self, agesList):
        '''
            Establish figure and plot data. See examples for "agesList".
        '''

        # Read and parse data
        self.__readData__(agesList)

        # Initialize figure
        self.__setupFig__()


    def __setupFig__(self):
        '''
            Setup initial figure.
        '''
        self.Fig = plt.figure(figsize=(10,10))
        self.ax = self.Fig.add_subplot(111)

        self.labelList = [] # keep track of labels


    def __readData__(self,agesList):
        '''
            Open list of age files, gather inputs as "data".
        '''
        # Read data from file
        with open(agesList,'r') as ageListFile:
            self.data = ageListFile.readlines()
            ageListFile.close()

        # Flip so top line is at the chart top
        self.data = self.data[::-1]


    def plotData(self, pdfScale=1.0, genericColor='k', genericAlpha=1.0):
        '''
            Work line by line to plot data based on datum type.
        '''

        # Initialize colors
        self.__initColors__(genericColor, genericAlpha)

        # Plot data one by one
        k = 0 # start counter

        for datum in self.data:
            # Format and parse datum type
            datum = datum.strip('\n') # remove trailing new line
            datum = datum.split(':')

            # Update label list; leave blank if label not provided
            try:
                self.labelList.append(datum[1])
            except:
                self.labelList.append('')

            # Base action on datum type
            datumType = datum[0]

            ## Plot line breaks
            # Plot simple line
            if datumType.lower() == 'line':
                self.ax.axhline(k, color = (0.7,0.75,0.8))

            # Plot dashed line
            elif datumType.lower() == 'dashedline':
                self.ax.axhline(k, color = (0.7,0.75,0.8), linestyle = '--')

            # Plot thicker line
            elif datumType.lower() == 'boldline':
                self.ax.axhline(k, color = (0.3,0.35,0.35))

            ## Assume age PDF
            else:
                # Load and format age datum
                xAge, pxAge = self.__formatAge__(datum, pdfScale)

                # Shift
                pxAge += k

                # Plot data using color from color table
                self.ax.fill(xAge, pxAge,
                    color = self.colors[datumType.lower()],
                    alpha = self.alphas[datumType.lower()])


            k += 1 # update counter


    def __formatAge__(self, datum, pdfScale):
        '''
            Load and format age data.
        '''
        # Load data from file
        assert len(datum) == 3, 'Cannot plot age. Type, label, and filepath must be specified.'
        ageFile = datum[2]
        ageData = np.loadtxt(ageFile)
        xAge = ageData[:,0]
        pxAge = ageData[:,1]

        # Zero-pad
        xAge = np.pad(xAge,(1,1),'edge')
        pxAge = np.pad(pxAge,(1,1),'constant')

        # Scale probability to 1.0 * scale factor
        pxAge = pdfScale*pxAge/pxAge.max()

        return xAge, pxAge


    def __initColors__(self, genericColor, genericAlpha):
        '''
            Color lookup table.
        '''
        self.colors = {}
        self.alphas = {}
        # Non-specific age
        self.colors['age'] = (0,0,0)
        self.alphas['age'] = 1.0
        # Prior age
        self.colors['prior'] = (0.6,0.6,0.6)
        self.alphas['prior'] = 1.0
        # Posterior age
        self.colors['posterior'] = (0,0,0)
        self.alphas['posterior'] = 1.0
        # Phase age
        self.colors['phaseage'] = 'm'
        self.alphas['phaseage'] = 1.0
        # Between age
        self.colors['betweenage'] = 'b'
        self.alphas['betweenage'] = 1.0
        # Silt age
        self.colors['silt'] = (0.843,0.867,0.235)
        self.alphas['silt'] = 1.0
        # Sandy silt age
        self.colors['sandysilt'] = (0.529,0.831,0.898)
        self.alphas['sandysilt'] = 1.0
        # Sand age
        self.colors['sand'] = (0.961,0.580,0.192)
        self.alphas['sand'] = 1.0
        # Gravel age
        self.colors['gravel'] = (0.922,0.129,0.184)
        self.alphas['gravel'] = 1.0
        # Generic PDF
        self.colors['generic'] = genericColor
        self.alphas['generic'] = genericAlpha


    def finalizeFig(self, title=None, xlabel='Age', labelRotation=0, outName=None):
        '''
            Finalize figure.
        '''
        # Y-labels
        ticks = np.arange(len(self.labelList))
        self.ax.set_yticks(ticks)
        self.ax.set_yticklabels(self.labelList, rotation=labelRotation)

        # X-labels
        self.ax.set_xlabel(xlabel)

        # Title
        if title: self.ax.set_title(title)

        # Finalize
        self.Fig.tight_layout()

        # Save to file
        if outName:
            savename = '{}.png'.format(outName)
            self.Fig.savefig(savename, dpi=600)
            print('Saved to: {}'.format(savename))



### MAIN ---
if __name__ == '__main__':
    ## Arguments and setup
    # Parser
    inps = cmdParser()

    # Plot ages
    ages = agePlot(agesList = inps.agesList)

    ages.plotData(pdfScale = inps.pdfScale,
        genericColor = inps.genericColor, genericAlpha = inps.genericAlpha)

    ages.finalizeFig(title = inps.title,
        xlabel = inps.xlabel,
        labelRotation = inps.labelRotation,
        outName = inps.outName)


    plt.show()