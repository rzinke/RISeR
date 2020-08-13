#!/usr/bin/env python3
'''
    ** MCMC Incremental Slip Rate Calculator **
    Plot ages on a single figure

    Rob Zinke 2019, 2020
'''

### IMPORT MODULES ---
import argparse
try:
    import yaml
except:
    print('Please install pyyaml'); exit()
import numpy as np
import matplotlib.pyplot as plt


### PARSER ---
Description = '''Plot ages on a single figure. This script is meant for display purposes only and
will not affect the slip rate calculations.

This script accepts a list of ages, encoded with the sample type, sample label, and path to file,
separated by colons. I.e.,
Sample Name: {\"property\": \"value\", ...}

For example:
Sample 1: {"file": "Sample1_age.txt", "dtype": "age"}

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
* Dashed Line
* Bold Line

The line elements do not require file specification. For example:
Strat Boundary: {\"dtype\": \"bold line\"}
'''

Examples = '''EXAMPLES

# Plot ages from age list
plotAges.py AgeList.yaml -x 'ages (ka)' -t 'Sample Examples' -r 10 -o AgePlotExample
'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    parser.add_argument(dest='ageList', type=str,
        help='YAML document with list of age files. See above for example.')
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
        Plot provided ages on a single figure.
    '''
    def __init__(self, ageList):
        '''
            Establish figure and plot data. See examples for "ageList".
        '''

        # Initialize figure
        self.__setupFig__()

        # Read and parse data
        self.__loadData__(ageList)


    def __setupFig__(self):
        '''
            Setup initial figure.
        '''
        self.Fig = plt.figure(figsize=(10,10))
        self.ax = self.Fig.add_subplot(111)


    def __loadData__(self,ageList):
        '''
            Load age data files specifed in YAML format.
            Each entry gives the datum name and optional parameters.
            For use with the plotAges function.
        '''
        with open(ageList,'r') as ageFile:
            # Parse data within file
            self.ageData = yaml.load(ageFile, Loader=yaml.FullLoader)

            self.dataNames = list(self.ageData.keys())[::-1]


    def plotData(self, pdfScale=1.0, genericColor='k', genericAlpha=1.0):
        '''
            Work line by line to plot data based on datum type.
        '''

        # Initialize colors
        self.__initColors__(genericColor, genericAlpha)

        # Plot data one by one
        k = 0 # start counter

        for key in self.dataNames:
            ## Gather and format data
            datum = self.ageData[key]
            properties = list(datum.keys())
            properties = [property.lower() for property in properties]

            # Determine datum type
            if not 'dtype' in properties:
                # Default data type = age
                dtype = 'age'
            else:
                dtype = datum['dtype'].lower().replace(' ','')


            ## Plot line breaks
            # Plot simple line
            if dtype == 'line':
                self.ax.axhline(k, color = (0.7,0.75,0.8))

            # Plot dashed line
            elif dtype == 'dashedline':
                self.ax.axhline(k, color = (0.7,0.75,0.8), linestyle = '--')

            # Plot thicker line
            elif dtype == 'boldline':
                self.ax.axhline(k, color = (0.3,0.35,0.35))


            ## Assume age PDF
            else:
                # Plot prior if listed
                if 'prior' in properties:
                    self.__plotAge__(k,
                        fname = datum['prior'],
                        pdfScale = pdfScale,
                        color = self.colors['prior'],
                        alpha = self.alphas['prior'])

                # Plot age datum
                if 'color' not in properties:
                    color = self.colors[dtype]
                else:
                    color = datum['color']

                if 'alpha' not in properties:
                    alpha = self.alphas[dtype]
                else:
                    alpha = datum['alpha']

                self.__plotAge__(k = k,
                    fname = datum['file'],
                    pdfScale = pdfScale,
                    color = color,
                    alpha = alpha)


            ## Update counter
            k += 1


    def __plotAge__(self,k,fname,pdfScale,color,alpha):
        '''
            Plot age datum as filled PDF.
        '''
        # Load and format age datum
        xAge, pxAge = self.__formatAge__(fname, pdfScale)

        # Shift
        pxAge += k

        # Plot datum
        self.ax.fill(xAge, pxAge, color=color, alpha=alpha)


    def __formatAge__(self, fname, pdfScale):
        '''
            Load and format age data.
        '''
        # Load data from file
        ageData = np.loadtxt(fname)
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
        ticks = np.arange(len(self.ageData))
        self.ax.set_yticks(ticks)
        self.ax.set_yticklabels(self.dataNames, rotation=labelRotation)

        # X-labels
        self.ax.set_xlabel(xlabel)

        # Title
        if title: self.ax.set_title(title)

        # Finalize
        self.Fig.tight_layout()

        # Save to file
        if outName:
            savename = '{}.pdf'.format(outName)
            self.Fig.savefig(savename, type='pdf')
            print('Saved to: {}'.format(savename))



### MAIN ---
if __name__ == '__main__':
    ## Arguments and setup
    # Parser
    inps = cmdParser()

    # Plot ages
    ages = agePlot(inps.ageList)

    ages.plotData(pdfScale = inps.pdfScale,
        genericColor = inps.genericColor, genericAlpha = inps.genericAlpha)

    ages.finalizeFig(title = inps.title,
        xlabel = inps.xlabel,
        labelRotation = inps.labelRotation,
        outName = inps.outName)


    plt.show()