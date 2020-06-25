#!/usr/bin/env python3
'''
    ** MCMC Incremental Slip Rate Calculator **
    Plot displacements on a single figure

    Rob Zinke 2019, 2020
'''

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt


### PARSER ---
Description = '''Plot displacements on a single figure. This script is meant for display purposes only and
will not affect the slip rate calculations.

This script accepts a list of displacements, encoded with the sample type, sample label, and path to file,
separated by colons. I.e.,
DatumType:[SampleName]:[FilePath]

For example:
Displacement:Offset1:Offset1.txt

shows that the sample is a Displacement datum, labelled Offset1, and the path to file Offset1.txt is
provided. The displacement sample type can be color coded by type or user specification:
* Displacement
* Generic

Generic displacement PDFs can also be specified by the user with the --generic-color and --generic-alpha
options.

Other plot elements are also available, including line separators for ease of visualization:
* Line
* DashedLine
* BoldLine

The line elements do not require file specification. For example:
Line:Strat Boundary
'''

Examples = '''EXAMPLES

# Plot displacements from displacement list
plotDisplacements.py DspList.txt -y 'displacements (m)' -t 'Measurement Examples' -r 80 -o DspPlotExample --generic-color b
'''

def createParser():
    import argparse
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)
    parser.add_argument(dest='dspsList', type=str,
        help='Text document with list of displacement files. See above for example.')
    parser.add_argument('-r','--label-rotation', dest='labelRotation', default=0, type=float,
        help='Label rotation')
    parser.add_argument('-t','--title', dest='title', default=None, type=str,
        help='Base for graph title.')
    parser.add_argument('-y','--ylabel', dest='ylabel', default='Displacement', type=str,
        help='Y-axis label')
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



### PLOT DISPLACEMENTS ---
class dspPlot:
    '''
        Plot provided displacements on a single figure
    '''
    def __init__(self, dspsList):
        '''
            Establish figure and plot data. See examples for "dspsList".
        '''

        # Read and parse data
        self.__readData__(dspsList)

        # Initialize figure
        self.__setupFig__()


    def __setupFig__(self):
        '''
            Setup initial figure.
        '''
        self.Fig = plt.figure(figsize=(10,10))
        self.ax = self.Fig.add_subplot(111)

        self.labelList = [] # keep track of labels


    def __readData__(self,dspsList):
        '''
            Open list of displacement files, gather inputs as "data".
        '''
        # Read data from file
        with open(dspsList,'r') as dspListFile:
            self.data = dspListFile.readlines()
            dspListFile.close()

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
                self.ax.axvline(k, color = (0.7,0.75,0.8))

            # Plot dashed line
            elif datumType.lower() == 'dashedline':
                self.ax.axvline(k, color = (0.7,0.75,0.8), linestyle = '--')

            # Plot thicker line
            elif datumType.lower() == 'boldline':
                self.ax.axvline(k, color = (0.3,0.35,0.35))

            ## Assume displacement PDF
            else:
                # Load and format displacement datum
                xDsp, pxDsp = self.__formatDsp__(datum, pdfScale)

                # Shift
                pxDsp += k

                # Plot data using color from color table
                self.ax.fill(pxDsp, xDsp,
                    color = self.colors[datumType.lower()],
                    alpha = self.alphas[datumType.lower()])


            k += 1 # update counter


    def __formatDsp__(self, datum, pdfScale):
        '''
            Load and format displacement data.
        '''
        # Load data from file
        assert len(datum) == 3, 'Cannot plot displacement. Type, label, and filepath must be specified.'
        dspFile = datum[2]
        dspData = np.loadtxt(dspFile)
        xDsp = dspData[:,0]
        pxDsp = dspData[:,1]

        # Zero-pad
        xDsp = np.pad(xDsp,(1,1),'edge')
        pxDsp = np.pad(pxDsp,(1,1),'constant')
        # print(xDsp); exit()

        # Scale probability to 1.0 * scale factor
        pxDsp = pdfScale*pxDsp/pxDsp.max()

        return xDsp, pxDsp


    def __initColors__(self, genericColor, genericAlpha):
        '''
            Color lookup table.
        '''
        self.colors = {}
        self.alphas = {}
        # Non-specific displacement
        self.colors['displacement'] = (0,0,0)
        self.alphas['displacement'] = 1.0
        # Generic PDF
        self.colors['generic'] = genericColor
        self.alphas['generic'] = genericAlpha


    def finalizeFig(self, title=None, ylabel='Displacement', labelRotation=0, outName=None):
        '''
            Finalize figure.
        '''
        # X-labels
        ticks = np.arange(len(self.labelList))
        self.ax.set_xticks(ticks)
        self.ax.set_xticklabels(self.labelList, rotation=labelRotation)

        # Y-labels
        self.ax.set_ylabel(ylabel)

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

    # Plot displacements
    displacements = dspPlot(dspsList = inps.dspsList)

    displacements.plotData(pdfScale = inps.pdfScale,
        genericColor = inps.genericColor, genericAlpha = inps.genericAlpha)

    displacements.finalizeFig(title = inps.title,
        ylabel = inps.ylabel,
        labelRotation = inps.labelRotation,
        outName = inps.outName)


    plt.show()