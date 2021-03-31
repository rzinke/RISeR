#!/usr/bin/env python3
'''
** RISeR Incremental Slip Rate Calculator **
View and compute statistics for a probability density function.

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from resultSaving import confirmOutputDir


### PARSER ---
# Command line parser
def createParser():
    import argparse
    parser = argparse.ArgumentParser(description='Quickly view and compute statistics for a probability density function (PDF)')
    parser.add_argument(dest='pdfFile', type=str,
        help='File with PDF data, first column = x, second column = px.')
    parser.add_argument('-t','--title', dest='title', default=None, type=str, 
        help='x-axis label')
    parser.add_argument('-x','--x-label', dest='xLabel', default=None, type=str, 
        help='x-axis label')
    parser.add_argument('--fill-color', dest='fillColor', type=str, default=None,
        help='PDF fill color.')
    parser.add_argument('--show-stats', dest='showStats', action='store_true', 
        help='Plot lines representing statistics.')
    parser.add_argument('-o','--out-name', dest='outName', type=str, default=None, 
        help='Save plot to outName.')
    return parser

def cmdParser(inps=None):
    parser = createParser()
    return parser.parse_args(inps)



### STATISTICS ---
def checkArea(x, px):
    '''
    Check that area has unity mass.
    '''
    # Calculate area
    P = np.trapz(px, x)

    # Report if area is not 1.0
    if np.abs(1-P) > 1E-6: print('WARNING: Mass is not unity')

class pdfStats:
    def __init__(self, x, px):
        '''
        Calculate the statistics of a PDF.
        '''
        # Extrema
        self.__findExtrema__(x)

        # Compute mass
        self.P = np.trapz(px, x)

        # Percentiles
        self.__findPercentiles__(x, px)

        # Find mode
        self.__findMode__(x, px)

        # Compute mean
        self.__findMean__(x, px)

    def __findExtrema__(self, x):
        '''
        Find min/max extrema.
        '''
        self.min = x.min()
        self.max = x.max()

    def __findPercentiles__(self, x, px):
        '''
        Find percentiles.
        '''
        # Compute CDF
        Px = cumtrapz(px, x, initial=0)

        # Inverse transform
        Icdf=interp1d(Px, x, kind='linear', bounds_error=False, fill_value=np.nan)

        # 50 percentile (median)
        self.median = Icdf(50/100)

        # 68.27 percentile
        self.lower68, self.upper68 = self.__calcPercentile__(x, px, Icdf, 68.27)

        # 95.45 percentile
        self.lower95, self.upper95 = self.__calcPercentile__(x, px, Icdf, 95.45)

        # 99.73 percentile
        self.lower99, self.upper99 = self.__calcPercentile__(x, px, Icdf, 99.73)

    def __calcPercentile__(self, x, px, Icdf, pct):
        '''
        Calculate the nth percentile using the inverse transform.
        '''
        lowerPct = Icdf((50-pct/2)/100)
        upperPct = Icdf((50+pct/2)/100)

        return lowerPct, upperPct

    def __findMode__(self, x, px):
        '''
        Find mode. Use median of index values if more than one maximum detected.
        '''
        modeValues = x[px==px.max()]
        self.mode = np.median(modeValues)

    def __findMean__(self, x, px):
        '''
        Find mean.
        '''
        dx = np.mean(np.diff(x))
        self.mean = np.sum(dx*x*px)


def reportStats(stats):
    '''
    Print statistics.
    '''
    # Print basic statistics
    print('''Mass: {P:.7f}
Min:       {min:.4f}
Lower 95%: {lower95:.4f}
Lower 68%: {lower68:.4f}
Mean:      {mean:.4f}
Median:    {median:.4f}
Mode:      {mode:.4f}
Upper 68%: {upper68:.4f}
Upper 95%: {upper95:.4f}
Max:       {max:.4f}'''.format(**vars(stats)))

    # Print mode +- range
    plus = stats.upper95 - stats.mode
    minus = stats.mode - stats.lower95
    print('\nMode +- 95% range:')
    print('{:.2f} +{:.2f} -{:.2f}'.format(stats.mode, plus, minus))



### PLOTTING ---
def plotPDF(x, px, stats, title=None, xLabel=None, fillColor=None, showStats=False):
    '''
    Plot a probability density function.
    Optionally plot statistics.
    '''
    # Spawn figure
    fig, ax = plt.subplots()

    # Plot PDF
    if fillColor is not None:
        # Plot filled
        ax.fill(x, px, color=fillColor, zorder=1, label='data')
    else:
        # Plot line
        ax.plot(x, px, color='k', linewidth=2, zorder=2, label='data')

    # Plot statistics if requested
    if showStats == True: plotStats(ax, x, px, stats)

    # Format plot
    ax.legend()
    ax.set_xlabel(xLabel)
    ax.set_ylabel('rel. prob.')
    ax.set_title(title)

    return fig, ax

def plotStats(ax, x, px, stats):
    '''
    Plot relevant statistical values on PDF plot.
    '''
    # Plot percentiles
    color99 = (0.9, 0.9, 0.9)
    plotField(ax, x, px, stats.lower99, stats.upper99, color=color99)

    color95 = (0.75, 0.75, 0.75)
    plotField(ax, x, px, stats.lower95, stats.upper95, color=color95)

    color68 = (0.5, 0.5, 0.5)
    plotField(ax, x, px, stats.lower68, stats.upper68, color=color68)

    # Plot central values
    plotIntersection(ax, x, px, stats.mean, color=(0.58, 0.00, 0.83), label='mean')
    plotIntersection(ax, x, px, stats.median, color='b', label='median')
    plotIntersection(ax, x, px, stats.mode, color='r', label='mode')

def plotField(ax, x, px, lower, upper, color):
    '''
    Plot a filled-in field.
    '''
    # Find values that meet conditions
    ndx = (x>=lower) & (x<=upper)

    # Pad values at edges
    xpad = np.pad(x[ndx], (1, 1), 'edge')
    pxpad = np.pad(px[ndx], (1, 1), 'constant')

    # Plot field
    ax.fill(xpad, pxpad, color=color, zorder=1)

def plotIntersection(ax, x, px, xValue, color='k', label=None):
    '''
    Plot the first instance of intersection with the curve.
    '''
    # Build interpolation
    intp = interp1d(x, px, kind='linear')

    # Determine value probability
    pxValue = intp(xValue)

    # Plot values
    ax.plot([xValue, xValue], [0, pxValue], color=color, label=label)
    ax.plot(xValue, pxValue, color=color, marker='.')



### SAVING ---
def saveFigure(fig, outName):
    '''
    Save figure to file.
    '''
    # Confirm output directory exists
    confirmOutputDir(inps.outName)

    # Format save name
    if outName[-4:] not in ['.png', '.PNG']: outName += '.png'

    # Save
    fig.savefig(outName, dpi=600)



### MAIN ---
if __name__ == '__main__':
    # Gather arguments
    inps = cmdParser()

    # Load and parse file
    PDF = np.loadtxt(inps.pdfFile)
    x = PDF[:,0]  # values
    px = PDF[:,1]  # probabilities

    # Check area = 1.0
    checkArea(x, px)

    # Compute statistics
    stats = pdfStats(x, px)

    # Report stats
    reportStats(stats)

    # Plot figure
    fig, ax = plotPDF(x, px, stats, 
        title=inps.title, xLabel=inps.xLabel,
        fillColor=inps.fillColor,
        showStats=inps.showStats)

    # Save if requested
    if inps.outName:
        saveFigure(fig, inps.outName)


    plt.show()