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
    parser.add_argument('--x-label', dest='xLabel', default=None, type=str, 
        help='x-axis label')
    parser.add_argument('--y-label', dest='yLabel', default='rel. probability', type=str, 
        help='y-axis label')
    parser.add_argument('--show-stats', dest='showStats', action='store_true', 
        help='Plot lines representing statistics.')
    parser.add_argument('-o','--out-name', dest='outName', type=str, default=None, 
        help='Save plot to outName.')
    return parser

def cmdParser(inps=None):
    parser = createParser()
    return parser.parse_args(inps)



### ANCILLARY FUNCTIONS ---
class pdfStats:
    # Initialize
    def __init__(self, x, px):
        ## Simple stats
        self.min = x.min()
        self.max = x.max()

        # Mode
        self.mode = x[px==px.max()][0]


        ## Percentiles
        # Compute CDF
        Px = cumtrapz(px, x, initial=0)

        # Interpolate percentiles
        Icdf=interp1d(Px, x, kind='linear', bounds_error=False, fill_value=np.nan)

        # Find percentiles
        confidenceInterval = 95.45
        self.lower = Icdf((50-confidenceInterval/2)/100)
        self.upper = Icdf((50+confidenceInterval/2)/100)
        self.median = Icdf(50/100)


        ## Mean
        dx = np.mean(np.diff(x))
        self.mean = np.sum(dx*x*px)



### MAIN ---
if __name__ == '__main__':
    # Gather arguments
    inps = cmdParser()

    # Load and parse file
    PDF = np.loadtxt(inps.pdfFile)
    x = PDF[:,0]  # values
    px = PDF[:,1]  # probabilities

    # Check area = 1.0
    P = np.trapz(px, x)
    if np.abs(1-P) > 1E-14:
        print('WARNING: Mass is not unity')

    # Compute statistics
    stats = pdfStats(x, px)

    # Report stats
    print('Mass: {}'.format(P))
    print('Min: {:.4f}'.format(stats.min))
    print('Lower 95.45%: {:.4f}'.format(stats.lower))
    print('Mean: {:.4f}'.format(stats.mean))
    print('Median: {:.4f}'.format(stats.median))
    print('Upper 95.45%: {:.4f}'.format(stats.upper))
    print('Max: {:.4f}'.format(stats.max))
    print('Mode 95%: {:.2f} +{:.2f} -{:.2f}'.format(stats.mode, stats.upper-stats.mode, stats.mode-stats.lower))


    # Plot figure
    fig, ax = plt.subplots()
    ax.plot(x, px, 'k', linewidth=2, zorder=2, label='data')
    if inps.title: ax.set_title(inps.title)
    if inps.xLabel: ax.set_xlabel(inps.xLabel)
    if inps.yLabel: ax.set_ylabel(inps.yLabel)

    # Label statistics if requested
    if inps.showStats is True:
        ax.axvline(stats.lower, color='b', alpha=0.25, zorder=1, label='95.45%')
        ax.axvline(stats.upper, color='b', alpha=0.25, zorder=1, label='95.45%')
        ax.axvline(stats.median, color='b', alpha=0.5, zorder=1, label='median')
        ax.axvline(stats.mode, color='r', alpha=0.5, zorder=1, label='mode')
        ax.axvline(stats.mean, color=(148/255,0/255,211/255), alpha=0.5, zorder=1, label='mean')
        ax.legend()

    # Save if requested
    if inps.outName:
        # Confirm output directory exists
        confirmOutputDir(inps.outName)

        # Format save name
        if inps.outName[-4:] in ['.png', '.PNG']: inps.outName = inps.outName[:-4]

        # Save
        fig.savefig('{:s}.png'.format(inps.outName), dpi=600)


    plt.show()