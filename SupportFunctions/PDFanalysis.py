'''
** RISeR Incremental Slip Rate Calculator **
Find inter-quantile range (IQR) or highest posterior density (HPD)
 of a probability density function (PDF).

Rob Zinke 2019-2021
'''

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
try:
    from scipy.integrate import cumulative_trapezoid as sp_cumulative_trapezoid
except ImportError:
    from scipy.integrate import cumtrapz as sp_cumulative_trapezoid


### PDF ANALYSIS CLASSES ---
## IQR method
class IQRpdf:
    '''
    Find values of PDF using inter-quantile range method (more stable, but biased toward skewed tail)
    INPUTS
        x is an array of evenly spaced values -- even spacing is important!
        px is an array of probabilities at those values
        confidence is the confidence interval, in percent (e.g., 95, 68)
    '''
    def __init__(self, x, px, confidence, outName=None, verbose=False):
        if verbose is True:
            print('Calculating interquantile range at {:.2f}% confidence limits'.format(confidence))

        # Determine bounds
        self.confidence = confidence  # record percent to object
        confidence /= 100  # percent to fraction
        lower = 0.5-confidence/2; upper = 0.5+confidence/2

        # Integrate to CDF
        P = np.trapz(px, x) # check that area = 1.0
        px /= P  # normalize
        Px = cumulative_trapezoid(px, x, initial=0)

        # Interpolate CDF to back-calculate values
        Icdf = interp1d(Px, x, kind='linear')
        lowerValue, median, upperValue = Icdf([lower, 0.5, upper])

        # Values within confidence range
        plotNdx = (x>=lowerValue) & (x<=upperValue)
        xIQR = x[plotNdx]; xIQR=np.pad(xIQR, (1,1), 'edge')
        pxIQR = px[plotNdx]; pxIQR = np.pad(pxIQR, (1,1), 'constant')

        # Outputs
        if verbose is True:
            print('\tLower value: {0:5f};\tUpper value: {1:5f}'.format(lowerValue, upperValue))

        self.x = x; self.px = px  # PDF
        self.Px = Px  # CDF
        self.Icdf = Icdf  # interpolation function
        self.median = median
        self.lowerValue = lowerValue
        self.upperValue = upperValue
        self.xIQR = xIQR  # x values within confidence range
        self.pxIQR = pxIQR  # px values within confidence range


    ## Plot input data function
    def plotInput(self):
        fig = plt.figure()
        axPDF = fig.add_subplot(121)
        axPDF.plot(self.x, self.px, color='k', linewidth=2)
        axPDF.plot([self.lowerValue, self.upperValue], [0,0], 'bx')
        axPDF.set_ylabel('rel prob')
        axCDF = fig.add_subplot(122)
        axCDF.plot(self.x, self.Px, color='k', linewidth=2)
        axCDF.plot([self.lowerValue, self.lowerValue], [self.lower, 0], color='b',zorder=1)
        axCDF.plot([self.upperValue, self.upperValue], [self.upper, 0], color='b',zorder=1)
        axCDF.set_xlabel('value'); axCDF.set_ylabel('integ\'d prob')


    ## Plot output function
    def plotOutput(self,outName=None):
        # Format ranges for plotting
        x = np.pad(self.x, (1,1), 'edge')
        px = np.pad(self.px, (1,1), 'constant')
        xIQR = np.pad(self.xIQR, (1,1), 'edge')
        pxIQR = np.pad(self.pxIQR, (1,1), 'constant')
        # Plot
        fig, ax = plt.subplots()
        ax.fill(x, px, color=(0.4,0.4,0.4), zorder=2)
        ax.fill(xIQR, pxIQR, color=(0.3,0.3,0.6), zorder=3)
        ax.set_xlabel('value'); ax.set_ylabel('ref prob')
        if outName:
            Fpdff.savefig(outName, dpi=600)


## HPD method
class HPDpdf:
    '''
    Find values of PDF using highest posterior density method (more representative of probable values)
    INPUTS
        x is an array of evenly spaced values -- even spacing is important!
        px is an array of probabilities at those values
        confidence is the confidence interval, in percent (e.g., 95, 68)
    '''
    def __init__(self, x, px, confidence, step_tolerance=1.05, outName=None, verbose=False):
        if verbose is True:
            print('Calculating highest posterior density at {}% confidence'.format(confidence))

        # Convert confidence bound from percent to fraction
        self.confidence = confidence  # record percent to object
        confidence /= 100

        # Confirm area = 1.0
        P = np.trapz(px, x)  # total area
        px /= P  # normalize

        # Check conditioning - points must be spaced approximately evenly
        xsteps = np.diff(x)  # spacing between sample points
        avestep = np.mean(xsteps)  # average spacing between sample points
        if np.sum(xsteps>step_tolerance*avestep) > 0:
            print('WARNING: Sample spacing must be approximately equal for HPD calculation')

        # Sort from highest-lowest probability
        sortNdx = np.argsort(px)  # indices lowest-highest
        sortNdx = sortNdx[::-1]  # flip highest-lowest

        # Sum probabilities until they reach the specified confidence limit
        xSort = x[sortNdx]; pxSort = px[sortNdx]  # place arrays in order
        PxSort = np.cumsum(pxSort)  # cumulative probability
        PxSort /= PxSort.max()  # normalize sum to 1.0

        # Find values that sum to confidence limit
        sortNdxRelevant = sortNdx[PxSort<=confidence]  # relevant index values
        xRelevant = x[sortNdxRelevant]  # x values
        pxRelevant = px[sortNdxRelevant]  # px values

        # Sort by x value
        sortBack = np.argsort(xRelevant)
        xRelevant = xRelevant[sortBack]
        pxRelevant = pxRelevant[sortBack]

        # Min and Max values
        lowestValue = xRelevant.min()
        highestValue = xRelevant.max()

        # Identify "clusters"
        nRelevant = len(xRelevant)  # number of relevant data points
        x_clusters = []; px_clusters = []  # empty lists
        breakpt = 0
        for i in range(1, nRelevant):
            if (xRelevant[i]-xRelevant[i-1]) > step_tolerance*avestep:
                x_clusters.append(xRelevant[breakpt:i])
                px_clusters.append(pxRelevant[breakpt:i])
                breakpt = i
        x_clusters.append(xRelevant[breakpt:])  # add final values
        px_clusters.append(pxRelevant[breakpt:])  # add final values

        # Outputs
        if verbose is True:
            nClusters=len(x_clusters)
            print('\tLowest value: {0:5f};\tHighest value: {1:5f}'.format(lowestValue, highestValue))
            print('\tNumber of clusters: {}'.format(nClusters))
            for n in range(nClusters):
                print('\t\tcluster {0}: {1:.2f}-{2:.2f}'.format(n, x_clusters[n].min(), x_clusters[n].max()))

        self.x = x; self.px = px  # PDF
        self.PxSort = PxSort
        self.mode = x[px==px.max()][0]  # most probable value
        self.xRelevant = xRelevant  # all x values within confidence range
        self.pxRelevant = pxRelevant  # all px values within confidence range
        self.lowestValue = lowestValue
        self.highestValue = highestValue
        self.nClusters = len(x_clusters)  # number of clusters
        self.x_clusters = x_clusters  # clusters of x values
        self.px_clusters = px_clusters  # clusters of px values


    ## Plot input data function
    def plotInput(self):
        nRelevant = len(self.xRelevant)
        fig, ax = plt.subplots()
        ax.plot(self.PxSort, color='k', linewidth=2)
        ax.set_xlabel('value'); ax.set_ylabel('integ\'d prob')
        ax.plot([nRelevant,nRelevant], [self.confidence/100,0],'b')
        ax.set_xticks([])
        ax.set_title('CDF')


    ## Plot output function
    def plotOutput(self, outName=None):
        x = np.pad(self.x, (1,1), 'edge')
        px = np.pad(self.px, (1,1), 'constant')
        fig, ax = plt.subplots()
        ax.fill(x, px,color=(0.4,0.4,0.4), zorder=2)
        # Plot clusters
        for i in range(self.nClusters):
            xHPD = self.x_clusters[i]
            xHPD = np.pad(xHPD, (1,1), 'edge')
            pxHPD = self.px_clusters[i]
            pxHPD = np.pad(pxHPD, (1,1), 'constant')
            ax.fill(xHPD, pxHPD, color=(0.3,0.3,0.6), zorder=3)
        ax.set_xlabel('value'); ax.set_ylabel('ref prob')
        ax.set_title('PDF')
        if outName:
            Fpdff.savefig(outName, dpi=600)



### PLOTTING FUNCTIONS ---
## Plot PDF as line
def plotPDF(x, px, title=None):
    Fpdf, axPDF = plt.subplots()
    axPDF.plot(x, px, color='k', linewidth=2)
    axPDF.set_xlabel('value'); axPDF.set_ylabel('rel prob')
    axPDF.set_title('PDF',title)
    return Fpdf, axPDF


## Plot PDF as filled polygon
def plotPDFfilled(x, px, title=None):
    xFill = np.pad(x, (1,1), 'edge')
    pxFill = np.pad(px, (1,1), 'constant')
    Fpdff, axPDFf = plt.subplots()
    axPDFf.fill(xFill, pxFill, color=(0.4,0.4,0.4), zorder=2)
    axPDFf.set_xlabel('value'); axPDFf.set_ylabel('ref prob')
    axPDFf.set_title('PDF',title)
    return Fpdff, axPDFf


## Plot cumulative density function
def plotCDF(x, Px, title=None):
    Fcdf, axCDF = plt.subplots()
    axCDF.plot(x, Px, color='k', linewidth=2)
    axCDF.set_xlabel('value'); axCDF.set_ylabel('integ\'d prob')
    axCDF.set_title('CDF',title)
    return Fcdf, axCDF



### SUPPORT FUNCTIONS ---
## Definition of gaussian functionn
def gauss(x, mu, sigma):
    y = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2);
    return y


## Gaussian kernel for filtering
def gauss_kernel(sigma):
    '''
    Return a kernel with Gaussian shape over 2-sigma range with
     width = sigma (samples)
    '''
    x = np.arange(-2*sigma, 2*sigma+1)
    y = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*(x/sigma)**2);
    return y


## Smooth PDF
def smoothPDF(x, px, ktype, kwidth):
    '''
    Smooth a PDF by convolution with a boxcar or Gaussian kernel.
    '''

    # Create kernel
    if ktype == 'boxcar':
        k = np.ones(kwidth)/kwidth
    elif ktype.lower() in ['gaussian','gauss']:
        k = gauss_kernel(kwidth)
    else:
        print('Choose from ktype = boxcar or gauss')

    # Smooth by convolution
    px = np.convolve(px, k, 'same')

    # Normalize area to 1.0
    P = np.trapz(px,x)
    px /= P

    return px


def cumulative_trapezoid(*args, **kwargs):
    """
    Function to allow cumulative trapezoid integration without updating
    scipy, or potentially having to upgrade python to 3.7
    """
    return sp_cumulative_trapezoid(*args, **kwargs)