import math
import statistics
import csv
import os

# MatPlotlib
from matplotlib import pylab
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib import rc

# Scientific libraries
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline

def plotLightCurve(filename, P, t0, toplim, botlim, outputFilename):
    magerr = []
    mag = []
    abjd = []

    with open('/Users/touatokuchi/Desktop/MSU/PlotsForPaper/' + filename) as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                mag.append(float(row[1]))
                magerr.append(float(row[2]))
                abjd.append(float(row[5]))
            except:
                print("error")


    aphase = []
    for i in range(len(abjd)):
        aphase.append(foldAt(abjd[i],P,T0=t0,getEpoch=False))

    mag2 = []
    aphase2 = []
    magerr2 = []
    for i in range(len(aphase)):
        mag2.append(mag[i])
        magerr2.append(magerr[i])
        aphase2.append(aphase[i]+1)
    mag = mag + mag2
    aphase = aphase + aphase2
    magerr = magerr + magerr2

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.ylabel(r'\textbf{\textit{V} Magnitude}')
    plt.xlabel(r'\textbf{Phase}')
    #plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xticks(np.arange(0, 2.01
    , 0.5))
    plt.scatter(aphase, mag, color = 'red', s = .5)
    plt.errorbar(aphase, mag, yerr=magerr, ecolor = 'black', lw = .6, capsize = 1, fmt = 'none')
    #plt.yticks(np.arange(14.75,15.75,.25))
    plt.yticks(np.arange(toplim, botlim, abs(botlim-toplim)/5))
    #ax.scatter(aphase, mag, color = 'red', s = .5)
    #ax.errorbar(aphase, mag, yerr=magerr, ecolor = 'black', lw = .3, capsize = .5, fmt = 'none')


    plt.gca().invert_yaxis()
    plt.savefig(outputFilename)
    plt.show()


def foldAt(time, period, T0=0.0, getEpoch=False):
    epoch = np.floor( (time - T0)/period )
    phase = (time - T0)/period - epoch
    if getEpoch:
        return phase, epoch
    return phase

def round_of_rating(number):
    """Round a number to the closest half integer.
    >>> round_of_rating(1.3)
    1.5
    >>> round_of_rating(2.6)
    2.5
    >>> round_of_rating(3.0)
    3.0
    >>> round_of_rating(4.1)
    4.0"""

    return round(number * 2) / 2

def main():
    #Stripe Variable at 311.42801, 0.4186
    P = 1.0750575
    t0 = 53505.42526
    plotLightCurve('StripeVariable.csv', P, t0, toplim=17.2, botlim=18, outputFilename = 'Figure8.pdf')


if __name__== "__main__":
    main()
