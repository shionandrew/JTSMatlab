## Specifically for final paper
# Author: Shion Andrew

import math
import statistics
import csv
import os

# MatPlotlib
from matplotlib import pylab
#import plotly.plotly as py
#import plotly.graph_objs as go
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib import rc

# Scientific libraries
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline

# Astropy and Gaia
'''from astroquery.gaia import Gaia
import astroquery
import keyring
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord'''
from astropy import stats

def plotCatalinaAmplitude():
    # plot baselinecurve
    dtypeBaseline = [('GMagnitude', float), ('RPMagnitude', float), ('BPMagnitude', float), ('GMagError', float), ('RPMagError', float), ('BPMagError', float)]
    finalValuesBaseline = []
    with open('/Users/touatokuchi/Desktop/MSU/KnownVariablesTest/RandomSample.csv') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                Gmagnitude = float(row[0])
                RPmagnitude = float(row[1])
                BPmagnitude = float(row[2])

                Gflux = float(row[3])
                GfluxError = float(row[4])
                GmagError = 1.09*GfluxError/Gflux

                RPflux = float(row[5])
                RPfluxError = float(row[6])
                RPmagError = 1.09*RPfluxError/RPflux

                BPflux = float(row[7])
                BPfluxError = float(row[8])
                BPmagError = 1.09*BPfluxError/BPflux

                finalValuesBaseline.append((Gmagnitude, RPmagnitude, BPmagnitude, GmagError, RPmagError, BPmagError))
            except ValueError:
                continue

    sampleStars = np.array(finalValuesBaseline, dtype = dtypeBaseline)

    dtype = [('Ras', float), ('Decs', float), ('period', float), ('Amplitude', float), ('GMagnitude', float), ('RPMagnitude', object), ('BPMagnitude', object), ('GMagError', float), ('RPMagError', object), ('BPMagError', object), ('SourceId', int)]
    finalValues = []

    totalCount = 0
    count = 0
    with open('/Users/touatokuchi/Desktop/MSU/KnownVariablesTest/CatalinaFull.csv') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                ra = float(row[0])
                dec = float(row[1])
                period = float(row[2])
                Vamp = float(row[3])
                sourceId = int(row[4])
                totalCount +=1

                Gflux = float(row[5])
                GfluxError = float(row[6])
                Gmagnitude = float(row[7])
                GmagError = 1.09*GfluxError/Gflux
                if Gmagnitude > 14 and Gmagnitude < 19.5:
                    count+=1
                BPflux = row[8]
                BPfluxError = row[9]
                BPmagnitude = row[10]
                try:
                    BPmagError = 1.09*float(BPfluxError)/float(BPflux)
                except:
                    BPmagError = 0

                RPflux = row[11]
                RPfluxError = row[12]
                RPmagnitude = row[13]
                try:
                    RPmagError = 1.09*float(RPfluxError)/float(RPflux)
                except:
                    RPMagError = 0

                finalValues.append((ra, dec, period, Vamp, Gmagnitude, RPmagnitude, BPmagnitude, GmagError, RPmagError, BPmagError, sourceId))

            except:
                  print(row[0])

    print(totalCount)
    print(count)
    print(len(finalValues))
    variableStars = np.array(finalValues, dtype=dtype)

    ## plot averages of each bin
    ampBin1 = variableStars[np.where(variableStars['Amplitude'] <= .2)]
    medAbsDev1 = stats.median_absolute_deviation(ampBin1['GMagError'])
    stdevMagError1 = 1.4826*medAbsDev1

    ampBin2 = variableStars[np.where((variableStars['Amplitude'] <= .4) & (variableStars['Amplitude'] >= .2 ))]
    medAbsDev2 = stats.median_absolute_deviation(ampBin2['GMagError'])
    stdevMagError2 = 1.4826*medAbsDev2

    ampBin3 = variableStars[np.where((variableStars['Amplitude'] <= .6) & (variableStars['Amplitude'] >= .4 ))]
    medAbsDev3 = stats.median_absolute_deviation(ampBin3['GMagError'])
    stdevMagError3 = 1.4826*medAbsDev3

    ampBin4 = variableStars[np.where((variableStars['Amplitude'] <= .8) & (variableStars['Amplitude'] >= .6 ))]
    medAbsDev4 = stats.median_absolute_deviation(ampBin4['GMagError'])
    stdevMagError4 = 1.4826*medAbsDev4

    ampBin5 = variableStars[np.where((variableStars['Amplitude'] <= 1) & (variableStars['Amplitude'] >= .8 ))]
    medAbsDev5 = stats.median_absolute_deviation(ampBin5['GMagError'])
    stdevMagError5 = 1.4826*medAbsDev5

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.ylabel(r'\textbf{\textit{G} Magnitude}')
    plt.xlabel(r'\textbf{\textit{G} Magnitude Error}')

    plt.axvline(x=statistics.mean(ampBin1['GMagError']), color = 'darkblue', ls = '--')
    plt.axvline(x=statistics.mean(ampBin2['GMagError']), color = 'blue', ls = '--')
    plt.axvline(x=statistics.mean(ampBin3['GMagError']), color = 'cyan', ls = '--')
    plt.axvline(x=statistics.mean(ampBin4['GMagError']), color = 'yellow', ls = '--')
    plt.axvline(x=statistics.mean(ampBin5['GMagError']), color = 'red', ls = '--')

    plt.errorbar(statistics.mean(ampBin1['GMagError']), statistics.mean(ampBin1['GMagnitude']), xerr=stdevMagError1, ecolor = 'black', capsize = 10)
    plt.errorbar(statistics.mean(ampBin2['GMagError']), statistics.mean(ampBin2['GMagnitude']), xerr=stdevMagError2, ecolor = 'black', capsize = 10)
    plt.errorbar(statistics.mean(ampBin3['GMagError']), statistics.mean(ampBin3['GMagnitude']), xerr=stdevMagError3, ecolor = 'black', capsize = 10)
    plt.errorbar(statistics.mean(ampBin4['GMagError']), statistics.mean(ampBin4['GMagnitude']), xerr=stdevMagError4, ecolor = 'black', capsize = 10)
    plt.errorbar(statistics.mean(ampBin5['GMagError']), statistics.mean(ampBin5['GMagnitude']), xerr=stdevMagError5, ecolor = 'black', capsize = 10)


    plt.scatter(sampleStars['GMagError'], sampleStars['GMagnitude'], color = 'black', s = 0.2)
    plt.scatter(variableStars['GMagError'], variableStars['GMagnitude'], c=variableStars['Amplitude'], cmap='jet', s = .2)
    cbar = plt.colorbar()
    cbar.set_label(r'\textbf{Amplitude (mag)}')
    ax = plt.gca()
    fig = plt.gcf()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(.005))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.set_xlim(0,.035)
    plt.gca().invert_yaxis()
    plt.clim(0, 1)
    plt.savefig('Figure2.pdf')
    plt.show()

    '''plt.ylabel(r'\textbf{BP mean magnitude}')
    plt.xlabel(r'\textbf{BP magnitude error}')

    plt.scatter(sampleStars['BPMagError'], sampleStars['BPMagnitude'], color = 'black', s = 0.2)
    plt.scatter(variableStars['BPMagError'], variableStars['BPMagnitude'], c=variableStars['Amplitude'], cmap='jet', s = .2)
    plt.colorbar()
    ax = plt.gca()
    fig = plt.gcf()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(.005))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.set_xlim(0,.09)
    plt.gca().invert_yaxis()
    plt.clim(0, 1)
    plt.show()'''


def main():
    plotCatalinaAmplitude()

if __name__== "__main__":
    main()
