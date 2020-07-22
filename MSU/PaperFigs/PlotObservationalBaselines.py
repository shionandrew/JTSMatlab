## Specifically for making prettier plots..
#
# Author: Shion Andrew
#

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
from astropy import stats

def plotObservationalBaselines1():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.ylabel(r'\textbf{\textit{G} Magnitude}')
    plt.xlabel(r'\textbf{\textit{G} Magnitude Error}')
    baselines = [50,100,200,300,400,500]

    for i in range(len(baselines)):
        getBaseline_photocount1(baselines[i])
        plt.clim(0,500)
    ax = plt.gca()
    fig = plt.gcf()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(.005))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.set_xlim(0,.04)
    ax.set_ylim(12,22)
    cbar = plt.colorbar()
    cbar.set_label(r'\textbf{Number of \textit{G} Observations}')
    plt.gca().invert_yaxis()
    plt.savefig('Figure3.pdf')
    plt.show()

def plotObservationalBaselines2():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.ylabel(r'\textbf{\textit{G} Magnitude}')
    plt.xlabel(r'\textbf{\textit{G} Magnitude Error}')
    baselines = [500,600,700,800,900]

    for i in range(len(baselines)):
        getBaseline_photocount1(baselines[i])
        plt.clim(500,900)
    ax = plt.gca()
    fig = plt.gcf()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(.005))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.set_xlim(0,.04)
    ax.set_ylim(12,22)
    cbar = plt.colorbar()
    cbar.set_label(r'\textbf{Number of \textit{G} Observations}')
    plt.gca().invert_yaxis()
    plt.savefig('Figure4.pdf')
    plt.show()

def getBaseline_photocount1(photocount):
    dtype = [('MagError', float), ('Magnitude', float), ('Ras', float), ('Decs', float), ('NumObs', int), ('SourceId', int)]
    values = []
    with open('/Users/touatokuchi/Desktop/MSU/PaperFigs/PhotObsG' + str(photocount) + '.csv') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                ra = float(row[0])
                dec = float(row[1])
                magnitude = float(row[2])
                flux = float(row[3])
                fluxError = float(row[4])
                numObs = int(row[5])
                sourceId = int(row[6])
                magError = 1.09*fluxError/flux
                values.append((magError, magnitude, ra, dec, numObs, sourceId))

            except:
                continue

    baselineData = np.array(values, dtype=dtype)
    plt.scatter(baselineData["MagError"], baselineData['Magnitude'], c=baselineData['NumObs'], cmap='plasma', s = .2)


def plotLowPhotObs():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.ylabel(r'\textbf{\textit{G} Magnitude}')
    plt.xlabel(r'\textbf{\textit{G} Magnitude Error}')

    dtype = [('MagError', float), ('Magnitude', float), ('Ras', float), ('Decs', float), ('NumObs', int), ('SourceId', int)]
    values = []
    with open('/Users/touatokuchi/Desktop/MSU/PaperFigs/LowestPhotObs.csv') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                ra = float(row[0])
                dec = float(row[1])
                magnitude = float(row[2])
                flux = float(row[3])
                fluxError = float(row[4])
                numObs = int(row[5])
                sourceId = int(row[6])
                magError = 1.09*fluxError/flux
                values.append((magError, magnitude, ra, dec, numObs, sourceId))

            except:
                continue

    baselineData = np.array(values, dtype=dtype)
    plt.scatter(baselineData["MagError"], baselineData['Magnitude'], color='black', s = .2)
    ax = plt.gca()
    fig = plt.gcf()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(.005))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.set_xlim(0,.04)
    ax.set_ylim(12,22)
    #cbar.set_label(r'\textbf{Number of \textit{G} Observations}')
    plt.gca().invert_yaxis()
    plt.savefig('Figure5.pdf')
    plt.show()


def main():
    plotObservationalBaselines1()
    plotObservationalBaselines2()
    plotLowPhotObs()
if __name__== "__main__":
    main()
