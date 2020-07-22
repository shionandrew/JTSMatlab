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

# Astropy and Gaia
from astropy import stats

def plotRPObservationalBaselines():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.ylabel(r'\textbf{\textit{RP} Magnitude}')
    plt.xlabel(r'\textbf{\textit{RP} Magnitude Error}')
    #baselines = [0,50,100,200,300,400,500]
    baselines = [10,15, 25, 50, 75, 100]

    for i in range(len(baselines)):
        getRPBaseline_photocount1(baselines[i])

    ax = plt.gca()
    fig = plt.gcf()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(.005))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.set_xlim(0,.04)
    ax.set_ylim(12,21)
    cbar = plt.colorbar()
    cbar.set_label(r'\textbf{Number of \textit{RP} Observations}')
    plt.gca().invert_yaxis()
    plt.savefig('Figure9.pdf')
    plt.show()


def getRPBaseline_photocount1(photocount, band):
    dtype = [('MagError', float), ('Magnitude', float), ('Ras', float), ('Decs', float), ('NumObs', int), ('SourceId', int)]
    values = []
    with open('/Users/touatokuchi/Desktop/MSU/HPCC/RP_Band/PhotObsRP' + str(photocount) + '.csv') as csvfile:
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
    plt.clim(10, 100)



def plotBPObservationalBaselines():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.ylabel(r'\textbf{\textit{BP} Magnitude}')
    plt.xlabel(r'\textbf{\textit{BP} Magnitude Error}')
    #baselines = [0,50,100,200,300,400,500]
    baselines = [10, 20,30,40,50,60,70, 80, 90, 100]

    for i in range(len(baselines)):
        getBPBaseline_photocount1(baselines[i])

    ax = plt.gca()
    fig = plt.gcf()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(.005))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.set_xlim(0,.04)
    ax.set_ylim(12,21)
    cbar = plt.colorbar()
    cbar.set_label(r'\textbf{Number of \textit{BP} Observations}')
    plt.gca().invert_yaxis()
    plt.savefig('Figure10.pdf')
    plt.show()

def getBPBaseline_photocount1(photocount):
    dtype = [('MagError', float), ('Magnitude', float), ('NumObs', int)]
    values = []
    with open('/Users/touatokuchi/Desktop/MSU/HPCC/BP_Band/PhotObsBP' + str(photocount) + '.csv') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                magnitude = float(row[0])
                flux = float(row[1])
                fluxError = float(row[2])
                numObs = int(row[3])
                magError = 1.09*fluxError/flux
                values.append((magError, magnitude, numObs))

            except:
                continue

    baselineData = np.array(values, dtype=dtype)
    plt.scatter(baselineData["MagError"], baselineData['Magnitude'], c=baselineData['NumObs'], cmap='plasma', s = .2)
    plt.clim(10,100)

def main():
    #plotRPObservationalBaselines()
    plotBPObservationalBaselines()

if __name__== "__main__":
    main()
