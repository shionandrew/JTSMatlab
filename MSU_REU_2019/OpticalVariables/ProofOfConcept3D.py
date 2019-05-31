#
# Author: Shion Andrew
#
# Optical Variables: plotting magnitude error as a function of magnitude and optical variance


from mpl_toolkits import mplot3d

import matplotlib.pyplot as plt
import csv
import os
import numpy as np

fig = plt.figure()
ax = plt.axes(projection='3d')


def plotStarsSouth(filename):
    x0 = [] #holds magError for stars of Vamp  < .2
    y0 = []

    x1 = [] #holds magError for stars of .2 < Vamp  < .4
    y1 = [] #holds magnitude for stars of .2 < Vamp < .4

    x2 = [] #holds stars of .3 < Vamp < .6
    y2 = [] #""

    x3 = [] #holds stars of Vamp > .6
    y3 = [] #""

    xRed = [] #holds stars in redback systems
    yRed = [] #""

    with open(filename) as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                #Calculate Error in Magnitude
                fluxError = float(row[35])
                flux = float(row[34])
                magError = 1.09*fluxError/flux
                Vamp = float(row[9])
                magnitude = float(row[36])


                if Vamp < .2:
        	   	  #append error in magnitude
                  x0.append(magError)
                  #Magnitude
                  y0.append(magnitude)

                elif Vamp >= .2 and Vamp < .4:
                    x1.append(magError)
                    y1.append(magnitude)

                elif Vamp >= .4 and Vamp < .6:
                          x2.append(magError)
                          y2.append(magnitude)
                elif Vamp >= .6:
                          x3.append(magError)
                          y3.append(magnitude)


                if magError > 0.05:
                    print(float(row[5]))
            except ValueError:
                  continue

        plt.scatter(x3,y3, s=9, color = 'blue')
        plt.scatter(x2,y2, s=2, color = 'yellow')
        plt.scatter(x1, y1, s=1, color = 'orange')
        plt.scatter(x0,y0, s = 2, color = 'red')
        plt.gca().invert_yaxis()
        plt.title("Catalina Southern Hemisphere Optical Variables")
        plt.xlabel('Magnitude Error')
        plt.ylabel('G Magnitude')
        plt.show()


def plotStarsNorth(filename):
    x0 = [] #holds magError for stars of Vamp  < .2
    y0 = []
    v0 = []

    x1 = [] #holds magError for stars of .2 < Vamp  < .4
    y1 = [] #holds magnitude for stars of .2 < Vamp < .4
    v1 = []

    x2 = [] #holds stars of .3 < Vamp < .6
    y2 = [] #""
    v2 = []

    x3 = [] #holds stars of Vamp > .6
    y3 = [] #""
    v3 = []

    xRed = [] #holds stars in redback systems
    yRed = [] #""

    with open(filename) as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                #Calculate Error in Magnitude
                fluxError = float(row[27])
                flux = float(row[26])
                magError = 1.09*fluxError/flux
                Vamp = float(row[6])
                magnitude = float(row[28])

                if Vamp < .12 or Vamp > 1 or magError > .035:
                    continue

                elif Vamp < .2:
        	   	  #append error in magnitude
                  x0.append(magError)
                  #Magnitude
                  y0.append(magnitude)
                  v0.append(Vamp)

                elif Vamp >= .2 and Vamp < .4:
                    x1.append(magError)
                    y1.append(magnitude)
                    v1.append(Vamp)

                elif Vamp >= .4 and Vamp < .6:
                          x2.append(magError)
                          y2.append(magnitude)
                          v2.append(Vamp)

                elif Vamp >= .6:
                          x3.append(magError)
                          y3.append(magnitude)
                          v3.append(Vamp)


                if magError > 0.05:
                    print(float(row[5]))

            except ValueError:
                  continue



        xdata = x0 + x1 + x2 + x3
        ydata = y0 + y1 + y2 + y3
        zdata = v0 + v1 + v2 + v3

        xdata = xdata[:1000]
        ydata = ydata[:1000]
        zdata = zdata[:1000]

        ax = plt.axes(projection='3d')
        ax.scatter3D( xdata, ydata, zdata, s=1, c=zdata, cmap='plasma');
        ax.set_xlabel('MagError')
        ax.set_ylabel('Magnitude')
        ax.view_init(elev=-177, azim=89)
        plt.show()

def main():
     plotStarsNorth("/Users/touatokuchi/Desktop/MSU_REU_2019/OpticalVariables/NorthStars.csv")

if __name__== "__main__":
  main()
