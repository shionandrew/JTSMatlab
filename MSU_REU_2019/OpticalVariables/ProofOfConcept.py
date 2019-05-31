#
# Author: Shion Andrew
#
# Optical Variables: plotting magnitude error as a function of magnitude and optical variance



import matplotlib.pyplot as plt
import csv
import os
import numpy as np


def plotStarsSouth(filename):
    v = []

    x0 = [] #holds magError for stars of Vamp  < .2
    y0 = []

    x1 = [] #holds magError for stars of .2 < Vamp  < .4
    y1 = [] #holds magnitude for stars of .2 < Vamp < .4

    x2 = [] #holds stars of .3 < Vamp < .6
    y2 = [] #""

    x3 = [] #holds stars of Vamp > .6
    y3 = [] #""

    x4 = []
    y4 = []

    x5 = []
    y5 = []

    x6 = []
    y6 = []

    x7 = []
    y7 = []

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

                #if magError >  .75:
                #    continue
                if Vamp < .2:
        	   	  #append error in magnitude
                  x0.append(magError)
                  #Magnitude
                  y0.append(magnitude)

                elif Vamp >= .2 and Vamp < .3:
                    x1.append(magError)
                    y1.append(magnitude)

                elif Vamp >= .3 and Vamp < .4:
                          x2.append(magError)
                          y2.append(magnitude)

                elif Vamp >= .4 and Vamp < .5:
                          x3.append(magError)
                          y3.append(magnitude)

                elif Vamp >= .5 and Vamp < .6:
                          x4.append(magError)
                          y4.append(magnitude)

                elif Vamp >= .6 and Vamp < .7:
                          x5.append(magError)
                          y5.append(magnitude)

                elif Vamp >= .7 and Vamp < .8:
                          x6.append(magError)
                          y6.append(magnitude)
                elif Vamp >= .8:
                          x7.append(magError)
                          y7.append(magnitude)


                #if magError > 0.05:
                #    print(float(row[5]))

            except ValueError:
                  continue



        x0av = sum(x0)/len(x0)
        y0av = sum(y0)/len(y0)

        x1av = sum(x1)/len(x1)
        y1av = sum(y1)/len(y1)

        x2av = sum(x2)/len(x2)
        y2av = sum(y2)/len(y2)

        x3av = sum(x3)/len(x3)
        y3av = sum(y3)/len(y3)

        x4av = sum(x4)/len(x4)
        y4av = sum(y4)/len(y4)

        x5av = sum(x5)/len(x5)
        y5av = sum(y5)/len(y5)

        x6av = sum(x6)/len(x6)
        y6av = sum(y6)/len(y6)

        x7av = sum(x7)/len(x7)
        y7av = sum(y7)/len(y7)

        plt.scatter(x7av,y7av, s=9, color = 'black')
        plt.scatter(x6av,y6av, s=9, color = 'purple')
        plt.scatter(x5av,y5av, s=9, color = 'blue')
        plt.scatter(x4av,y4av, s=9, color = 'Aqua')
        plt.scatter(x3av,y3av, s=9, color = 'green')
        plt.scatter(x2av,y2av, s=10, color = 'yellow')
        plt.scatter(x1av, y1av, s=10, color = 'orange')
        plt.scatter(x0av,y0av, s = 10, color = 'red')

#        plt.scatter(x7, y7, s=1, color = 'black')
#        plt.scatter(x6,y6, s=1, color = 'purple')
#        plt.scatter(x5,y5, s=1, color = 'blue')
#        plt.scatter(x4,y4, s=1, color = 'Aqua')
#        plt.scatter(x3,y3, s=1, color = 'green')
#        plt.scatter(x2,y2, s=1, color = 'yellow')
#        plt.scatter(x1,y1, s=1, color = 'orange')
#        plt.scatter(x0,y0, s = 1, color = 'red')

#        plt.scatter(np.log10(x7), np.log10(y7), s=1, color = 'black')
#        plt.scatter(np.log10(x6), np.log10(y6), s=1, color = 'purple')
#        plt.scatter(np.log10(x4), np.log10(y4), s=1, color = 'Aqua')
#        plt.scatter(np.log10(x3), np.log10(y3), s=1, color = 'green')
#        plt.scatter(np.log10(x2), np.log10(y2), s=1, color = 'yellow')
#        plt.scatter(np.log10(x1), np.log10(y1), s=1, color = 'orange')
#        plt.scatter(np.log10(x0), np.log10(y0), s = 1, color = 'red')

#        plt.scatter(np.log2(x7), y7, s=1, color = 'black')
#        plt.scatter(np.log2(x6), y6, s=1, color = 'purple')
#        plt.scatter(np.log2(x5), y5, s=1, color = 'blue')
#        plt.scatter(np.log2(x4), y4, s=1, color = 'Aqua')
#        plt.scatter(np.log2(x3), y3, s=1, color = 'green')
#        plt.scatter(np.log2(x2), y2, s=1, color = 'yellow')
#        plt.scatter(np.log2(x1), y1, s=1, color = 'orange')
#        plt.scatter(np.log2(x0), y0, s = 1, color = 'red')

#        plt.scatter(x7, np.log10(y7), s=1, color = 'black')
#        plt.scatter(x6, np.log10(y6), s=1, color = 'purple')
#        plt.scatter(x5, np.log10(y5), s=1, color = 'blue')
#        plt.scatter(x4, np.log10(y4), s=1, color = 'Aqua')
#        plt.scatter(x3, np.log10(y3), s=1, color = 'green')
#        plt.scatter(x2, np.log10(y2), s=1, color = 'yellow')
#        plt.scatter(x1, np.log10(y1), s=1, color = 'orange')
#    plt.scatter(x0, np.log10(y0), s = 1, color = 'red')


        plt.title("Catalina Southern Hemisphere Optical Variables")
        plt.xlabel('Mag Error')
        plt.ylabel('Magnitude')
        plt.show()


def plotStarsNorth(filename):
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
                fluxError = float(row[27])
                flux = float(row[26])
                magError = 1.09*fluxError/flux
                Vamp = float(row[6])
                magnitude = float(row[28])


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
        plt.title("Catalina Northern Hemisphere Optical Variables")
        plt.xlabel('Magnitude Error')
        plt.ylabel('G Magnitude')
        plt.show()

def main():
     plotStarsSouth("/Users/touatokuchi/Desktop/MSU_REU_2019/OpticalVariables/SouthStars.csv")
     plotStarsNorth("/Users/touatokuchi/Desktop/MSU_REU_2019/OpticalVariables/NorthStars.csv")

if __name__== "__main__":
  main()
