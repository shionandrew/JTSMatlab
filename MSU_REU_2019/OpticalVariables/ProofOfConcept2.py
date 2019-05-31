#
# Author: Shion Andrew
#
# Optical Variables: plotting magnitude error as a function of magnitude and optical variance



import matplotlib.pyplot as plt
import csv
import os


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

                if magnitude < 10:
                    y0.append(Vamp)
                    x0.append(magError)

                elif magnitude < 11:
                    y1.append(Vamp)
                    x1.append(magError)
                elif magnitude < 12:
                    y2.append(Vamp)
                    x2.append(magError)

                elif magnitude < 13:
                    y3.append(Vamp)
                    x3.append(magError)

            except ValueError:
                  continue


        x0av = sum(x0)/len(x0)
        y0av = sum(y0)/len(y0)

        x1av = sum(x1)/len(x1)
        y1av = sum(y1)/len(y1)

        x2av = sum(x2)/len(x2)
        y2av = sum(y2)/len(y2)

        x3av = sum(x2)/len(x2)
        y3av = sum(y2)/len(y2)

        plt.scatter(x3av,y3av, s=9, color = 'blue')
        plt.scatter(x2av,y2av, s=10, color = 'yellow')
        plt.scatter(x1av, y1av, s=10, color = 'orange')
        plt.scatter(x0av,y0av, s = 10, color = 'red')
        #plt.gca().invert_yaxis()
        plt.title("Catalina Southern Hemisphere Optical Variables")
        plt.xlabel('Magnitude Error')
        plt.ylabel('V amp')
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
     #plotStarsNorth("/Users/touatokuchi/Desktop/MSU_REU_2019/OpticalVariables/NorthStars.csv")

if __name__== "__main__":
  main()
