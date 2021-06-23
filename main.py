import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
import numpy as np
import scipy as sp
import statistics
import argparse
import random
import math

class Model(object):
    '''
    Main object that will conduct all modeling, calculating, and plotting
    of the amplicon formation and adsorption to the oil-water interface during
    emulsion LAMP reactions

    input args:

    '''

    def __init__(self):
        self.setSeed()

    def setSeed(self):
        random.seed(1234)

if __name__ == '__main__':
    # Set the random seed to make data reproducable
    random.seed(1234)

    # Import droplet diameters from non-amplified emulsion droplets
    # containing LAMP reaction mixtures
    diams = pd.read_csv('dropletDiameters.csv', delimiter='\n')
    diams = pd.DataFrame(diams)

    # Create data column representing the corresponding volume of each
    # droplet from the diameter
    microRadius = (diams["Diameter (um)"] / 2)       # get radius in um
    standardRadius = microRadius / (10 ** 6)    # get radius in meters
    centiRadius = standardRadius * 100          # radius in centimeters
    volume = (4/3)*(centiRadius**3)*math.pi     # volume in cm^3 or mL
    picoVolume = volume * (10**9)               # volume in pL

    diams["Volume (pL)"] = round(picoVolume, 3)

    # Create a data column representing the corresponding surface area
    # of each droplet from the diameter
    nanoRadius = (diams["Diameter (um)"] / 2) * 1000    # get radius in nm
    surfaceArea = 4 * (nanoRadius ** 2) * math.pi       # surface area (nm^3)

    diams["Surface Area (nm^3)"] = round(surfaceArea, 3)

    # Tell the user the average and median volume/surface area of the droplets
    avgVol = statistics.mean(diams["Volume (pL)"].tolist())
    avgSA  = statistics.mean(diams["Surface Area (nm^3)"].tolist())
    medVol = statistics.median(diams["Volume (pL)"].tolist())
    medSA  = statistics.median(diams["Surface Area (nm^3)"].tolist())
    print(f'Average volume of measured droplets (in pL):         {avgVol}')
    print(f'Median volume of measured droplets (in pL):          {medVol}')
    print(f'Average surface area of measured droplets (in nm^3): {avgSA}')
    print(f'Median surface area of measured droplets (in nm^3):  {medSA}')

    # sb.set_theme(style="darkgrid")

    # plt.figure()
    # ax = sb.boxplot(x=temp)
    # plt.show()

    parser = argparse.ArgumentParser(description='Input for eLAMP monitoring')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')