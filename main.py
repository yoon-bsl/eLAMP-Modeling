import matplotlib.pyplot as plt
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

    def __init__(self, dataFile='dropletDiameters.csv'):
        self.dataFile = dataFile

        self.setSeed()

    def run(self):
        pass

    def setSeed(self):
        random.seed(1234)

    def dropletSizeAnalysis(self):
         # Import droplet diameters from non-amplified emulsion droplets
        # containing LAMP reaction mixtures
        diams = pd.read_csv(self.dataFile, delimiter='\n')
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
        surfaceArea = 4 * (nanoRadius ** 2) * math.pi       
        # surface area (nm^3)

        diams["Surface Area (nm^3)"] = round(surfaceArea, 3)

        # Print the various statistics:
        avgVol = statistics.mean(diams["Volume (pL)"].tolist())
        avgSA  = statistics.mean(diams["Surface Area (nm^3)"].tolist())
        medVol = statistics.median(diams["Volume (pL)"].tolist())
        medSA  = statistics.median(diams["Surface Area (nm^3)"].tolist())
        print(f'Average volume of measured droplets (in pL):         {avgVol}')
        print(f'Median volume of measured droplets (in pL):          {medVol}')
        print(f'Average surface area of measured droplets (in nm^3): {avgSA}')
        print(f'Median surface area of measured droplets (in nm^3):  {medSA}')

    def calculateGrowthConstant(self, Dt):
        '''
        Given a doubling time (in seconds) for LAMP, the function returns the
        growth constant for the reaciton

        input:
        Dt = doubling time of the reaction from literature/specifications (s)

        output:
        k  = growth constant (s^-1)
        '''
        return 1/Dt

    def calcaulteDiffusivity(self, bp):
        '''
        Given the base pair length of a 
        '''
        pass

    def calculate(self, initial, k, D):
        '''
        Main amplicon creation/diffusion calculation. Uses same model as in
        Ulep et al. published in July 2019. Combines Fick's diffusion 
        equation with exponential growth of LAMP amplicons to model amplicon
        adsorption to the oil-water interface.
        '''
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Model Emulsion LAMP amplicon creation and adsorption'
    )

    parser.add_argument('-i', '--data', required=False, 
                        help='Droplet data file name')
    # parser.add_argument('-c', '--column', required=True,
    #                 help='Image column containing sample colors')
    # parser.add_argument('-f', '--fristclass', nargs='+', required=True,
    #                 help='List containing sample #s of the first sample type')
    # parser.add_argument('-s', '--secondclass', nargs='+', required=True,
    #                 help='List containing sample #s of the second sample type')
    # parser.add_argument('-t', '--thirdclass', nargs='+', required=True,
    #                 help='List containing sample #s of the third sample type')
    # parser.add_argument('-l', '--lastclass', nargs='+', required=True,
    #                 help='List containing sample #s of the fourth sample type')

    args = parser.parse_args()

    model = Model(
        args.data       
        )
    model.run()
