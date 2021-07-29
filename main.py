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

    def __init__(self, Dt, initialC, bplength, type, 
                            dataFile='dropletDiameters.csv'):

        self.dataFile = dataFile

        self.Dt      = Dt
        self.initial = initialC
        self.length  = bplength

        if type.lower() not in ['avg', 'med']:
            print('"Type" input parameter must be avg or med')
            quit()
        else:
            self.type = type

        self.run()

    def run(self):
        self.dropletSizeAnalysis()

        self.k = self.calculateGrowthConstant(self.Dt)
        self.diff = self.calcaulteDiffusivity(self.length)

        saturated = 0.0
        t = 0
        while saturated < 100:
            conc = self.calculateDiffusion(t, self.initial, self.k, self.diff)
            # TODO: calculate # of amplicons given above concentration over area
            # TODO: calculate SA of the total number of amplicons above

            if self.type == 'avg':
                pass
                # TODO: calculate saturation of the interface (% coverage)
            else:
                pass
                # TODO: calculate saturation of the interface (% coverage)
            # TODO: save the saturation to a list for plotting
        
        # TODO: plot the saturation over time 

    def setSeed(self):
        '''
        If needed, will set the rng seed to ensure reproducibility of the
        model's results over multiple runs
        '''
        random.seed(1234)

    def dropletSizeAnalysis(self):
        '''
        Method takes the given droplet diameter data file and calculates the 
        average and median volumes/surface areas of the droplets for use in
        the model's given parameters
        '''
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

        diams["Surface Area (nm^2)"] = round(surfaceArea, 3)

        # Print the various statistics:
        avgVol = statistics.mean(diams["Volume (pL)"].tolist())
        avgSA  = statistics.mean(diams["Surface Area (nm^3)"].tolist())
        medVol = statistics.median(diams["Volume (pL)"].tolist())
        medSA  = statistics.median(diams["Surface Area (nm^3)"].tolist())
        print(f'Average volume of measured droplets (in pL):         {avgVol}')
        print(f'Median volume of measured droplets (in pL):          {medVol}')
        print(f'Average surface area of measured droplets (in nm^2): {avgSA}')
        print(f'Median surface area of measured droplets (in nm^2):  {medSA}')

        self.avgSA  = avgSA
        self.medSA  = medSA
        self.avgVol = avgVol
        self.medVol = medVol

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
        Given the base pair length of an amplicon, will return the diffusivity
        of said amplicon. Based on Lukacs et al. from January 2000 in the
        Journal of Biological Chemistry.

        input:
        bp = base pair length of the amplicons

        output:
        diff = diffusivity of the amplicon (in cm^2/s)
        '''
        return ( 4.9 * ( 10 ** -6 ) ) * ( bp ** -0.72 )

    def calculateDiffusion(self, t, initial, k, D):
        '''
        Main amplicon creation/diffusion calculation. Uses same model as in
        Ulep et al. published in July 2019 in Scientific Reports. 
        Combines Fick's diffusion equation with exponential growth of LAMP 
        amplicons to model amplicon adsorption to the oil-water interface.
        '''
        return (2 * initial) * (math.exp(k * t )) * math.sqrt((D*t)/math.pi)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Model Emulsion LAMP amplicon creation and adsorption'
    )

    parser.add_argument('-i', '--input', required=False, 
                        help='Droplet data file name')
    parser.add_argument('-d', '--doublingTime', required=True,
                    help='Doubling time of the LAMP reaction (in seconds)')
    parser.add_argument('-c', '--concentration', required=True,
                    help='Initial concentration of target (in copy #/cm^3)')
    parser.add_argument('-l', '--length', required=True,
                    help='Base pair length of the primary target')
    parser.add_argument('-t', '--type', required=True,
                help='Type of statistic to use for droplet size (avg or med)')

    args = parser.parse_args()

    model = Model(
        args.doublingTime,
        args.concentration,
        args.length,
        args.type,
        dataFile=args.input       
        )
