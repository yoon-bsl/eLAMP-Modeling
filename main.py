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

    def __init__(self, dataFile, Dt, bplength, type, eq):

        self.dataFile = dataFile

        self.Dt      = int(Dt)
        self.length  = int(bplength)

        if type.lower() not in ['avg', 'med']:
            print('"Type" input parameter must be avg or med')
            quit()
        else:
            self.type = type.lower()

        if eq.lower() not in ['exp', 'log']:
            print('"Equation" input parameter must be exp or log')
            quit()
        else:
            self.eq = eq.lower()

        if self.eq == 'exp':
            self.runExp()
        else:
            self.runLog()

    def runExp(self):
        self.dropletSizeAnalysis()

        if self.type == 'avg':
            self.initial = 1 / self.avgVol
        elif self.type == 'med':
            self.initial = 1 / self.medVol

        print(f'Initial concentration (copy per cm^3): {self.initial}')

        self.k = self.calculateGrowthConstant(self.Dt)
        self.diff = self.calcaulteDiffusivity(self.length)

        print(f'Growth constant (seconds^-1): {self.k}')
        print(f'Diffusivity (cm^2 / s): {self.diff}')

        saturated = 0.0
        t = 0
        saturationLevels = []
        time = []
        print('----------Modeling----------')
        while saturated < 100:
            print(f'Time (seconds): {t}')
            conc = self.calculateExpDiffusion(t, 
                                            self.initial, 
                                            self.k, 
                                            self.diff)
            print(f'Concentration: {conc}')

            amplicons = self.convertConcToAmplicons(conc, self.type)
            print(f'Number of amplicons: {amplicons}')
            
            SA = self.convertAmpliconToArea(amplicons, self.length)
            print(f'Surface area of amplicons: {SA}')

            saturated = self.calculateSaturation(SA, self.type)
            print(f'Current saturation: {saturated}')

            saturationLevels.append(saturated)
            time.append(t/60)
    
            t += 1
            if t > 600:
                break
            print('-----------')
        
        plt.plot(time[:-1], saturationLevels[:-1])
        plt.xlabel('Time (min)')
        plt.ylabel('Droplet Saturation (%)')
        plt.show()

    def runLog(self):
        '''
        '''
        pass

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

        diams["Volume (cm^3)"] = round(volume, 20)

        # Create a data column representing the corresponding surface area
        # of each droplet from the diameter
        nanoRadius = (diams["Diameter (um)"] / 2) * 1000    # get radius in nm
        surfaceArea = 4 * (nanoRadius ** 2) * math.pi       
        # surface area (nm^3)

        diams["Surface Area (nm^2)"] = round(surfaceArea, 3)

        # Print the various statistics:
        avgVol = statistics.mean(diams["Volume (cm^3)"].tolist())
        avgSA  = statistics.mean(diams["Surface Area (nm^2)"].tolist())
        medVol = statistics.median(diams["Volume (cm^3)"].tolist())
        medSA  = statistics.median(diams["Surface Area (nm^2)"].tolist())
        print(f'Average volume of measured droplets (in mL):         {avgVol}')
        print(f'Median volume of measured droplets (in mL):          {medVol}')
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

    def convertAmpliconToArea(self, amplicons, bp):
        '''
        Convert the base pair length of amplicons to hydrophobic area in nm^2
        '''
        mw = (2 * (bp * 607.4)) + 157.9 # MW in Da
        radius = 0.066 * (mw ** (1/3))  # Radius in nm
        area = math.pi * (radius ** 2)  # Hydrophobic area in nm^2
        return (area * amplicons)

    def convertConcToAmplicons(self, conc, type):
        '''
        Given the concentration of amplicons per unit area, will convert it to
        # of amplicons
        '''
        if type == 'avg':
            amplicons = (conc * (10**-14)) * self.avgSA
        elif type == 'med':
            amplicons = (conc * (10**-14)) * self.medSA
        return amplicons

    def calculateSaturation(self, SA, type):
        '''
        '''
        if type == 'avg':
            return (SA / self.avgSA)
        elif type == 'med':
            return (SA / self.medSA)

    def calculateExpDiffusion(self, t, initial, k, D):
        '''
        Exponential amplicon creation/diffusion calculation. Uses same model as 
        in Ulep et al. published in July 2019 in Scientific Reports. 
        Combines Fick's diffusion equation with exponential growth of LAMP 
        amplicons to model amplicon adsorption to the oil-water interface.
        '''
        return (2 * initial) * (math.exp(k * t)) * math.sqrt((D*t)/math.pi)

    def calculateLogDiffusion(self, t, initial, c, a, D, b=0.6):
        '''
        Logistic amplicon creation/diffusion calculation. Combines Fick's 
        diffusion equation with logistic growth model of LAMP amplicon creation
        to model amplicon adsorption to the oil-water interface.
        '''
        log = (2 * c * initial) / (1 + a * (math.exp(-b * t)))
        diff = math.sqrt((D * t) / math.pi)
        return (log * diff)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Model Emulsion LAMP amplicon creation and adsorption'
    )

    parser.add_argument('-i', '--input', required=False, 
                        default='dropletDiameters.csv',
                        help='Droplet data file name')
    parser.add_argument('-d', '--doublingTime', required=False,
                    default='24',
                    help='Doubling time of the LAMP reaction (in seconds)')
    parser.add_argument('-l', '--length', required=True,
                    help='Base pair length of the primary target')
    parser.add_argument('-t', '--type', required=True,
                help='Type of statistic to use for droplet size (avg or med)')
    parser.add_argument('-e', '--equation', required=True,
            help='Type of equation used for diffusion modeling (exp or log)')

    args = parser.parse_args()

    model = Model(
        args.input,
        args.doublingTime,
        args.length,
        args.type,
        args.equation     
        )
