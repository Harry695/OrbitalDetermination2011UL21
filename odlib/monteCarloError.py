import numpy as np
from odlib.OrbitalBody import OrbitalBody
import numpy.random as random
from math import *
from odlib.ObervationInput import ObservationInput


def getMonteCarloErrorBounds(obs, errorArr):
    """
    Parameters:
        obs: an array of 3 ObservationInput objects.
        errorArr: a 2D array of [raError, decError] in degrees for each observation.
    Returns a list of tuples of (mean, std) of each orbital element.
    """
    # random.seed(3)
    N = 1000

    mean = OrbitalBody.fromObservations(obs)
    meanOEs = mean.getAllOrbitalElements()
    
    distribution = []
    count = 0
    nonConvergentCount = 0
    for i in range(N):
        # if input("Continue?") != "":
        #     break
        count += 1
        print("Monte Carlo", count)

        # print("ra0 mean", degrees(obs[0].ra))
        ra0 = random.normal(obs[0].ra, errorArr[0][0])
        ra1 = random.normal(obs[1].ra, errorArr[1][0])
        ra2 = random.normal(obs[2].ra, errorArr[2][0])
        dec0 = random.normal(obs[0].dec, errorArr[0][1])
        dec1 = random.normal(obs[1].dec, errorArr[1][1])
        dec2 = random.normal(obs[2].dec, errorArr[2][1])
        
        # input ra & dec to normal distribution ones
        obsInput = [ObservationInput(obs[0].time_Jd, ra0, dec0, obs[0].earthSunVector),
                 ObservationInput(obs[1].time_Jd, ra1, dec1, obs[1].earthSunVector),
                 ObservationInput(obs[2].time_Jd, ra2, dec2, obs[2].earthSunVector)]
        
        # debug
        # for e in obsInput:
        #     print(e)

        # append associated orbital elements
        body = OrbitalBody.fromObservations(obsInput)
        if body == None:
            nonConvergentCount += 1
            continue
        # body.printAllOrbitalElements() # debug
        distribution.append(body.getAllOrbitalElements())
    
    print("count nonconvergent", nonConvergentCount)
    mean.printAllOrbitalElements()
    # finished generating, start analyzing
    distribution = np.transpose(distribution) # 2d array of [[a], [e], [i], [omega], [w]]
    oeStats = []
    for oe in distribution:
        oeStats.append((np.mean(oe), np.std(oe))) # get mean and std of each of the distribution of oe
    return oeStats