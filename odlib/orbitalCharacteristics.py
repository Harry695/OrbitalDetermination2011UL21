import numpy as np
from math import *
from odlib.mathUtils import *
from odlib.constants import Constants

def getAngularMomentum(posVec, velVec):
    return np.cross(posVec, velVec)

def getSemimajorAxis(posVec, velVec): # works
    posMag = magnitude(posVec)
    return 1 / (2 / posMag - np.dot(velVec, velVec) / Constants.MU)

def getEccentricity(posVec, velVec): # works
    a = getSemimajorAxis(posVec, velVec)
    return sqrt(1 - magnitude(np.cross(posVec, velVec)) ** 2 / (Constants.MU * a))

def getInclination(posVec, velVec): # no worky
    hVec = getAngularMomentum(posVec, velVec)
    return acos(hVec[2] / magnitude(hVec))

def getLongitudeOfAscendingNode(posVec, velVec): #work?
    hVec = getAngularMomentum(posVec, velVec)
    hMag = magnitude(hVec)
    i = getInclination(posVec, velVec)
    return quadrantCorrection(hVec[0] / (hMag * sin(i)), -hVec[1] / (hMag * sin(i)))

def getArgumentOfPerihelion(posVec, velVec): # no worky
    posMag = magnitude(posVec)
    a = getSemimajorAxis(posVec, velVec)
    e = getEccentricity(posVec, velVec)
    i = getInclination(posVec, velVec)
    omega = getLongitudeOfAscendingNode(posVec, velVec)
    hVec = getAngularMomentum(posVec, velVec)
    hMag = magnitude(hVec)
    sinU = posVec[2] / (posMag * sin(i))
    cosU = (posVec[0] * cos(omega) + posVec[1] * sin(omega)) / hMag
    u = quadrantCorrection(sinU, cosU)
    sinMu = a * (1 - e ** 2) * np.dot(posVec, velVec) / (e * hMag * posMag)
    cosMu = (a * (1 - e ** 2) / hMag - 1) / e
    mu = quadrantCorrection(sinMu, cosMu)
    return u - mu