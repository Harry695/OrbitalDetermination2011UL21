import numpy as np
from math import *
from odlib.mathUtils import *
from odlib.constants import Constants


def getAngularMomentum(posVec, velVec):
    return np.cross(posVec, velVec)


def getSemimajorAxis(posVec, velVec):  # works
    posMag = magnitude(posVec)
    return 1 / (2 / posMag - np.dot(velVec, velVec) / Constants.MU)


def getEccentricity(posVec, velVec):  # works
    a = getSemimajorAxis(posVec, velVec)
    return sqrt(1 - magnitude(np.cross(posVec, velVec)) ** 2 / (Constants.MU * a))


def getInclination(posVec, velVec):  # works
    """
    Returns inclination in radians.
    """
    hVec = getAngularMomentum(posVec, velVec)
    return acos(hVec[2] / magnitude(hVec))


def getLongitudeOfAscendingNode(posVec, velVec):  # works
    """
    Returns longitude of ascending node in radians.
    """
    hVec = getAngularMomentum(posVec, velVec)
    hMag = magnitude(hVec)
    i = getInclination(posVec, velVec)
    return radians(quadrantCorrection(hVec[0] / (hMag * sin(i)), -hVec[1] / (hMag * sin(i))))


def getArgumentOfPerihelion(posVec, velVec):  # works
    """
    Returns argument of perihelion in radians.
    """
    posMag = magnitude(posVec)
    a = getSemimajorAxis(posVec, velVec)
    e = getEccentricity(posVec, velVec)
    i = getInclination(posVec, velVec)
    omega = getLongitudeOfAscendingNode(posVec, velVec)
    hVec = getAngularMomentum(posVec, velVec)
    hMag = magnitude(hVec)
    sinU = posVec[2] / (posMag * sin(i))
    # print(sinU)
    cosU = (posVec[0] * cos(omega) + posVec[1] * sin(omega)) / posMag
    # print(cosU)
    u = radians(quadrantCorrection(sinU, cosU))
    sinMu = a * (1 - e**2) * np.dot(posVec, velVec) / (e * hMag * posMag)
    # print(sinMu)
    cosMu = (a * (1 - e**2) / posMag - 1) / e
    # print(cosMu)
    # mu = quadrantCorrection(sinMu, cosMu)
    mu = getTrueAnomaly(posVec, velVec)
    # print("U", u)
    # print("Mu", mu)
    w = u - mu
    if w < 0:
        w += 2 * np.pi
    return w


def getTrueAnomaly(posVec, velVec):
    """
    Returns true anomaly in radians.
    """
    posMag = magnitude(posVec)
    h = getAngularMomentum(posVec, velVec)
    hMag = magnitude(h)
    a = getSemimajorAxis(posVec, velVec)
    e = getEccentricity(posVec, velVec)
    sinMu = a * (1 - e**2) * np.dot(posVec, velVec) / (e * hMag * posMag)
    # print(sinMu)
    cosMu = (a * (1 - e**2) / posMag - 1) / e
    # print(cosMu)
    return radians(quadrantCorrection(sinMu, cosMu))


def getMeanAnomaly(posVec, velVec):
    """
    Returns mean anomaly in radians.
    """
    posMag = magnitude(posVec)
    a = getSemimajorAxis(posVec, velVec)
    e = getEccentricity(posVec, velVec)
    mu = getTrueAnomaly(posVec, velVec)
    E = acos((1 - posMag / a) / e) if mu < np.pi else 2 * np.pi - acos((1 - posMag / a) / e)
    # print("E at observation", E)
    # print("M at observation", E - e * sin(E))
    return E - e * sin(E)


def getEccentricAnomaly(M, ecc):
    return solveUsingNewtonsMethod(
        M, lambda E: E - ecc * sin(E) - M, lambda E: 1 - ecc * cos(E)
    )


def getLastPerihelionTime(a, M_rad, t_Gd):
    # print("M initial", M_rad)
    # print("last perihelion", t_Gd - M_rad * (a**1.5))
    return t_Gd - M_rad * (a**1.5)
