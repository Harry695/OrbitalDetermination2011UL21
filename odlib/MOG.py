import numpy as np
from odlib.mathUtils import magnitude
from odlib.constants import Constants

def getFAndGConstants(tau1, tau3, posVec2, velVec2):
    posMag2 = magnitude(posVec2)
    u = 1 / posMag2**3.
    z = np.dot(posVec2, velVec2) / posMag2**2.
    q = np.dot(velVec2, velVec2) / posMag2**2. - u

    f1 = 1 - (u * tau1**2) / 2. + (u * z * tau1**3) / 2. + (3 * u * q - 15 * u * z**2 + u**2) * tau1**4 / 24.
    g1 = tau1 - (u * tau1**3) / 6. + (u * z * tau1**4) / 4.
    f3 = 1 - (u * tau3**2) / 2. + (u * z * tau3**3) / 2. + (3 * u * q - 15 * u * z**2 + u**2) * tau3**4 / 24.
    g3 = tau3 - (u * tau3**3) / 6. + (u * z * tau3**4) / 4.
    return f1, f3, g1, g3

def getcConstants(f1, f3, g1, g3):
    c1 = g3 / (f1 * g3 - g1 * f3)
    c3 = -g1 / (f1 * g3 - g1 * f3)
    return c1, c3

def getdConstants(f1, f3, g1, g3):
    d1 = -f3 / (f1 * g3 - f3 * g1)
    d3 = f1 / (f1 * g3 - f3 * g1)
    return d1, d3

def getScalarEquationDConstants(p1Dir, p2Dir, p3Dir, earthSunVector):
    """
    WARNING: Returns D contants for the same earth-sun vector. Probably need to transpose matrix.
    """
    Dj1 = np.dot(np.cross(earthSunVector, p2Dir), p3Dir)
    Dj2 = np.dot(np.cross(p1Dir, earthSunVector), p3Dir)
    Dj3 = np.dot(p1Dir, np.cross(p2Dir, earthSunVector))
    return [Dj1, Dj2, Dj3]

def getD0(p1Dir, p2Dir, p3Dir):
    return np.dot(p1Dir, np.cross(p2Dir, p3Dir))

def getDistances(cArr, D0, DArr): # c2 should be -1
    """
    cArr = [c1, c2, c3]
    DArr = [[d11, d12, d13],
            [d21, d22, d23],
            [d31, d32, d33]]
    """
    rhoList = []
    for i in range(3):
        rhoList.append(cArr[0] * DArr[i][0] + cArr[1] * DArr[i][1] + cArr[2] * DArr[i][2]) / (cArr[i] * D0)
    return rhoList

def lightSpeedCorrection(rhoMag, t):
    return t - rhoMag / Constants.LIGHT_SPEED_M_PER_S

def getRhoDirection(ra, dec):
    return np.array([[np.cos(ra) * np.cos(dec)],
                     [np.sin(ra) * np.cos(dec)],
                     [np.sin(dec)]])

def getTauArr(t1_Jd, t2_Jd, t3_Jd):
    """
    Parameter units are in Julian days!!!
    Returns an array of [tau1, tau3, tau0]
    """
    return np.array([t1_Jd - t2_Jd, t3_Jd - t2_Jd, t3_Jd - t1_Jd]) / Constants.DAY_IN_GAUSSIAN_DAY # convert to Gd
