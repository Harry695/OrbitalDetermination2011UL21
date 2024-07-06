import numpy as np

def getAngularMomentum(posVec, velVec):
    return np.cross(posVec, velVec)