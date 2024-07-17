from math import*
import numpy as np
from enum import Enum

class Axis(Enum):
    X = "X"
    Y = "Y"
    Z = "Z"

def magnitude(vec): # TODO: ACTUALLY TEST  
    vec = np.array(vec)
    vecSq = vec * vec
    return sqrt(np.sum(vecSq))

def quadrantCorrection(sin, cos): #TODO: test edge cases more
    thetaS = degrees(asin(sin))
    thetaC = degrees(acos(cos))
    if (sin > 0 and cos >= 0): #Q1
        return thetaC
    elif (sin > 0 and cos < 0): #Q2
        return thetaC
    elif (sin < 0 and cos < 0): #Q3
        return 180 - thetaS
    elif (sin < 0 and cos > 0): #Q4
        return thetaS + 360
    else:
        return 0

def vectorRotation(vec, axis, angle_rad): # active rotation
    vec = np.matrix([[vec[0]], [vec[1]], [vec[2]]]) #convert 1x3 to 3x1 for matrix multiplication
    match axis:
        case Axis.X:
            rotMatrix = np.array([[1, 0, 0],
                                   [0, cos(angle_rad), -sin(angle_rad)],
                                   [0, sin(angle_rad), cos(angle_rad)]])
        case Axis.Y:
            rotMatrix = np.array([[cos(angle_rad), 0, sin(angle_rad)],
                                  [0, 1, 0],
                                  [-sin(angle_rad), 0, cos(angle_rad)]])
        case Axis.Z:
            rotMatrix = np.array([[cos(angle_rad), -sin(angle_rad), 0],
                                  [sin(angle_rad), cos(angle_rad), 0],
                                  [0, 0, 1]])
    res = rotMatrix @ vec
    return np.array([res[0, 0], res[1, 0], res[2, 0]])

def solveUsingNewtonsMethod(xGuess, f, fPrime):
    # catch f'(xGuess) = 0
    if fPrime(xGuess) == 0:
        xGuess += 0.1
    # update xGuess
    xNew = xGuess - f(xGuess) / fPrime(xGuess)
    print(xNew) # debug
    if abs(xNew - xGuess) > 1.0e-10:
        xNew = solveUsingNewtonsMethod(xNew, f, fPrime) # try again if not within bound yet
    return xNew