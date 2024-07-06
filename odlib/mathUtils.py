from math import*

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
    