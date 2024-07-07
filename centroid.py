from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

class BackgroundMethods(Enum):
    MEAN = "Mean"
    MEDIAN = "Median"
    FIRST_TO_THIRD_QUARTILE_MEAN = "First to third quartile mean"
    MODE = "Mode"

def findCentroid(file, bkgdMethod, targetX, targetY, r=3, inSkyR=5, outSkyR=9):
    if bkgdMethod not in BackgroundMethods:
        raise ValueError(f"The value is not in {BackgroundMethods}")
    #open file
    #define background
    #average all pixel in background
    #substract bkgd
    #calculate xcm
    #calculate ycm
    #calculate xsigma
    #calculate ysigma
    # return xcm, ycm, xsigma, ysigma

def getStDev(arr):
    pass

centroid_x, centroid_y, uncert_x, uncert_y = findCentroid("sampleimage.fits", BackgroundMethods.MEAN, 351, 154, 3, 5, 9)

# centroid_x, centroid_y, uncert_x, uncert_y = findCentroid("sampleimage.fits", 459, 397, 2)

if abs(centroid_x - 350.7806) < 0.1 and abs(centroid_y - 153.5709) < 0.1:
    print("centroid calculation CORRECT")
else:
    print(
        "centroid calculation INCORRECT, expected (350.7806, 153.5709), got ({}, {})".format(centroid_x, centroid_y))
    

