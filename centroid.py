from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from math import sqrt

class BackgroundMethods(Enum):
    MEAN = "Mean"
    MEDIAN = "Median"
    FIRST_TO_THIRD_QUARTILE_MEAN = "First to third quartile mean"
    MODE = "Mode"

def findCentroid(file, bkgdMethod, targetX, targetY, r=3, inSkyR=5, outSkyR=9):
    if bkgdMethod not in BackgroundMethods:
        raise ValueError(f"The value {bkgdMethod} is not in {[e.value for e in BackgroundMethods]}")

    # open file
    data = fits.getdata(file)
    # debug
    # plt.imshow(data, origin="lower", vmin=2.0e3, vmax=2.8e3)

    # sum up background
    bkgdSum = 0
    bkgdCounter = 0
    for row in range(len(data)):
        for col in range(len(data[0])):
            distance = getDistance(row, col, 0, 0)
            if (distance > inSkyR and distance < outSkyR):
                bkgdCounter += 1
                bkgdSum += data[row, col]
    
    # average all pixel in background
    bkgdAvg = bkgdSum / bkgdCounter
    print("Background average", bkgdAvg)

    # define aperture
    aperture = data[targetX - r : targetX + r + 1, targetY - r : targetY + r + 1]
    # print(targetX, r, targetX - r, targetX + r + 1) # debug
    # print(targetY, r, targetY - r, targetY + r + 1)
    # substract bkgd
    np.subtract(aperture, bkgdAvg)
    # print(np.shape(aperture)) #debug
    # calculate sum of aperture
    apSum = np.sum(aperture)

    # calculate xcm 
    xSumArr = np.sum(aperture, axis=0)
    xCoordsArr = np.arange(targetX - r, targetX + r + 1)
    # print(xCoordsArr) #debug
    sumX = xCoordsArr * xSumArr
    xcm = sumX.sum() / apSum
    # calculate ycm
    YSumArr = np.sum(aperture, axis=1)
    yCoordsArr = np.arange(targetY - r, targetY + r + 1)
    sumY = yCoordsArr * YSumArr
    ycm = sumY.sum() / apSum

    # calculate xsigma
    xSigma = getStDev(xSumArr, xCoordsArr, xcm, apSum)
    # calculate ysigma
    ySigma = getStDev(YSumArr, yCoordsArr, ycm, apSum)
    #end
    # plt.show()
    return xcm, ycm, xSigma, ySigma

def getStDev(valueArr, coordsArr, center, sum):
    return sqrt(np.sum(valueArr * (coordsArr - center)**2) / (sum**2 - sum))

def getDistance(x, y, x0, y0):
    return sqrt((x - x0)**2 + (y - y0)**2)

plt.gray()

#REMEMBER TO SUSBTRACT 1,1 FROM DS9 COORDS
centroid_x, centroid_y, uncert_x, uncert_y = findCentroid(
    "centroid_sample.fits", BackgroundMethods.MEAN, 351, 154, r=4, inSkyR=7, outSkyR=12)

# centroid_x, centroid_y, uncert_x, uncert_y = findCentroid("sampleimage.fits", 459, 397, 2)

if abs(centroid_x - 350.7806) < 0.1 and abs(centroid_y - 153.5709) < 0.1:
    print("centroid calculation CORRECT")
else:
    print(
        "centroid calculation INCORRECT, expected (350.7806, 153.5709), got ({}, {})".format(
            centroid_x, centroid_y
        )
    )
