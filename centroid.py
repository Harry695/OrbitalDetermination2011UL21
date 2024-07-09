from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from math import sqrt

class BackgroundMethods(Enum):
    MEAN = "Mean"
    MEDIAN = "Median"

def findCentroid(file, bkgdMethod, targetX, targetY, r=3, inSkyR=5, outSkyR=9):
    """
    Checks:
        Enforce targetX and targetY to integers
        Constrain background method to valid options
        Ensure valid radii
    """
    # prevent floats
    targetX = int(targetX)
    targetY = int(targetY)
    # make sure background methods is what is designated
    if bkgdMethod not in BackgroundMethods:
        raise ValueError(f"The value {bkgdMethod} is not in {[e.value for e in BackgroundMethods]}")
    # make sure radii are valid
    if inSkyR < r or r < 0 or inSkyR < 0 or inSkyR < 0:
        raise ValueError("The radii are not valid!!!!")

    #flip x and y to correct for swapping of x and y to rows and cols
    temp = targetX
    targetX = targetY
    targetY = temp

    # open file
    data = fits.getdata(file)
    # debug
    # plt.imshow(data, origin="lower", vmin=2.0e3, vmax=2.8e3)

    # sum up background
    data2 = data.copy()
    bkgdCounter = 0
    for row in range(targetX - outSkyR, targetX + outSkyR + 1):
        for col in range(targetY - outSkyR, targetY + outSkyR + 1):
            distance = getDistance(row, col, targetX, targetY)
            if (distance < inSkyR or distance > outSkyR):
                data2[row, col] = 0
            else: 
                bkgdCounter += 1
    bkgd = data2[targetX - outSkyR : targetX + outSkyR + 1, targetY - outSkyR : targetY + outSkyR + 1]
    # average all pixel in background
    bkgdAvg = 0
    if bkgdMethod == BackgroundMethods.MEAN:
        bkgdAvg = int(np.sum(bkgd) / bkgdCounter)
    elif bkgdMethod == BackgroundMethods.MEDIAN:
        bkgdAvg = np.median(bkgd)
        print(np.median(bkgd))

    print("Background average", bkgdAvg)

    # define aperture
    aperture = data[targetX - r : targetX + r + 1, targetY - r : targetY + r + 1]
    # substract bkgd
    aperture = np.clip(aperture, bkgdAvg, 65535)
    print(aperture)
    aperture = np.subtract(aperture, bkgdAvg)
    # circular aperture
    for i in range(len(aperture)):
        for j in range(len(aperture[0])):
            distance = getDistance(i, j, len(aperture)/2.0-0.7, len(aperture[0])/2.0-1)
            if distance > r:
                aperture[i, j] = 0
    # print(targetX, r, targetX - r, targetX + r + 1) # debug
    # print(targetY, r, targetY - r, targetY + r + 1)
    # print(np.shape(aperture)) #debug
    print(aperture)
    # aperture[aperture < 0] = 0  
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
    print(aperture)
    print(xCoordsArr)
    return ycm, xcm, xSigma, ySigma

def getStDev(valueArr, coordsArr, center, sum):
    return sqrt(np.sum(valueArr * (coordsArr - center)**2) / (sum**2 - sum))

def getDistance(x, y, x0, y0):
    return sqrt((x - x0)**2 + (y - y0)**2)

plt.gray()

#REMEMBER TO SUSBTRACT 1,1 FROM DS9 COORDS
centroid_x, centroid_y, uncert_x, uncert_y = findCentroid(
    "centroid_sample.fits", BackgroundMethods.MEDIAN, 351, 154, r=3, inSkyR=6, outSkyR=12)
print("Sigma X", uncert_x)
print("Sigma Y", uncert_y)

# centroid_x, centroid_y, uncert_x, uncert_y = findCentroid("sampleimage.fits", 459, 397, 2)

if abs(centroid_x - 350.7806) < 0.1 and abs(centroid_y - 153.5709) < 0.1:
    print(f"centroid calculation CORRECT: {centroid_x}, {centroid_y}")
else:
    print(
        "centroid calculation INCORRECT, expected (350.7806, 153.5709), got ({}, {})".format(
            centroid_x, centroid_y
        )
    )
