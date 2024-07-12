import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from math import sqrt, log10

from centroid import BackgroundMethods

def aperturePhotometry(file, bkgdMethod, targetX, targetY, r=3, inSkyR=6, outSkyR=12):
    # prevent floats
    targetX = int(targetX)
    targetY = int(targetY)
    
    # make sure background methods is what is designated
    if bkgdMethod not in BackgroundMethods:
        raise ValueError(
            f"The value {bkgdMethod} is not in {[e.value for e in BackgroundMethods]}"
        )
    
    # make sure radii are valid
    if inSkyR < r or r < 0 or inSkyR < 0 or inSkyR < 0:
        raise ValueError("The radii are not valid!!!!")
    
    # open file
    data = fits.getdata(file)
    plt.imshow(data, origin="lower", vmin=1.1e3, vmax=2.8e3)

    # define background
    dataCopy = data.copy()
    bkgd = dataCopy[  # slice background out
        targetY - outSkyR : targetY + outSkyR + 1,
        targetX - outSkyR : targetX + outSkyR + 1,
    ]

    # create mask
    c, v = np.ogrid[-outSkyR : outSkyR + 1, -outSkyR : outSkyR + 1]
    bkgdMask1 = c**2 + v**2 <= outSkyR**2
    bkgdMask2 = c**2 + v**2 >= inSkyR**2
    
    # apply mask
    bkgd[~bkgdMask1] = 0
    bkgd[~bkgdMask2] = 0
    # plt.imshow(bkgd)
    
    # average all pixel in background
    bkgdAvg = 0
    if bkgdMethod == BackgroundMethods.MEAN:
        bkgdAvg = int(np.sum(bkgd) / np.count_nonzero(bkgd))
    elif bkgdMethod == BackgroundMethods.MEDIAN:
        bkgdAvg = np.median(bkgd)

    print("Background average", bkgdAvg)
    # testing
    bkgdAvg = 1243

    # define aperture
    aperture = data.copy()[targetY - r : targetY + r + 1, targetX - r : targetX + r + 1]

    # substract bkgd
    aperture = np.clip(aperture, bkgdAvg, 65535)
    aperture = np.subtract(aperture, bkgdAvg)
    
    # circular mask on aperture
    y, x = np.ogrid[-r : r + 1, -r : r + 1]
    apertureMask = x**2 + y**2 <= r**2
    aperture[~apertureMask] = 0

    # calculate sum of aperture
    # plt.imshow(aperture, origin="lower")
    signal = np.sum(aperture)
    
    # calculate xcm
    xSumArr = np.sum(aperture, axis=1)
    xCoordsArr = np.arange(targetX - r, targetX + r + 1)
    sumX = xCoordsArr * xSumArr
    xcm = sumX.sum() / signal
    
    # calculate ycm
    YSumArr = np.sum(aperture, axis=0)
    yCoordsArr = np.arange(targetY - r, targetY + r + 1)
    sumY = yCoordsArr * YSumArr
    ycm = sumY.sum() / signal

    print(xcm, ycm)
    # recenter aperture
    xcm = int(xcm)
    ycm = int(ycm)
    y, x = np.ogrid[-r : r + 1, -r : r + 1]
    apertureMask = x**2 + y**2 <= r**2
    newAperture = data[ycm - r : ycm + r + 1, xcm - r : xcm + r + 1]
    newAperture[~apertureMask] = 0
    plt.imshow(newAperture)
    # newAperture = np.subtract(np.clip(newAperture, bkgdAvg, None), bkgdAvg)

    # photometry variables
    readNoise = 11 # e-
    darkCurrent = 10 # e-
    gain = 1
    nPix = np.count_nonzero(newAperture)
    # nPix = 77
    signal = gain * np.sum(newAperture) - bkgdAvg * nPix
    print("signal", signal)
    print("nPix", nPix)
    nBkgd = np.count_nonzero(bkgd)
    print("nBkgd", nBkgd)
    bkgd = np.clip(bkgd, bkgdAvg, 65535)
    # background = gain * np.sum(np.subtract(bkgd, bkgdAvg))
    # print(np.subtract(bkgd, bkgdAvg))
    # print("background", background)
    rho2 = readNoise ** 2 + gain ** 2 / 12
    # print(rho2)
    ab = 1 + nPix / nBkgd
    print("ab", ab)
    # print(nPix * ab * (bkgdAvg + darkCurrent))
    snr = signal / sqrt(signal + nPix * ab * (bkgdAvg + darkCurrent) + nPix * ab * rho2)
    # snr = sqrt(signal) / sqrt(1 + nPix * ab * ((bkgdAvg + darkCurrent + rho2) / signal))
    print("SNR", snr)
    mInst = -2.5 * log10(signal)
    print("mInst", mInst)
    uncertInst = 1.0875 / snr
    print("uncertInst", uncertInst)

aperturePhotometry("aptest.fit", BackgroundMethods.MEAN, 490, 293, r=5, inSkyR=8, outSkyR=13)
plt.gray()
plt.show()