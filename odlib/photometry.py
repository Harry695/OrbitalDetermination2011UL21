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
    
    # average all pixel in background
    bkgdAvg = 0
    if bkgdMethod == BackgroundMethods.MEAN:
        bkgdAvg = np.nanmean(bkgd[bkgd!=0])
    elif bkgdMethod == BackgroundMethods.MEDIAN:
        bkgdAvg = np.nanmedian(bkgd[bkgd!=0])
    # plt.imshow(bkgd)

    # print("Background average", bkgdAvg)

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

    # print(f"Centroid at {xcm :2f}, {ycm: 2f}")

    # recenter aperture
    xcm = int(xcm)
    ycm = int(ycm)
    y, x = np.ogrid[-r : r + 1, -r : r + 1]
    apertureMask = x**2 + y**2 <= r**2
    newAperture = data[ycm - r : ycm + r + 1, xcm - r : xcm + r + 1]
    newAperture[~apertureMask] = 0
    # plt.imshow(newAperture)
    # newAperture = np.subtract(np.clip(newAperture, bkgdAvg, None), bkgdAvg)

    # photometry variables
    readNoise = 11 # e-
    darkCurrent = 10 # e-
    gain = 1
    nPix = np.count_nonzero(newAperture)
    # nPix = 77
    signal = gain * np.sum(newAperture) - bkgdAvg * nPix
    # print("signal", signal)
    # print("nPix", nPix)
    nBkgd = np.count_nonzero(bkgd)
    # print("nBkgd", nBkgd)
    # bkgd = np.clip(bkgd, bkgdAvg, 65535)
    # background = gain * np.sum(np.subtract(bkgd, bkgdAvg))
    # print(np.subtract(bkgd, bkgdAvg))
    # print("background", background)
    rho2 = readNoise ** 2 + gain ** 2 / 12
    # print(rho2)
    ab = 1 + nPix / nBkgd
    # print("ab", ab)
    # print(nPix * ab * (bkgdAvg + darkCurrent))
    snr = signal / sqrt(signal + nPix * ab * (bkgdAvg + darkCurrent) + nPix * ab * rho2)
    # print("SNR", snr)
    # print("Signal uncertainty", signal / snr)
    mInst = -2.5 * log10(signal)
    # print("mInst", mInst)
    uncertInst = 1.0875 / snr
    # print("uncertInst", uncertInst)

    return mInst, uncertInst

def differentialPhotometry(starFile, image, targetX, targetY):
    # parse data
    data = [e.strip() for e in open(starFile).readlines()]
    dataMat = []
    for obj in data:
        dataMat.append(obj.split())
    dataMat = np.array(dataMat)
    dataMat = dataMat.astype(float)
    # print(dataMat) # debug

    # calculate mInst
    mInstArr = []
    for star in dataMat:
        mInst = aperturePhotometry(image, BackgroundMethods.MEDIAN, star[0], star[1], 5, 8, 13)[0]
        mInstArr.append(mInst)
    # print(mInstArr) # debug

    # calculate c
    transposedMat = np.transpose(dataMat)
    print(transposedMat[2])
    cList = np.subtract(transposedMat[2], np.array(mInstArr))
    # print(cList) # debug
    cSum = np.sum(cList)
    # print(cSum) # debug
    c = cSum / len(dataMat)
    print("C value:", c) # debug

    # calculate dm
    dm = np.std(cList)
    print("dm value:", dm) # debug

    # calculate catalog magnitude of target star
    output = aperturePhotometry(image, BackgroundMethods.MEDIAN, targetX, targetY, 5, 8, 13)
    mCat = output[0] + c
    print("Catalogue magnitude of requested star:", mCat) # debug
    uncert = output[1]
    print("uncert", uncert) # debug

    # error calculation
    error = sqrt(dm ** 2 + uncert ** 2)
    print("Catalogue magnitude error: ", error)

    return mCat, error

def differentialPhotometryLOBF(signal, vMag, asteroid_signal):
    """
    Inputs are arrays of signal and vMag of the surrounding stars.
    """
    x = np.log10(signal)
    y = vMag
    #find line of best fit
    m, b = np.polyfit(x, y, 1)
    #add points to plot
    plt.scatter(x, y)
    #add line of best fit to plot
    plt.plot(x, m*x+b)
    #compute asteroid's V magnitude
    Vmag_asteroid = m*np.log10(asteroid_signal) + b
    return Vmag_asteroid

# aperturePhotometry("aptest.fit", BackgroundMethods.MEDIAN, 490, 293, r=5, inSkyR=8, outSkyR=13)
differentialPhotometry("stars_test.txt", "diff_phot.fit", 546, 327)
plt.gray()
# plt.show()
differentialPhotometryLOBF()