import numpy as np
from math import sqrt
from odlib.conversions import *

file = "LSPRtestinput1.txt"


def parseData(file):
    # load data
    data = [e.strip() for e in open(file).readlines()]
    dataMat = []
    for obj in data:
        dataMat.append(obj.split())
    dataMat = np.array(dataMat).transpose()  # easier to do stuff to rows to column

    # numberify all values
    for i in range(len(dataMat[2])):  # ra
        ra = [float(e) for e in dataMat[2, i].split(":")]
        dataMat[2, i] = HMStoDeg(ra[0], ra[1], ra[2])  # convert to decimal deg
    for i in range(len(dataMat[3])):  # dec
        dec = [float(e) for e in dataMat[3, i].split(":")]
        dataMat[3, i] = DMStoDeg(dec[0], dec[1], dec[2])  # convert to decimal deg
    dataMat = dataMat.astype(float)  # make the whole array floats

    return dataMat


def findPlateConstants(dataMat):
    # print(dataMat) # debug

    # form matrix variables
    n = len(dataMat[0])
    # print("n", n) #debug
    sumX, sumY, sumRA, sumDEC = np.sum(dataMat, axis=1)
    # print("sumX", sumX) #checked
    # print("sumY", sumY) #checked
    # print("sumRA", sumRA)
    # print("sumDEC", sumDEC)
    sumX2 = np.sum(dataMat[0] ** 2)
    # print("sumX2", sumX2) #checked
    sumY2 = np.sum(dataMat[1] ** 2)
    # print("sumY2", sumY2) #checked
    sumXY = np.sum(dataMat[0] * dataMat[1])
    # print("sumXY", sumXY)
    sumRAX = np.sum(dataMat[0] * dataMat[2])
    # print("sumRAX", sumRAX)
    sumRAY = np.sum(dataMat[1] * dataMat[2])
    # print("SumRAY", sumRAY)
    sumDECX = np.sum(dataMat[0] * dataMat[3])
    # print("SumDECX", sumDECX)
    sumDECY = np.sum(dataMat[1] * dataMat[3])
    # print("SumDECY", sumDECY)

    # building matrices
    matXY = np.linalg.inv(
        np.matrix([[n, sumX, sumY], [sumX, sumX2, sumXY], [sumY, sumXY, sumY2]])
    )
    # print("matXY", matXY) #debug
    matRA = np.matrix([[sumRA], [sumRAX], [sumRAY]])
    # print("matRA", matRA) #debug
    matDEC = np.matrix([[sumDEC], [sumDECX], [sumDECY]])
    # print("matDEC", matDEC) #debug

    # multiply matrices
    mat1 = matXY @ matRA
    mat2 = matXY @ matDEC
    # debug
    # for e in mat1:
    # print(e)
    # for e in mat2:
    # print(e)

    # assigning variables TODO: make less horrible
    b1, a11, a12 = np.array(mat1)
    b2, a21, a22 = np.array(mat2)
    b1 = b1[0]
    b2 = b2[0]
    a11 = a11[0]
    a12 = a12[0]
    a21 = a21[0]
    a22 = a22[0]

    return b1, b2, a11, a12, a21, a22


def findUncert(dataMat, b1, b2, a11, a12, a21, a22):
    dataMat = dataMat.astype(float) #security check
    
    # print(dataMat)
    n = len(dataMat[0])

    # RA
    raSum = 0
    for i in range(len(dataMat[0])):
        raSum += (dataMat[2][i] - b1 - a11 * dataMat[0][i] - a12 * dataMat[1][i]) ** 2
    sigmaRA = sqrt(raSum / (n - 3)) * 3600

    # dec
    decSum = 0
    for i in range(len(dataMat[0])):
        decSum += (dataMat[3][i] - b2 - a21 * dataMat[0][i] - a22 * dataMat[1][i]) ** 2
    sigmaDEC = sqrt(decSum / (n - 3)) * 3600

    return sigmaRA, sigmaDEC


dataMat = parseData(file)
b1, b2, a11, a12, a21, a22 = findPlateConstants(dataMat)
print("Plate constants are", b1, b2, a11, a12, a21, a22, sep="\n")
print("Uncertainties are", findUncert(dataMat, b1, b2, a11, a12, a21, a22))
