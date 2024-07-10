import numpy as np
from odlib.conversions import*

# parse data
file = "LSPRtestinput1.txt"
data = [e.strip() for e in open(file).readlines()]
dataMat = []
for obj in data:
    dataMat.append(obj.split())
dataMat = np.array(dataMat).transpose() # easier to do stuff to rows to column

print(data) # debug
print(dataMat) # debug

# numberify all values
for i in range(len(dataMat[2])): # ra
    ra = [float(e) for e in dataMat[2, i].split(":")]
    dataMat[2, i] = HMStoDeg(ra[0], ra[1], ra[2])
for i in range(len(dataMat[3])): # dec
    dec = [float(e) for e in dataMat[3, i].split(":")]
    dataMat[3, i] = DMStoDeg(dec[0], dec[1], dec[2]) 
dataMat = dataMat.astype(float)
print(dataMat) #debug

# form matrix variables
n = len(dataMat[0]) 
print("n", n) #debug
sumX, sumY, sumRA, sumDEC = np.sum(dataMat, axis=1)
print("sumX", sumX) #checked
print("sumY", sumY) #checked
print("sumRA", sumRA)
print("sumDEC", sumDEC)
sumX2 = np.sum(dataMat[0]**2)
print("sumX2", sumX2) #checked
sumY2 = np.sum(dataMat[1]**2)
print("sumY2", sumY2) #checked
sumXY = np.sum(dataMat[0] * dataMat[1])
print("sumXY", sumXY)
sumRAX = np.sum(dataMat[0] * dataMat[2])
print("sumRAX", sumRAX)
sumRAY = np.sum(dataMat[1] * dataMat[2])
print("SumRAY", sumRAY)
sumDECX = np.sum(dataMat[0] * dataMat[3])
print("SumDECX", sumDECX)
sumDECY = np.sum(dataMat[1] * dataMat[3])
print("SumDECY", sumDECY)

# building matrices
matXY = np.linalg.inv(np.matrix([[n, sumX, sumY],
                   [sumX, sumX2, sumXY],
                   [sumY, sumXY, sumY2]]))
print("matXY", matXY)
matRA = np.matrix([[sumRA],
                   [sumRAX],
                   [sumRAY]])
print("matRA", matRA)
matDEC = np.matrix([[sumDEC],
                    [sumDECX],
                    [sumDECY]])
print("matDEC", matDEC)

# multiply matrices
mat1 = matXY @ matRA
mat2 = matXY @ matDEC
# debug
for e in mat1:
    print(e)
for e in mat2:
    print(e)

#assigning variables
b1, a11, a12 = mat1
b2, a21, a22 = mat2