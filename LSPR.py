import numpy as np
from odlib.conversions import*

# parse data
file = "LSPRtestinput1.txt"
data = [e.strip() for e in open(file).readlines()]
dataMatrix = []
for obj in data:
    dataMatrix.append(obj.split())
dataMatrix = np.array(dataMatrix).transpose() # easier to do stuff to rows to column

print(data) # debug
print(dataMatrix) # debug

# numberify all values
dataMatrix[0] = [float(e) for e in dataMatrix[0]]# x
dataMatrix[1] = [float(e) for e in dataMatrix[1]]# y
for i in range(len(dataMatrix[2])): # ra
    ra = [float(e) for e in dataMatrix[2, i].split(":")]
    dataMatrix[2, i] = HMStoDeg(ra[0], ra[1], ra[2])
for i in range(len(dataMatrix[3])): # dec
    dec = [float(e) for e in dataMatrix[3, i].split(":")]
    dataMatrix[3, i] = DMStoDeg(dec[0], dec[1], dec[2]) 

print(dataMatrix) #debug

# form matrix variables
n = len(dataMatrix[0]) 
print(n)

