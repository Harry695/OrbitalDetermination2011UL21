import numpy as np

file = "LSPRtestinput1.txt"
data = open(file).readlines()
data = [e.strip() for e in data]
dataMatrix = []
for obj in data:
    dataMatrix.append(obj.split())
dataMatrix = np.array(dataMatrix)

print(data) # debug
print(dataMatrix) #debug


