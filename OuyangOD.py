from odlib.OrbitalBody import OrbitalBody
from odlib.ObervationInput import ObservationInput
from odlib.constants import*
from odlib.conversions import*
from odlib.mathUtils import*
from odlib.orbitalCharacteristics import*
from odlib.MOG import*
import numpy as np

file = open("OuyangInput.txt").readlines() # [time_Jd, ra_hr, dec_hr, R_x, R_y, R_z]
obs = []
for line in file:
    data = line.strip().split()
    ra_decimal = np.asfarray(data[1].split(":"))
    # print(ra_decimal)
    dec_decimal = np.asfarray(data[2].split(":"))
    # print(dec_decimal)
    sunEarthVec = np.asfarray(data[3:6])
    obs.append(ObservationInput(float(data[0]), HMSArrToDeg(ra_decimal), DMSArrToDeg(dec_decimal), sunEarthVec))

for e in obs:
    print(e)