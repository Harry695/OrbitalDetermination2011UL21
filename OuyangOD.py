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
    earthSunVec = np.asfarray(data[3:6])
    # print("input earthSunVec", earthSunVec)
    obs.append(ObservationInput(float(data[0]), HMSArrToDeg(ra_decimal), DMSArrToDeg(dec_decimal), earthSunVec))
# for e in obs:
#     print(e)

UL21 = OrbitalBody.fromObservations(obs)
UL21.printAllOrbitalElements()

correct = OrbitalBody(np.array([0.0557307, -1.12581648, 0.0550659177]),
                      np.array([0.019095149, 0.00559168, 0.0003274835]) * Constants.DAY_IN_GAUSSIAN_DAY, 
                      2460125.708333333)
print("correct r2", np.array([0.0557307, -1.12581648, 0.0550659177]))
print("correct r2dot", np.array([0.019095149, 0.00559168, 0.0003274835]) * Constants.DAY_IN_GAUSSIAN_DAY)
print("correct orbital elements")
correct.printAllOrbitalElements()

"""
true orbital elements of test case
a = 2.30430
e = 0.54791
i = 3.2412deg
omega = 213.218deg
w = 98.1439deg
"""

"""
true r2 and r2dot of test case
r2 [0.0557307, -1.12581648, 0.0550659177]
r2dot [0.019095149, 0.00559168, 0.0003274835]
"""

"""
true orbital elements of 2011UL21
a = 2.123041010121101
e = 0.6536184976531976
i = 34.88591511692778
omega = 275.3169229735104
w = 285.0014058281747
"""