from odlib.constants import*
from odlib.conversions import*
from odlib.mathUtils import*
from odlib.orbitalCharacteristics import*
import numpy as np
import math

file = open("OuyangInput.txt").readlines()
qf15 = []
for i in range(len(file)-2): # len(file)-2 to get rid of annoying stuff at the end. See the file.
    item = file[i].strip()
    qf15.append(float(item))

# test asteroid 2002 QF15
posVec = np.array([qf15[0], qf15[1], qf15[2]])
velVec = np.array([qf15[3], qf15[4], qf15[5]]) * Constants.DAY_IN_GAUSSIAN_DAY
print("Specific Angular momentum:", getAngularMomentum(posVec, velVec))
print("Semimajor Axis:", getSemimajorAxis(posVec, velVec))
print("Eccentricity:", getEccentricity(posVec, velVec))
print("Inclination:", degrees(getInclination(posVec, velVec)))
print("Longitude of ascending node:", getLongitudeOfAscendingNode(posVec, velVec))
print("Argument of perihelion:", getArgumentOfPerihelion(posVec, velVec))

#DEBUG
# print("HMStoDeg test", HMStoDeg(12, 3, 5.3)) #expected 180.7720833
# print("DMStoDeg test", DMStoDeg(-13, 45, 23.45)) #expected -13.75651389
# deg = 0
# print("Quadrant correction for", deg, quadrantCorrection(sin(math.radians(deg)), cos(math.radians(deg))))
# print("RAdecimalToHMS test 180.7720833:")
# printHMS(180.7720833)
# print("DecDecimalToDMS test 180.7720833:")
# printDMS(-13.75651389)
# print("VectorRotation test:\n", vectorRotation(np.array([-0.23, 0.887, 0.34]), Axis.Y, 45))