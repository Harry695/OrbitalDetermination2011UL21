from odlib.constants import*
from odlib.conversions import*
from odlib.mathUtils import*
from odlib.orbitalCharacteristics import*
import numpy as np

file = open("OuyangInput.txt").readlines()
qf15 = []
for i in range(len(file)-2): # len(file)-2 to get rid of annoying stuff at the end. See the file.
    item = file[i].strip()
    qf15.append(float(item))

# test asteroid 2002 QF15
posVec = np.array([qf15[0], qf15[1], qf15[2]])
velVec = np.array([qf15[3], qf15[4], qf15[5]]) * Constants.DAY_IN_GAUSSIAN_DAY

h = getAngularMomentum(posVec, velVec)
print("Specific Angular momentum:", h)
a = getSemimajorAxis(posVec, velVec)
print("Semimajor Axis:", a, f"with {(a - 1.05671892483881) / 1.05671892483881 :%} error")
e = getEccentricity(posVec, velVec)
print("Eccentricity:", e, f"with {(e - .3442798363212599) / .3442798363212599 :%} error")
i = degrees(getInclination(posVec, velVec))
print("Inclination:", i, f"with {(i - 25.15375982051783) / 25.15375982051783 :%} error")
omega = getLongitudeOfAscendingNode(posVec, velVec)
print("Longitude of ascending node:", omega, f"with {(omega - 236.2502480179119) / 236.2502480179119 :%} error")
w = getArgumentOfPerihelion(posVec, velVec)
print("Argument of perihelion:", w, f"with {(w - 255.5316184725226) / 255.5316184725226 :%} error")
M = getMeanAnomaly(posVec, velVec) 
print("Mean anomaly:", M, f"with {(M - 1.397389542591620E+02) / 1.397389542591620E+02 :%} error")

time_Jd = gregorianToJulianDay(2018, 7, 13, 6, 0, 0)
print("Julian date:", time_Jd)
lastPerihelionTime_Jd = getLastPerihelionTime(a, M, time_Jd / Constants.DAY_IN_GAUSSIAN_DAY) * Constants.DAY_IN_GAUSSIAN_DAY
print("Last perihelion time in Jd:", lastPerihelionTime_Jd)

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