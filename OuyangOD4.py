from odlib.OrbitalBody import OrbitalBody
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
print("Semimajor Axis:", a, f"with {(a - 1.056800055578855E+00) / 1.056800055578855E+00 :%} error")
e = getEccentricity(posVec, velVec)
print("Eccentricity:", e, f"with {(e - 3.442331103932664E-01) / 3.442331103932664E-01 :%} error")
i = degrees(getInclination(posVec, velVec))
print("Inclination:", i, f"with {(i - 2.515525144198502E+01) / 2.515525144198502E+01 :%} error")
omega = degrees(getLongitudeOfAscendingNode(posVec, velVec))
print("Longitude of ascending node:", omega, f"with {(omega - 2.362379803959657E+02) / 2.362379803959657E+02 :%} error")
w = degrees(getArgumentOfPerihelion(posVec, velVec))
print("Argument of perihelion:", w, f"with {(w - 2.555046145241498E+02) / 2.555046145241498E+02 :%} error")
M = degrees(getMeanAnomaly(posVec, velVec))
print("Mean anomaly:", M, f"with {(M - 1.404194574969259E+02) / 1.404194574969259E+02 :%} error")

time_Jd = gregorianToJulianDay(2018, 7, 14, 0, 0, 0)
print("Julian date:", time_Jd)
lastPerihelionTime_Jd = getLastPerihelionTime(a, M, time_Jd / Constants.DAY_IN_GAUSSIAN_DAY) * Constants.DAY_IN_GAUSSIAN_DAY
print("Last perihelion time in Jd:", lastPerihelionTime_Jd)


asteroid = OrbitalBody(posVec, velVec, 2458313.5)
asteroid.printAllOrbitalElements()
print("Ecliptic coords:", asteroid.getEclipticCoords(2458333.5 / Constants.DAY_IN_GAUSSIAN_DAY))
radec = asteroid.getRAandDEC(2458333.5 / Constants.DAY_IN_GAUSSIAN_DAY)
printHMS(radec[0])
print(f"RA error is {(radec[0] - HMStoDeg(17, 42, 21.20)) / HMStoDeg(17, 42, 21.20) :%}")
printDMS(radec[1])
print(f"DEC error is {(radec[1] - DMStoDeg(31, 52, 28.1)) / DMStoDeg(31, 52, 28.1) :%}")