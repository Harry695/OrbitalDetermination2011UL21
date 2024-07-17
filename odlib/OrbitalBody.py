from odlib.orbitalCharacteristics import *
from odlib.mathUtils import *

class OrbitalBody:
    def __init__(self, posVec, velVec, jD) -> None:
        """
        Units are AU, Gd, and radians. 
        jD is the time of obervation in Julian days. These arguments are enough to completely determine an asteroid's orbit.
        """
        self.a = getSemimajorAxis(posVec, velVec)
        self.e = getEccentricity(posVec, velVec)
        self.i = getInclination(posVec, velVec)
        self.omega = getLongitudeOfAscendingNode(posVec, velVec)
        self.w = getArgumentOfPerihelion(posVec, velVec)
        self.M = getMeanAnomaly(posVec, velVec)
        self.t0 = getLastPerihelionTime(
            self.a, self.M, jD / Constants.DAY_IN_GAUSSIAN_DAY
        )

    def getCurrentMeanAnomaly(self, time_gD):
        """Returns the current mean anomaly in radians"""
        return (time_gD - self.t0) / self.a**1.5

    def getEclipticCoords(self, time_gD):
        M = self.getCurrentMeanAnomaly(time_gD)
        print("M", M)
        E = getEccentricAnomaly(M, self.e) # works
        orbitalPos = np.array([
                        self.a * cos(E) - self.a * self.e,
                        self.a * sqrt(1 - self.e**2) * sin(E),
                        0])
        wRot = vectorRotation(orbitalPos, Axis.Z, self.w)
        iRot = vectorRotation(wRot, Axis.X, self.i)
        omegaRot = vectorRotation(iRot, Axis.Z, self.omega) # now in ecliptic coords
        return omegaRot

    def printAllOrbitalElements(self):
        print("-" * 10, "Orbital elements", "-" * 10)
        print("a:", self.a)
        print("e:", self.e)
        print("i", degrees(self.i))
        print("omega:", degrees(self.omega))
        print("w:", degrees(self.w))
        print("-" * 20)

