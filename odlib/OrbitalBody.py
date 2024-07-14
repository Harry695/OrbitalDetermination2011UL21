from odlib.orbitalCharacteristics import *
from odlib.mathUtils import *

class OrbitalBody:
    def __init__(self, posVec, velVec, jD) -> None:
        """
        Units are AU and Gd
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
        """Returns the current mean anomaly in degrees"""
        return degrees((time_gD - self.t0) / self.a**1.5)

    def getEclipticCoords(self, time_gD):
        M = self.getCurrentMeanAnomaly(self, time_gD)
        E = getEccentricAnomaly(M, self.e)
        orbitalPos = np.array([
                        self.a * cos(E) - self.a * self.e,
                        self.a * sqrt(1 - self.e**2) * sin(E),
                        0])
