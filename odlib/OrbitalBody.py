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
        """
        Returns the current mean anomaly in radians.
        """
        # print("time_gd", time_gD)
        # print("t0", self.t0)
        # print("n", 1 / self.a**1.5)
        return (time_gD - self.t0) / self.a**1.5

    def getEclipticCoords(self, time_gD):
        M = self.getCurrentMeanAnomaly(time_gD)
        # print("M", M)
        E = getEccentricAnomaly(M, self.e)  # works
        orbitalPos = np.array(
            [
                self.a * cos(E) - self.a * self.e,
                self.a * sqrt(1 - self.e**2) * sin(E),
                0,
            ]
        )
        # print("Orbital vec", orbitalPos)
        wRot = vectorRotation(orbitalPos, Axis.Z, self.w)
        iRot = vectorRotation(wRot, Axis.X, self.i)
        omegaRot = vectorRotation(iRot, Axis.Z, self.omega)  # now in ecliptic coords
        return omegaRot

    def getEarthBodyVector(self, time_gD):
        eclCoords = self.getEclipticCoords(time_gD)
        eqCoords = vectorRotation(eclCoords, Axis.X, radians(Constants.EARTH_TILT_DEG))
        # print("eq coords", eqCoords)
        # print("Earth-sun", Constants.EARTH_SUN_2458333_HALF)
        rho = eqCoords + Constants.EARTH_SUN_2458333_HALF  # hard coded date
        return rho

    def getRAandDEC(self, time_gD):
        """
        Returns RA and DEC of the body in degrees.
        """
        rho = self.getEarthBodyVector(time_gD)
        # print("rho", rho)
        rhoMag = magnitude(rho)
        dec = asin(rho[2] / rhoMag)
        # print("DEC", dec)
        sinRA = rho[1] / (rhoMag * cos(dec))
        cosRA = rho[0] / (rhoMag * cos(dec))
        return quadrantCorrection(sinRA, cosRA), degrees(dec)

    def printAllOrbitalElements(self):
        print("-" * 10, "Orbital elements", "-" * 10)
        print("a:", self.a)
        print("e:", self.e)
        print("i", self.i)
        print("omega:", degrees(self.omega))
        print("w:", degrees(self.w))
        print("M", self.M)
        print("-" * 20)
