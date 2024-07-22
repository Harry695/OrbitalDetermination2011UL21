from odlib.orbitalCharacteristics import *
from odlib.mathUtils import *
from odlib.ObervationInput import ObservationInput
from odlib.MOG import *

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

    @classmethod
    def fromObservations(cls, obs):
        """
        obs: a list of 3 ObservationInput objects.
        """
        if not isinstance(obs[0], ObservationInput):
            raise ValueError("ERROR: please input ObservationInput objects!")
        
        # unpack obervations
        obs1, obs2, obs3 = obs

        # calculate tau's
        TAU1, TAU3, TAU0 = getTauArr(obs1.time_Jd, obs2.time_Jd, obs3.time_Jd) # np array unpacking - possibly sketchy
        
        # build rhoHats
        P1DIR = getRhoDirection(obs1.ra, obs1.dec)
        P2DIR = getRhoDirection(obs2.ra, obs2.dec)
        P3DIR = getRhoDirection(obs3.ra, obs3.dec)

        # initial guesses
        c1Guess = TAU3 / TAU0
        c3Guess = -TAU1 / TAU0
        
        # get Ds and distance
        DARR = [getScalarEquationDConstants(P1DIR, P2DIR, P3DIR, obs1.earthSunVector), 
                getScalarEquationDConstants(P1DIR, P2DIR, P3DIR, obs2.earthSunVector),
                getScalarEquationDConstants(P1DIR, P2DIR, P3DIR, obs3.earthSunVector)]
        DARR = np.transpose(DARR) # !IMPORTANT! transpose DArr for input
        D0 = getD0(P1DIR, P2DIR, P3DIR)
        p1Mag, p2Mag, p3Mag = getDistances([c1Guess, -1, c3Guess], D0, DARR) 

        # find p's
        p1Vec = p1Mag * P1DIR
        p2Vec = p2Mag * P2DIR
        p3Vec = p3Mag * P3DIR

        # find rVec's
        posVec1 = p1Vec - obs1.earthSunVector
        posVec2 = p2Vec - obs2.earthSunVector
        posVec3 = p3Vec - obs3.earthSunVector

        # guess r2dot
        r12Dot = -(posVec2 - posVec1) / TAU1
        r23Dot = (posVec3 - posVec2) / TAU3
        velVec2 = (TAU3 / TAU0) * r12Dot - (TAU1 / TAU0) * r23Dot

        # ready for iterative process
        f1, f3, g1, g3 = getFAndGConstants(TAU1, TAU3, posVec2, velVec2)
        c1, c3 = getcConstants(f1, f3, g1, g3)
        p1Mag, p2Mag, p3Mag = getDistances([c1, -1, c3], D0, DARR) # names sketchy
        # find rVec's
        # find r2Dot
        # restart

            
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
