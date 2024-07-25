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
        # print("class r2", posVec)
        # print("class r2dot", velVec)
        self.t_gD = jD / Constants.DAY_IN_GAUSSIAN_DAY
        self.a = getSemimajorAxis(posVec, velVec)
        # print("class a", self.a)
        # print("h:", getAngularMomentum(posVec, velVec))
        self.e = getEccentricity(posVec, velVec)
        self.i = getInclination(posVec, velVec)
        self.omega = getLongitudeOfAscendingNode(posVec, velVec)
        self.w = getArgumentOfPerihelion(posVec, velVec)
        self.M = getMeanAnomaly(posVec, velVec)
        self.t0 = getLastPerihelionTime(
            self.a, self.M, jD / Constants.DAY_IN_GAUSSIAN_DAY
        )

    @classmethod
    def fromObservations(cls, obs): # TODO: implement more than 3 obs
        """
        obs: a list of 3 ObservationInput objects.
        """
        if not isinstance(obs[0], ObservationInput):
            raise ValueError("ERROR: please input ObservationInput objects!")
        
        MOG_THRESHOLD = 1.0e-12
        ITERATION_CAP = 500

        # unpack obervations
        obs1, obs2, obs3 = obs

        # calculate tau's
        tau1, tau3, tau0 = getTauArr(obs1.time_Jd, obs2.time_Jd, obs3.time_Jd) # np array unpacking - possibly sketchy
        # print("taus", tau1, tau3, tau0)

        # build rhoHats
        P1DIR = getRhoDirection(radians(obs1.ra), radians(obs1.dec))
        P2DIR = getRhoDirection(radians(obs2.ra), radians(obs2.dec))
        P3DIR = getRhoDirection(radians(obs3.ra), radians(obs3.dec))
        # print("pDir", P1DIR, P2DIR, P3DIR, sep="\n")

        # compute D constants
        # DARR = [getScalarEquationDConstants(P1DIR, P2DIR, P3DIR, obs1.earthSunVector), 
        #         getScalarEquationDConstants(P1DIR, P2DIR, P3DIR, obs2.earthSunVector),
        #         getScalarEquationDConstants(P1DIR, P2DIR, P3DIR, obs3.earthSunVector)]
        # DARR = np.transpose(DARR) # !IMPORTANT! transpose DArr for input

        # alt
        D11, D21, D31 = getScalarEquationDConstants(P1DIR, P2DIR, P3DIR, obs1.earthSunVector)
        D12, D22, D32 = getScalarEquationDConstants(P1DIR, P2DIR, P3DIR, obs2.earthSunVector)
        D13, D23, D33 = getScalarEquationDConstants(P1DIR, P2DIR, P3DIR, obs3.earthSunVector)
        DARR = [[D11, D12, D13],
                [D21, D22, D23],
                [D31, D32, D33]]
        # print("dArr", DARR)
        D0 = getD0(P1DIR, P2DIR, P3DIR)
        # print("D0", D0)
        
        # initial guesses
        c1Guess = tau3 / tau0
        c3Guess = -tau1 / tau0
        # print("initial c constants", c1Guess, c3Guess)
        
        # calculate distance
        p1Mag, p2Mag, p3Mag = getDistances([c1Guess, -1, c3Guess], D0, DARR) 
        # print("inital pmag", p1Mag, p2Mag, p3Mag)

        # find p's
        p1Vec = p1Mag * P1DIR
        p2Vec = p2Mag * P2DIR
        p3Vec = p3Mag * P3DIR

        # find rVec's
        posVec1 = p1Vec - obs1.earthSunVector
        posVec2 = p2Vec - obs2.earthSunVector
        posVec3 = p3Vec - obs3.earthSunVector
        # print("initial r1", posVec1)
        # print("inital r3", posVec3)

        # guess r2dot
        r12Dot = -(posVec2 - posVec1) / tau1
        r23Dot = (posVec3 - posVec2) / tau3
        # print("r12", r12Dot)
        # print("r23", r23Dot)
        velVec2 = (tau3 / tau0) * r12Dot - (tau1 / tau0) * r23Dot

        print("inital r2 guess", posVec2)
        print("inital r2dot guess", velVec2)

        count = 0
        while True:
            count += 1
            # print(count)
            # time correction
            newt1_Jd = obs1.lightSpeedCorrection(p1Mag)
            newt2_Jd = obs2.lightSpeedCorrection(p2Mag)
            newt3_Jd = obs3.lightSpeedCorrection(p3Mag)
            # print("light speed corrected t", newt1_Jd, newt2_Jd, newt3_Jd)
            tau1, tau3, tau0 = getTauArr(newt1_Jd, newt2_Jd, newt3_Jd) # np array unpacking - possibly sketchy

            # ready for iterative process
            f1, f3, g1, g3 = getFAndGConstants(tau1, tau3, posVec2, velVec2)
            c1, c3 = getcConstants(f1, f3, g1, g3)
            d1, d3 = getdConstants(f1, f3, g1, g3)

            # calculate distance
            p1Mag, p2Mag, p3Mag = getDistances([c1, -1, c3], D0, DARR) # names sketchy
            # print("rhomags", p1Mag, p2Mag, p3Mag)

            # find p's
            p1Vec = P1DIR * p1Mag
            p2Vec = P2DIR * p2Mag
            p3Vec = P3DIR * p3Mag

            # find rVec's
            newPosVec1 = p1Vec - obs1.earthSunVector
            newPosVec2 = p2Vec - obs2.earthSunVector
            newPosVec3 = p3Vec - obs3.earthSunVector
            
            # find r2Dot
            newVelVec2 = d1 * newPosVec1 + d3 * newPosVec3

            # if count <= 5:
                # print("iterative r2", posVec2)
                # print("iterative r2dot", velVec2)

            # determine when to end loop
            if (abs(newPosVec2[0] - posVec2[0]) < MOG_THRESHOLD and 
                abs(newPosVec2[1] - posVec2[1]) < MOG_THRESHOLD and
                abs(newPosVec2[2] - posVec2[2]) < MOG_THRESHOLD):
                posVec2 = newPosVec2
                velVec2 = newVelVec2
                print("final r2", posVec2)
                break
            else: # if not within threshold, continue
                posVec2 = newPosVec2
                velVec2 = newVelVec2
            # if (abs(p2MagNew - p2Mag) < MOG_THRESHOLD):
            #     break
            # else:
            #     p1Mag = p1MagNew
            #     p2Mag = p2MagNew
            #     p3Mag = p3MagNew
            # print("")

            # nonconvergence check
            if (count > ITERATION_CAP):
                print("NONCONVERGING MOG")
                return None

        # convert from equatorial to ecliptic
        posVec2 = vectorRotation(posVec2, Axis.X, -radians(Constants.EARTH_TILT_DEG)) 
        velVec2 = vectorRotation(velVec2, Axis.X, -radians(Constants.EARTH_TILT_DEG))
        print("\necliptic r2", posVec2)
        print("ecliptic r2dot", velVec2)
        print("a", getSemimajorAxis(posVec2, velVec2))

        # return
        # print("start class construction with", posVec2, velVec2, newt2_Jd)
        return cls(posVec2, velVec2, newt2_Jd)

            
    def getCurrentMeanAnomaly(self, time_gD):
        """
        Returns the current mean anomaly in radians.
        """
        # print("time_gd", time_gD)
        # print("t_gD", self.t_gD)
        # print("M", self.M)
        # print("change in M", (time_gD - self.t_gD) / self.a**1.5)
        return (self.M + (time_gD - self.t_gD) / self.a**1.5) % (2 * np.pi)

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
        rho = eqCoords + Constants.EARTH_SUN_2458333_HALF  # TODO: replace hard-coded earthsunvec
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
        print("i", degrees(self.i))
        print("omega:", degrees(self.omega))
        print("w:", degrees(self.w))
        print("M", degrees(self.M))
        print("-" * 40)
    
    def getAllOrbitalElements(self):
        """
        Returns an array of [a, e, i, omega, w] in AU, Gd, and degrees.
        """
        return np.array([self.a, self.e, degrees(self.i), degrees(self.omega), degrees(self.w)])
