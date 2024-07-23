from math import radians
from odlib.constants import Constants


class ObservationInput:
    def __init__(self, time_Jd, ra_deg, dec_deg, earthSunVector) -> None:
        self.time_Jd = time_Jd
        self.ra = radians(ra_deg)
        self.dec = radians(dec_deg)
        self.earthSunVector = earthSunVector
    
    def getTimeGd(self):
        return self.time_Jd / Constants.DAY_IN_GAUSSIAN_DAY
    
    def lightSpeedCorrection(self, rhoMag):
        return self.time_Jd - rhoMag / Constants.LIGHT_SPEED_M_PER_S

    def __str__(self) -> str:
        return f'''\nObservation at {self.time_Jd} Julian day:
                 RA: {self.ra}
                 DEC: {self.dec}
                 Sun-Earth vector: {self.earthSunVector}\n'''