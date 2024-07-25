from math import degrees, radians
from odlib.constants import Constants


class ObservationInput:
    def __init__(self, time_Jd, ra_deg, dec_deg, earthSunVector) -> None:
        self.time_Jd = time_Jd
        self.ra = ra_deg
        self.dec = dec_deg
        self.earthSunVector = earthSunVector
    
    def getTimeGd(self):
        return self.time_Jd / Constants.DAY_IN_GAUSSIAN_DAY
    
    def lightSpeedCorrection(self, rhoMag):
        return self.time_Jd - rhoMag / Constants.LIGHT_SPEED_AU_DAY

    def __str__(self) -> str:
        return f'''\nObservation at {self.time_Jd} Julian day:
                 RA: {self.ra}
                 DEC: {self.dec}
                 Sun-Earth vector: {self.earthSunVector}\n'''