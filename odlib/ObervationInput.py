from odlib.constants import Constants


class ObservationInput:
    def __init__(self, time_Jd, ra, dec, sunEarthVector) -> None:
        self.time_Jd = time_Jd
        self.ra = ra
        self.dec = dec
        self.sunEarthVector = sunEarthVector
    
    def getTimeGd(self):
        return self.time_Jd / Constants.DAY_IN_GAUSSIAN_DAY
    
    def __str__(self) -> str:
        return f'''\nObservation at {self.time_Jd} Julian day:
                 RA: {self.ra}
                 DEC: {self.dec}
                 Sun-Earth vector: {self.sunEarthVector}\n'''