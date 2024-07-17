import numpy as np

class Constants:
    M_IN_AU = 1.496e11
    DAY_IN_GAUSSIAN_DAY = 58.13244087 # Clearly the period is now not in a usual unit like seconds, days, or years. We call this new unit of time the Gaussian day. The Earth's period is 1 year = 2PI Gaussian days. Put another way, 1 Gaussian day is the time it takes for a hypothetical object on a perfectly circular orbit around the Sun with radius of 1 AU to move 1 radian around its orbit. 
    EARTH_MASS_KG = 5.97219e24
    SUN_MASS_KG = 1.989e30
    MU = 1 # in AU and Gaussian days
    EARTH_TILT_DEG = 23.4384987711 

    EARTH_SUN_2458333_HALF = np.array([-6.574011189521245E-01, 7.092445973782825E-01, 3.074588267894852E-01])