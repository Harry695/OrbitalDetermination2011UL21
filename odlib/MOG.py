import numpy as np

from odlib.mathUtils import magnitude

def getFAndGConstants(tau1, tau3, posVec2, velVec2):
    posMag2 = magnitude(posVec2)
    u = 1 / posMag2**3.
    z = np.dot(posVec2, velVec2) / posMag2**2.
    q = np.dot(velVec2, velVec2) / posMag2**2. - u

    f1 = 1 - (u * tau1**2) / 2. + (u * z * tau1**3) / 2. + (3 * u * q - 15 * u * z**2 + u**2) * tau1**4 / 24.
    g1 = tau1 - (u * tau1**3) / 6. + (u * z * tau1**4) / 4.
    f3 = 1 - (u * tau3**2) / 2. + (u * z * tau3**3) / 2. + (3 * u * q - 15 * u * z**2 + u**2) * tau3**4 / 24.
    g3 = tau3 - (u * tau3**3) / 6. + (u * z * tau3**4) / 4.
    return f1, f3, g1, g3