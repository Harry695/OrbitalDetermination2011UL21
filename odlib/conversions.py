def HMStoDeg(h, m, s):
    m = s / 60 + m
    h = m / 60 + h
    return h / 24 * 360

def DMStoDeg(d, m, s):
    m = s / 60 + m
    m = m*-1 if d < 0 else m
    return m / 60 + d