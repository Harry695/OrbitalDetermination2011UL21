from math import*

def HMStoDeg(h, m, s):
    m = s / 60 + m
    h = m / 60 + h
    return h / 24 * 360

def DMStoDeg(d, m, s):
    m = s / 60 + m
    m = m*-1 if d < 0 else m # convert m to sign of d for addition
    return m / 60 + d

def printHMS(ra):
    hoursDecimal = ra / 360 * 24
    h = floor(hoursDecimal)
    hoursDecimal -= h
    m = floor(hoursDecimal * 60)
    hoursDecimal -= m/60
    s = round(hoursDecimal * 3600, 2)
    print(f'''RA: {h}h {m}' {s}"''')

def printDMS(dec):
    d = -1 if dec < 0 else 1
    dec = abs(dec)
    d = d * floor(dec)
    dec -= abs(d)
    dec = abs(dec)
    m = floor(dec * 60)
    dec -= m/60
    s = round(dec * 3600, 2)
    print(f'''DEC: {d}deg {m}' {s}"''')

def gregorianToJulianDay(year, month, day, hour, minute, second):
    """
    Time is in UT!
    """
    jDay = 367 * year - (7 * (year + (month + 9) // 12)) // 4 + (275 * month) // 9 + day + 1721013.5 + (hour + minute / 60 + second / 3600) / 24 
    return jDay