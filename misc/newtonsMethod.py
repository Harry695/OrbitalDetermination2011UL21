from math import *

def solveUsingNewtonsMethod(xGuess, f, fPrime):
    # catch f'(xGuess) = 0
    if fPrime(xGuess) == 0:
        xGuess += 0.1
    # update xGuess
    xNew = xGuess - f(xGuess) / fPrime(xGuess)
    print(xNew) # debug
    if abs(xNew - xGuess) > 1.0e-18:
        xNew = solveUsingNewtonsMethod(xNew, f, fPrime) # try again if not within bound yet
    return xNew

# testing functions
def f1(x):
    return 2 * x ** 2 - 5 * x - 3

def f1Prime(x):
    return 4 * x - 5

print(solveUsingNewtonsMethod(-1, f1, f1Prime))