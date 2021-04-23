#EGEN KODE

"""
Prosjekt del 1
"""
import numpy as np
import numba as nb
import ast2000tools as ast
from scipy import integrate as int

def probfunc(x, mu, sigma):
    """
    The probability function with its 2 function parameters
    and the variable x.
    Here mu = the mean and sigma = standard deviation of the
    Gaussian distribution we are examining.
    """
    PF = 1/(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5*((x - mu)/sigma)**2)
    return PF

def P(a, b, mu, sigma):
    """
    A function that evaluates the probabilty of
    finding x, within the interval [a, b] by using
    the Gaussian probability function.

    We use scipy.integrate module and the Simpsonsmethod
    as this is more efficient and also relatively accurate way
    of numerical integration.
    """
    interv = np.linspace(a, b, 1000)
    P = int.simps(probfunc(interv, mu, sigma), interv)
    return P
"""
Task a2
The function P(a<x<b) is an integral that will give a value in [0, 1]
which shows the probability that a given, arbitrary x, will exist inside the interval
given by [a, b].
"""

#Task a3
"""
Here we simplify by using both sigma and mu = 1.
We will still see that the expected values of 68%, 95% and 99.7%
will be printed to our terminal regardless of parametervalues.
"""
sigma = 1
mu = 1
print(P(mu-sigma, mu+sigma, mu, sigma))
print(P(mu-2*sigma, mu+2*sigma, mu, sigma))
print(P(mu-3*sigma, mu+3*sigma, mu, sigma))

#TASK a4
"""
Se skriftark for utledning
"""
