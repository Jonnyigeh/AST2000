#EGENKODE
import matplotlib.pyplot as plt
import numpy as np
import numba as nb
import scipy.integrate as scint
import scipy.constants as scc

mh = 1.00784*scc.m_u
N = 1e5
T = 3000

vx = np.linspace(-2.5, 2.5, int(N))*1e4
@nb.jit
def Pvx(vx, m, T):
    """
    The function will return the probability
    density for velocity on a given range vx
    """
    P = np.sqrt(m/(2*np.pi*scc.k*T))*np.exp(-0.5*(m*vx**2)/(scc.k*T))
    return P

@nb.jit
def Pv(v, m, T):
    """
    Absolute velocity distribution
    """
    P = (m/(2*np.pi*scc.k*T))**(1.5)*np.exp(-0.5*m*v**2/(scc.k*T))*4*np.pi*v**2
    return P


probvx = scint.simps(Pvx(vx, mh, T), np.linspace(5,30,int(N))*1e3)
print(probvx)   # The probability that an arbitrary particle has a velocity
                # in the interval.

print(probvx*N)   # This will be a measure of the number of particles with
                  # a velocity inside the interval we used to define prob.
v = np.linspace(0, 3, int(N))*1e4
plt.plot(v, Pv(v, mh, T))
plt.ticklabel_format(axis = "both", style = "sci", scilimits=(0,0))
plt.plot(vx, Pvx(vx, mh, T))
plt.show()
