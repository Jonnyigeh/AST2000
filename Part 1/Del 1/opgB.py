# EGEN KODE

"""
Imports
"""
import numpy as np
import numba as nb
import scipy.integrate as scint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants as scc
import random as rn
import time

"""
Constants
"""
L = 1e-6                                        # Length of each side in the box, [m]
N = 100                                        # Number of particles
T = 3*1e3                                       # Temperature                     [J]
mh = 1.00784*scc.m_u                            # Mass of hydrogenatom            [kg]
m = 2*mh                                        # Mass of H2 molecule             [kg]
vavg = 6.3*1e3
dt = 1e-9                                       # Timestep                        [s]
timesteps = 1000                                 # Number of timesteps

"""
Functions
"""

def box(N, L, T, m):
    """
    This will create our box, with N particles distributed
    at random positions within our LxLxL cube with origin at the center(0,0,0).
    Also every particle will be given a velocity following
    a Gaussian distribution with parameteres mu and sigma.

    We define two arrays, with nested Nx3 matrices that contain position- and velocity vectors
    for each of the N particles.
    Returns the two arrays.
    """
    sigma = np.sqrt(scc.k*T/m)
    mu = np.sqrt(8*sigma/np.pi)
    position = np.zeros((N, timesteps, 3))
    velocity = np.zeros_like(position)
    for i in range(N):
        position[i][0] = (rn.uniform(-L/2, L/2), rn.uniform(-L/2, L/2), rn.uniform(-L/2, L/2))
        velocity[i] = (rn.gauss(mu, sigma), rn.gauss(mu, sigma), rn.gauss(mu, sigma))
    return position, velocity



def simulation(dt):
    """
    This function will simulate the particles moving inside the box
    by setting it up as a differential equation knowing that r = d/dt(v) with a constant timestep dt.
    We assume that the acceleration of the molecules is 0. (no interaction between molecules, and elastic wall collisions)
    We'll use ForwardEuler to compute the differential equation.
    Will also return time spent computing.
    """
    t1 = time.time()
    r, v = box(N,L,T,m)
    for i in range(N):
        for j in range(timesteps-1):                            # 1000 timesteps
            for k in range(3):
                if abs(r[i][j][k]) > abs(L/2):                  # Checks if a given particle is outside the box (collision)
                    v[i][j+1][k]*(-1)                           # If so, will flip the velocity component of the dimension that collided.
                else:
                    pass
            r[i][j+1] = r[i][j] + v[i][j+1]*dt                  # Each particle will move v*dt per timestep
    t2 = time.time()
    return r, v, t2-t1

if __name__ == "__main__":
    x, v, t = simulation(dt)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for i in range(0, 100, 10):
        ax.plot(x[i][0],x[i][1],x[i][2])
    plt.show()
    breakpoint()
