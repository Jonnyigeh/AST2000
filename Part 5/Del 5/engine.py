# EGEN KODE

"""
Imports
"""
import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt
import scipy.constants as scc
import random as rn
import time
import ast2000tools as ast

"""
Constants
"""
L = 1e-6                                        # Length of each side in the box, [m]
L_e = L*0.8                                    #lenght of the exhaust hole       [m]
N = 1000                                      # Number of particles
T = 3*1e3                                       # Temperature                     [J]
mh = 1.00784*scc.m_u                            # Mass of hydrogenatom            [kg]
m = 2*mh                                        # Mass of H2 molecule             [kg]
dt_box_simulation = 1e-12                                       # Timestep                        [s]
timesteps = 1000                                  # Number of timesteps
dry_mass = 1100                                 # Mass of rocket                    [kg]
planet_rad = 7744.8*1e3                             # Radius planet                     [m]
planet_mass = 5.18367161e-06*1.9891*1e30        # Mass of planet                    [kg]
planet_rotvel = np.array((0, planet_rad*2*np.pi/(1.02244013*24*3600)))  # Rotational velocity of planet     [m/s]
"""
Functions
"""

def box(N, L, T, m):
    """
    This will create our box, with N particles distributed
    at random positions within our LxL cube with origin at the center(0,0,0).
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
    v_z = []                                                    #List for all z-components of particle velocities that escaped through the exhaust hole.
    num_esc = 0                                                 #Escaped particle tracker
    sigma = np.sqrt(scc.k*T/m)
    mu = np.sqrt(8*sigma/np.pi)
    p = 0                                                    # Total Momentum
    epsilon = 1e-5                                           # The epsilon
    for i in range(N):
        for j in range(timesteps-1):                         # 1000 timesteps
            if  r[i][j][2] < (-L/2 + epsilon) and abs(r[i][j][0]) < L_e/2 and abs(r[i][j][1]) < L_e/2:
                """
                In our if-statement, we have three conditions that must be fulfilled in order
                for the particle to escape the exhaust hole.
                i) the particle position must be on xy-plane
                ii) x-component of the vector must be between -L*0.05 and L*0.05, this is the dimension of the exhaust hole.
                iii) y-component of the vector must be between -L*0.05 and L*0.05, this is the dimension of the exhaust hole.

                Because of the timesteps, z-component of position vector "never" gets to zero value, it simply "jumps" over it.
                In order to solve this problem, we introduce an epsilon value, a value that is aproximated to be near zero.
                In this way, we can aproximate when the z-component of the vector is at xy-plane.
                """
                #print(f"{num_esc} particles just escaped!")
                v_z.append(abs(v[i][j][2]))
                r[i][j+1] = (rn.uniform(-L/2, L/2), rn.uniform(-L/2, L/2), rn.uniform(-L/2, L/2))
                v[i][j+1] = (rn.gauss(mu, sigma), rn.gauss(mu, sigma), rn.gauss(mu, sigma))
                """
                When we detect that our particle just escaped throught the exhaust hole, we simply reposition
                it at random place in our box with random (with respect to Gaussian distribution) velocity.
                """
                num_esc += 1
                #p += (m*abs(v[i][j][2]))
                """"
                We want our rocket to go upwards in z-axis.
                Therefore, we are interested in the velocities of the particles that are paralel to
                the rocket acceleration direction.
                From theory, we know that the momentum (p) is proportional to the mass of
                the object (in our case the particle) time the velocity.

                """
            else:
                pass
            for k in range(3):
                if abs(r[i][j][k]) > abs(L/2):                  # Checks if a given particle is outside the box (collision)
                    v[i][j+1][k]*(-1)                           # If so, will flip the velocity component of the dimension that collided.
                else:
                    pass
            r[i][j+1] = r[i][j] + v[i][j+1]*dt                  # Each particle will move v*dt per timestep
    t2 = time.time()
    p = m*np.mean(v_z)                                          # Find the mean value of p = mv
    Force = num_esc*p/(dt*timesteps)                            # Total force exerted over our dt*timestep, will be assumed constant
    fuel_loss = num_esc*m                                       # Loss of fuel by all particles leaving the tank per 10-9 s
    print("====================================================================================")
    print(f"Total escaped particles per second: {num_esc*1e9}")
    print(f"The Change in momentum: {p*num_esc*1e9}")
    print(f"Net Force exerted by escaped particles: {Force} N")
    print(f"Total fuel loss per second: {fuel_loss*1e9} kg")
    print(f"=============FINISHED IN {t2 - t1} SECONDS=================================")
    return r, v, t2-t1, fuel_loss*1e9, Force, np.mean(v_z)

def thrust_fuelloss(thrust_by_one_box, number_of_boxes, fuel_loss_per_s_one_box):
    """
    Calculates total thrust from superposition of all
    boxes on our rocket.
    Assume they all produce equal thrust
    Return thrust from entire engine + fuel consumed per second
    """
    thrust = thrust_by_one_box*number_of_boxes
    total_consumption = fuel_loss_per_s_one_box*number_of_boxes
    return thrust, total_consumption
def boost(initial_rocket_mass, engine_thrust, fuel_consumption, change_in_vel):
    """
    Force = change in momentum
    We are looking at dt = 1 s (per second), and all the parameteres are given in per sec
    """
    F = engine_thrust
    FC = fuel_consumption/10
    M_r = initial_rocket_mass
    delta_v = change_in_vel
    i = 0
    while (F * i) <= (M_r * delta_v):
        i += 1

    fuel_consumed = FC * i
    print(f"We burned {fuel_consumed}kg of fuel to boost our rocket for a change in velocity = {change_in_vel}!")
    return fuel_consumed


def numb_boxes(L):
    """
    Return the maximum number of boxes with
    sidelength L that can fit on the rocket.
    """
    area = 4                # [m]
    numb = int(4/L)
    return numb*numb

if __name__ == "__main__":
    r, v, t, loss, thrust, mean_v = simulation(dt_box_simulation)
    number_of_boxes = numb_boxes(L)
    engine_thrust, fuel_cons_per_s = thrust_fuelloss(thrust, number_of_boxes, loss)

    breakpoint()
