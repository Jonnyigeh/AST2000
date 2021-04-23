# EGEN KODE

"""
Imports
"""
import numpy as np
import matplotlib.pyplot as plt
import ast2000tools.constants as c
import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import random
seed = 76408 #Seed for our solarsystem
system = SolarSystem(seed)


class Planetery_orbits():
    def __init__(self):
        self.timesteps = 10000
        self.dt = float(1/self.timesteps)
        self.m_s = system.star_mass #The mass of our sun
        self.G = 4 * np.pi**2



    def information(self):
        system.print_info()
        print(f"Rotation period of the planets 6 {system.rotational_periods[6]}")
        print(f"The surface temperature of the sun is {system.star_temperature}K")
        print(f"My solar system has {system.star_mass} solar mass with radius of {system.star_radius} km")
        for i in range(system.number_of_planets):
            print(f"Planet {i} has semi-major axis of {system.semi_major_axes[i]}AU, and eccentricity of {system.eccentricities[i]}AU")



    def analytical(self,N,i):
        f = np.linspace(0, 2*np.pi, N)
        for i in range(7):
            x = system.semi_major_axes[i] * np.cos(f)
            y = (system.semi_major_axes[i] * np.sqrt(1 - system.eccentricities[i]**2) * np.sin(f))
        plt.plot(x,y,"--", label=f"planet {i}")
        plt.plot(0,0,"*")
        plt.title("Analytical solution for planet orbits ")
        plt.legend()
        plt.xlabel("x-axis in AU")
        plt.ylabel("y-axis in AU")


    def numerical(self, x0, v0, T, N):

        def a(i):
            r = x[i]
            return (-self.G*self.m_s)*(r/np.linalg.norm(r)**3)

        dt = (T/N)
        t = np.zeros(N, float)
        v = np.zeros((N,2), float)
        x = np.zeros((N,2),float)
        t[0] = 0
        x[0] = x0
        v[0] = v0
        a_i = a(0)
        for i in range(N-1):
            t[i+1] = t[i] + dt
            x[i+1] = x[i] + v[i]*dt + 0.5*a_i*dt**2
            a_iplus1 = a(i)
            v[i+1] = v[i] +  0.5*(a_i + a_iplus1)*dt
            a_i = a_iplus1


        plt.plot(x[:,0], x[:,1])
        plt.plot(0,0,"*")
        plt.title("Numerical solution for planet orbit")
        plt.legend()
        plt.xlabel("x-axis in AU")
        plt.ylabel("y-axis in AU")

    def wobbling(self, x_s0, v_s0, x0, v0, T, N):


        dt = (T/N)
        t = np.zeros(N, float)
        self.v_s = np.zeros((N,2), float)
        x_s = np.zeros((N,2),float)
        t[0] = 0
        x_s[0] = x_s0
        self.v_s[0] = v_s0



        v = np.zeros((N,2), float)
        x = np.zeros((N,2),float)
        x[0] = x0
        v[0] = v0
        def a_s(i):
            r = x[i] - x_s[i]
            return (self.G*system.masses[5])*(r/np.linalg.norm(r)**3)
        def a(i):
            r = x_s[i] - x[i]
            return (self.G*self.m_s)*(r/np.linalg.norm(r)**3)
        a_i_s = a_s(0)
        a_i = a(0)

        for i in range(N-1):
            t[i+1] = t[i] + dt
            x[i+1] = x[i] + v[i]*dt + 0.5*a_i*dt**2
            x_s[i+1] = x_s[i] + self.v_s[i]*dt + 0.5*a_i_s*dt**2
            a_iplus1 = a(i)
            a_iplus1_s = a_s(i)
            v[i+1] = v[i] +  0.5*(a_i + a_iplus1)*dt
            self.v_s[i+1] = self.v_s[i] +  0.5*(a_i_s + a_iplus1_s)*dt
            a_i = a_iplus1
            a_i_s = a_iplus1_s

        plt.subplot(1, 2, 1)
        plt.plot(x[:,0], x[:,1], label = "planet")
        plt.plot(x_s[:,0], x_s[:,1], label="sun")
        plt.legend()
        plt.xlabel("x-axis in AU")
        plt.ylabel("y-axis in AU")
        plt.subplot(1, 2, 2)
        plt.plot(x[:,0], x[:,1], label = "planet")
        plt.plot(x_s[:,0], x_s[:,1], label="sun")
        plt.title("zoomed into the sun")
        plt.legend()
        plt.xlabel("x-axis in AU")
        plt.ylabel("y-axis in AU")
        plt.show()

        plt.plot(t,self.v_s[:,0])
        plt.title("velocity of the sun in x-direction")
        plt.legend()
        plt.xlabel("time in AU")
        plt.ylabel("Velocity in AU/yr")

        plt.show()

    def construct_curve(self):
        t_0 = 30000
        mean = 0
        vmax = 26
        sigma = 0.2*vmax
        støy = np.zeros(t_0)
        for i in range(t_0-1):
            støy[i] = random.gauss(mean,sigma)
        t = np.linspace(0,t_0,30000)
        v_r = vmax*np.cos(t*2*np.pi/10957)
        self.v_obs = v_r + støy

        def function(vstar, P, t0):
            return vstar*np.cos(2*np.pi/P * (t - t0))

        plt.plot(t, self.v_obs, label="observed")
        plt.plot(t, v_r, "r", label="actuall")
        plt.plot(t, function(26.1111, 11000, 11000), "g" , label="least squared")


        plt.title("The Radial Velocity Cure")
        plt.xlabel("Observation Time in days")
        plt.ylabel("Radial Velocity in km/sec")
        plt.legend()
        plt.show()

    def solution(self):
        N = 10
        v_star = np.linspace(15, 35, N)                # km/s
        P = np.linspace(10000, 13000, N)
        t0 = np.linspace(9000, 11000, N)
        T = 30000
        def approxcurve(vstar, P, t0):
            t = np.linspace(0, T, T)
            def vr(t):
                return vstar*np.cos(2*np.pi/P * (t - t0))
            vr = vr(t)
            S = sum((self.v_obs[i] - vr[i])**2 for i in range(T))
            return S
        vr = self.v_obs
        S = 1000000
        vs = 0
        Po = 0
        t0o = 0
        for i in v_star:
            for j in P:
                for k in t0:
                    Snew = (approxcurve(i, j, k))
                    if Snew < S:
                        S = Snew
                        vs = i; Po = j; t0o = k
        def function(vstar, P, t0):
            return vstar*np.cos(2*np.pi/P * (t - t0))
        t = np.linspace(0,T, T)
        plt.plot(t, function(vs, Po, t0o))
        plt.show()


    def upgrade(self, i1, i2, i3, T, N):
        x10 = np.array([system.initial_positions[0,i1], system.initial_positions[1,i1]])
        x20 = np.array([system.initial_positions[0,i2], system.initial_positions[1,i2]])
        x30 = np.array([system.initial_positions[0,i3], system.initial_positions[1,i3]])

        v10 = np.array([system.initial_velocities[0,i1], system.initial_velocities[1,i1]])
        v20 = np.array([system.initial_velocities[0,i2], system.initial_velocities[1,i2]])
        v30 = np.array([system.initial_velocities[0,i3], system.initial_velocities[1,i3]])


        dt = (T/N)
        t = np.zeros(N, float)
        t[0] = 0

        x1 = np.zeros((N,2),float)
        x2 = np.zeros((N,2),float)
        x3 = np.zeros((N,2),float)
        xs = np.zeros((N,2),float)


        v1 = np.zeros((N,2), float)
        v2 = np.zeros((N,2), float)
        v3 = np.zeros((N,2), float)
        vs = np.zeros((N,2),float)


        x1[0] = x10
        x2[0] = x20
        x3[0] = x30
        xs[0] = np.array([0,0])


        v1[0] = v10
        v2[0] = v20
        v3[0] = v30
        vs[0] = np.array([0,0])



        def a_s(i, j):
            r = j - xs[i]
            return (self.G*system.masses[5])*(r/np.linalg.norm(r)**3)
        def a1(i):
            r = xs[i] - x1[i]
            return (self.G*self.m_s)*(r/np.linalg.norm(r)**3)
        def a2(i):
            r = xs[i] - x2[i]
            return (self.G*self.m_s)*(r/np.linalg.norm(r)**3)
        def a3(i):
            r = xs[i] - x3[i]
            return (self.G*self.m_s)*(r/np.linalg.norm(r)**3)


        a_s = a_s(0, x1[0]) + a_s(0, x2[0]) + a_s(0, x3[0])
        a1 = a1(0)
        a2 = a2(0)
        a3 = a3(0)

        for i in range(N-1):
            t[i+1] = t[i] + dt

            x1[i+1] = x1[i] + v1[i]*dt + 0.5*a1*dt**2
            x2[i+1] = x2[i] + v2[i]*dt + 0.5*a2*dt**2
            x3[i+1] = x3[i] + v3[i]*dt + 0.5*a3*dt**2
            xs[i+1] = xs[i] + vs[i]*dt + 0.5*a_s*dt**2

            a1_iplus1 = a1(i)
            a2_iplus1 = a2(i)
            a3_iplus1 = a3(i)
            a_s_iplus1 = a_s(i)

            v1[i+1] = v1[i] +  0.5*(a1 + a1_iplus1)*dt
            v2[i+1] = v2[i] +  0.5*(a2 + a2_iplus1)*dt
            v3[i+1] = v3[i] +  0.5*(a3 + a2_iplus1)*dt
            vs[i+1] = vs[i] +  0.5*(a_s + a_s_iplus1)*dt

            a1 = a1_iplus1
            a2 = a2_iplus1
            a3 = a3_iplus1
            a_s = a_s_iplus1





















if __name__ == '__main__':
    # x0 = np.array([system.initial_positions[0,0], system.initial_positions[1,0]])
    # v0 = np.array([system.initial_velocities[0,0],system.initial_velocities[1,0]])
    # x_s0 = np.array([0,0])
    # v_s0 = np.array([0,0])
    inst = Planetery_orbits()
    inst.information()
    # inst.analytical(10000, 0)
    # inst.numerical(x0, v0, 1, 1000000)
    #inst.wobbling(x_s0, v_s0, x0, v0, 0.5, 100000)
    #inst.construct_curve()
    #inst.solution()
    # inst.upgrade(0,1,1, 1, 1000)
