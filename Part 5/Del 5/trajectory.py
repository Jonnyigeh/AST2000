import ast2000tools.constants as atc
import numpy as np
import matplotlib.pyplot as plt
import sys

from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
from ast2000tools.shortcuts import SpaceMissionShortcuts
from interpolation import PlanetaryOrbits


class Trajectory:
    def __init__(self, initial_time, initial_pos, initial_vel, simulation_time, dt):
        """
        Input:
        - initial_time [yr]
        - initial_pos [AU]
        - initial_vel [AU/yr]
        - timesteps = number of timesteps
        - dt = step length of timesteps
        Output:
        - final time [yr]
        - final position [AU]
        - final velocity [AU/yr]
        """
        self.system = SolarSystem(76408)
        self.mission = SpaceMission(76408)
        self.initial_time = initial_time
        self.initial_pos = initial_pos
        self.initial_vel = initial_vel
        self.simulation_time = simulation_time
        self.dt = dt
        self.G = atc.G_sol # I AU , yr og solmasser
        self.planet_masses = self.system.masses #* atc.m_sun # m_sun er solens masse, fordi solarsystem oppgir i solmasser
        self.M_s = self.system.star_mass #* atc.m_sun


    def calculate_trajectory(self):
        """
        Uses Leapfrog integration method to solve differentail equation
        for motion of the shuttle traverseing the solarsystem.
            Forces acting on the shuttle from all planets + star
            Uses the interpolated positions of the planets to calculate the forces using
                Newtons law of gravity.
        """
        t = np.arange(self.initial_time, self.simulation_time, self.dt)
        x = np.zeros((len(t), 2))                               # Posisjon sonde
        v = np.zeros_like(x)                                    # Hastighet sonde
        a = np.zeros_like(x)
        p6 = np.zeros_like(x)
        instance = PlanetaryOrbits()
        self.interpolynomial = list(np.zeros((7, 2)))
        for i in range(7):
            self.interpolynomial[i] = instance.interpolation(i)
        # def r_(j,i):
        #     return np.array(self.planet_distances[0][j][i], self.planet_distances[1][j][i])
        def r_(j,i):
            return np.array((self.interpolynomial[j][0](i), self.interpolynomial[j][1](i)))
        def F(i):
            F_tot = np.zeros(2)
            for j in range(7):
                F_tot += self.G * self.planet_masses[j] * (x[i] - r_(j,t[i])) / (np.linalg.norm(x[i] - r_(j, t[i])) ** 3)
            F_tot += -self.G * self.M_s * x[i] / (np.linalg.norm(x[i]) ** 3)
            return F_tot

        x[0] = self.initial_pos
        v[0] = self.initial_vel
        for i in range(len(t)-1):
            a[i] = F(i)
            x[i+1] = x[i] + v[i] * self.dt + 0.5 * a[i] * self.dt ** 2
            a_1 = F(i+1)
            v[i+1] = v[i] + 0.5 * (a[i] + a_1) * self.dt
        return x, v, a, t


if __name__ == "__main__":
    instance = Trajectory(2.11840704, np.array((-0.214116, -0.111399)), np.array((4.49143,-6.20743)), 15, 0.001)
    x, v, a, t = instance.calculate_trajectory()
    plt.plot(x[:,0], x[:,1]) #, v, t, a, t)
    plt.legend(["Position"])
    plt.axis([0, 20, 0, 30])
    plt.show()
