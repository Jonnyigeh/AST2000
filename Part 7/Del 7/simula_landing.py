#EGEN KODE
import numpy as np
import ast2000tools.constants as atc
from atmosphere import Model
from ast2000tools.solar_system import SolarSystem

class LandingSequence:
    def __init__(self, A, x, v):
        """
        Parameters: - Cross-sectional area A
                    - Initial position x
                    - Initial velocity v
        """
        sys = SolarSystem(76408)
        self.Cd = 1
        self.A = A                                                  # Crossectional area
        self.planet_radius = sys.radii[6] * 1000                    # radius i meter
        self.M = sys.masses[6] * atc.m_sun
        self.rot_period = sys.rotational_periods[6] * 24 * 3600     # Rotationsperiode T i sekunder
        self.omega0 = 2 * np.pi / self.rot_period
        self.initial_position = x
        self.intiial_velocity = v

    def airdrag(self, vdrag, area = self.A):
        Fd = 0.5 * rho * self.Cd * self.A * vdrag ** 2
        return Fd

    def gravity(self, r, m = 90):
        Fg = atc.G * m * self.M / np.linalg.norm(r) ** 2
        return Fg

    def landing(self):
        """
        Atmosfærisk hastighet i theta, og radiell vektorer
            Vanlig hastighet i kartesisk, må omgjøres

        Foreløbi sketch av differensialligning, trenge bare
            å få nåe plot me kan tolka for poenguttelling i guess
        """
        dt = 1
        t = 0
        xpos = self.initial_position
        vel = self.initial_velocity
        while np.linalg.norm(xpos) > self.planet_radius:
            atm_vel = self.omega0 * np.linalg.norm(xpos) * np.array([0,1])
            vdrag = vel - atm_vel
            a = (airdrag(vdrag) - gravity(xpos)) / m
            vel += a * dt
            xpos += vel * dt
            t += dt



if __name__ == "__main__":
