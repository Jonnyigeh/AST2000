# EGEN KODE
import ast2000tools.constants as atc
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
from scipy.interpolate import interp1d




class Model:
    def __init__(self):
        system = SolarSystem(76408)
        self.gamma = 1.4
        self.T_surface = 350                                    # Deduced in part 3 by look at luminosity from sun [K]
        self.planet_mass = system.masses[6] * atc.m_sun         # Mass in kg
        self.planet_radius = system.radii[6] * 1000             # Radii in meters
        self.mu = 17                                            # Ei linja H2O, to linje CH4
        self.hydromass = 1.00784 * 1.661e-27                    # Mass of hydrogencore, in kg
        self.rho_surface = system.atmospheric_densities[6]      # Overflatetetthet i kg/m^3

    def Temperature(self, height):
        Temp = self.temperature_profile(height)
        return Temp

    def Density(self, height):
        return self.density_profile(height)

    def gravity(self, height):
        return atc.G * self.planet_mass / (self.planet_radius + height) ** 2

    def create_profiles(self):
        dr = 1                                                  # Small dr [m]
        r = [0]
        height = 150000
        rho = [self.rho_surface]
        temp = [self.T_surface]
        while r[-1] < height:
            if temp[-1] > self.T_surface / 2:                                                                                               # Adiabatic layer
                dTdr = -(self.gamma - 1) * self.mu * self.hydromass * self.gravity(r[-1]) / (atc.k_B * self.gamma)
                T = temp[-1] + dTdr * dr
                drhodr = - rho[-1] / temp[-1] * dTdr - rho[-1] / temp[-1] * self.gravity(r[-1]) * self.mu * self.hydromass / atc.k_B
            else:
                drhodr = - rho[-1] / (self.T_surface / 2) * self.gravity(r[-1]) * self.mu * self.hydromass / atc.k_B                        # Isothermal Layer
                T = self.T_surface / 2
            p = rho[-1] + drhodr * dr

            temp.append(T)
            rho.append(p)
            r.append(r[-1] + dr)
        rhoarray = np.array(rho)
        temparray = np.array(temp)
        self.density_profile = interp1d(np.array(r), rhoarray)#, kind="nearest")
        self.temperature_profile = interp1d(np.array(r), temparray, kind="quadratic")




if __name__ == "__main__":
    instance = Model()
    instance.create_profiles()
    r = np.arange(0,150000, 1)
    plt.subplot(121)
    plt.plot(r, instance.Temperature(r) / instance.T_surface)
    plt.title("Temperatureprofile")
    plt.xlabel("Height [m]")
    plt.ylabel("Temperature/Tmax [K]")
    plt.subplot(122)
    plt.plot(r, instance.Density(r) / instance.rho_surface)
    plt.title("Densityprofile")
    plt.ylabel("Density/rho_max [kg/m^3]")
    plt.xlabel("Height [m]")
    plt.tight_layout()
    plt.show()

    # x = np.arange(0, 100000, 1)
    # plt.plot(x, instance.Temperature(x))
    # plt.show()
