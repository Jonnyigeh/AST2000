#EGEN KODE
import ast2000tools.constants as atc
import numpy as np
from ast2000tools.solar_system import SolarSystem

sys = SolarSystem(76408)
rho_0 = sys.atmospheric_densities[6]            # surface density in kg/m^3
M = sys.masses[6] * atc.m_sun                   # planet mass in kg
R = sys.radii[6] * 1000                         # planet radius in meters
m = 90                                          # lander mass in kg
G = atc.G                                       # gravitational constant in SI units
v_t = 3                                         # soft-landing threshold velocity in m/s
C_d = 1                                         # Constant

def cross_area():
    A = 2 * G * m * M / (R ** 2 * v_t ** 2 * rho_0 * C_d)
    return A

print(cross_area())
breakpoint()
