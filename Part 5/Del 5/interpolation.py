import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
from ast2000tools.shortcuts import SpaceMissionShortcuts


class PlanetaryOrbits:
    """
    Interpolerer verdier fra snarvei fra 0 til 20 år
    """
    def __init__(self):
        self.planets = ( SpaceMissionShortcuts(SpaceMission(76408),
                        [1925]).compute_planet_trajectories(np.linspace(0, 20, 10000)) )


    def interpolation(self, planet_index):
        t = np.linspace(0, 20, 10000)
        self.fx = interp1d(t, self.planets[0][planet_index])
        self.fy = interp1d(t, self.planets[1][planet_index])
        return (self.fx, self.fy)

if __name__ == "__main__":
    instance = PlanetaryOrbits()
    t = np.linspace(0, 20, 10000)
    fx6, fy6 = instance.interpolation(6)
    fx, fy = instance.interpolation(0)
    plt.plot(fx6(t),fy6(t),fx(t),fy(t))
    breakpoint()
