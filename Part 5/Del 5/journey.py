    # EGEN KODE
import numpy as np
import matplotlib.pyplot as plt
import ast2000tools.constants as atc
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts
from trajectory import Trajectory
from interpolation import PlanetaryOrbits

class Journey:
    def __init__(self, initial_year, initial_position, initial_velocity, amount_of_fuel = 1000, time_of_transfer = 0.05939681):
        self.current_year = initial_year
        self.current_position = initial_position
        self.current_velocity = initial_velocity
        self.current_fuel = amount_of_fuel
        self.time_of_destination = np.array([initial_year + time_of_transfer])
        self.dt = 0.00001
        self.time_of_journey = np.arange(initial_year, initial_year + 0.6, self.dt)
        self.t = self.time_of_journey
        self.planetary_positions = list(np.zeros((7, 2)))
        for i in range(7):
            self.planetary_positions[i] = PlanetaryOrbits().interpolation(i)
        self.system = SolarSystem(76408)


    def find_distance_to_destination(self):
        """
        Finner avstand til destinasjon, og retning
        """
        # Finding vector from sun to planet 0 and planet 6, and distance
        self.current_year_index = np.where(self.time_of_journey == self.current_year)
        self.r0 = np.array((self.planetary_positions[0][0](self.t[self.current_year_index]), self.planetary_positions[0][1](self.t[self.current_year_index])))
        self.r6 = np.array((self.planetary_positions[6][0](self.t[self.current_year_index]), self.planetary_positions[6][1](self.t[self.current_year_index])))
        # Finding vector, and unitvector from shuttle to planet 6
        self.roc_to_planet = np.zeros_like(self.current_position)
        self.roc_to_planet[0] = self.r6[0] - self.current_position[0]
        self.roc_to_planet[1] = self.r6[1] - self.current_position[1]
        self.l = np.linalg.norm(self.current_position) * np.sqrt(self.system.masses[6] / (10 * self.system.star_mass))

    def injection(self):
        vstable = np.sqrt(atc.G_sol * self.system.masses[6] / np.linalg.norm(self.roc_to_planet))
        theta = np.arccos(np.dot(np.array([1,0]), -self.roc_to_planet) / np.linalg.norm(self.roc_to_planet))
        e_theta = - np.sin(theta) * np.array([1,0]) + np.cos(theta) * np.array([0,1])
        e_r = np.cos(theta) * np.array([1,0]) + np.sin(theta) * np.array([0,1])
        v_pm = vstable * e_theta
        delta_vinj = v_pm - self.current_velocity
        self.deltav2 = self.Vaphe * e_theta - self.current_velocity
        self.boost(self.current_fuel + 1100, deltav = delta_vinj)




    def boost(self, initial_rocket_mass, deltav, engine_thrust = 9.518 * 10 **4, fuel_consumption = 3.389):
        """
        Force = change in momentum
        We are looking at dt = 1 s (per second), and all the parameteres are given in per sec
        """
        F = engine_thrust
        FC = fuel_consumption
        M_r_init = initial_rocket_mass               # [kg]
        tmp = self.current_velocity
        change_vel = deltav
        delta_vms = np.linalg.norm(change_vel) * 4740.57172               # omgjøres til m/s
        self.current_velocity = tmp + deltav
        i = 0
        M_r = M_r_init
        while (F * i) <= (M_r_init * delta_vms):
            i += 1
            M_r -= FC * i
        fuel_consumed = FC * i
        print(f"We burned {fuel_consumed}kg of fuel to boost our rocket for a change in velocity = {delta_vms}!")
        self.current_fuel -= fuel_consumed


    def part_of_voyage(self, number_of_timesteps):
        """
        Finds the trajectory from current year to current year + numberoftimestep * dt
        """
        instance = Trajectory(self.current_year, self.current_position,
                                        self.current_velocity,self.current_year+number_of_timesteps*self.dt, self.dt)
        x, v, a, t = instance.calculate_trajectory()
        self.instance = instance
        self.current_year = t[-1]
        self.current_position = x[-1]
        self.current_velocity = v[-1]
        self.find_distance_to_destination()
        if np.linalg.norm(self.roc_to_planet) < np.linalg.norm(self.l):
            print("We've reached orbit!!")

        tplot = np.linspace(0, 1, 10000)
        plt.plot(instance.interpolynomial[6][0](self.current_year), instance.interpolynomial[6][1](self.current_year), "k*") # Vi henter ut klasseobjekt PlanterayOrbit fra Trajectory koden
        plt.plot(instance.interpolynomial[0][0](self.t[0]), instance.interpolynomial[0][1](self.t[0]), "k*")
        plt.plot(x[-1,0], x[-1,1], "r*")
        plt.plot(instance.interpolynomial[0][0](tplot),instance.interpolynomial[0][1](tplot)
                        ,instance.interpolynomial[6][0](tplot), instance.interpolynomial[6][1](tplot))

        plt.plot(x[:,0], x[:,1], "r--") #, v, t, a, t)
        plt.legend(["Planet 0", "Planet 6", "Shuttle"])
        plt.title("Hohmann transfer")
        plt.show()

    def Hohmann(self):
        self.mu = atc.G_sol * self.system.star_mass
        r1 = self.system.semi_major_axes[0]
        r2 = self.system.semi_major_axes[6]
        alpha = np.pi * (1 - 1 / (2 * np.sqrt(2)) * np.sqrt((r1/r2) + 1) **3 )
        positions = self.planetary_positions
        theta = np.arccos(np.dot(np.array([1,0]), np.array((positions[0][0](self.t[0]), positions[0][1](self.t[0])))) \
                                            / np.linalg.norm(np.array((positions[0][0](self.t[0]), positions[0][1](self.t[0])))))

        e_theta = - np.sin(theta) * np.array([1,0]) + np.cos(theta) * np.array([0,1])
        self.T_transfer = 0.5 * ( (r1 + r2) / 2) ** (1.5)
        self.a_transfer = (r1 + r2) / 2
        self.P_transfer = np.sqrt( (4 * np.pi**2 * self.a_transfer**3)/ self.mu ) # yr
        self.Vperi = ( 2 * np.pi * self.a_transfer / self.P_transfer ) * np.sqrt( 2 * self.a_transfer / r1 - 1)
        Vp = np.sqrt(2 * self.mu * (r2 / (r1 * (r1 + r2)))) * e_theta
        Va = r1 / r2 * Vp
        self.Vaphe = r1 / r2 * self.Vperi
        self.V_p0 = (2 * np.pi * r1) / (self.system.rotational_periods[0] / 365)        # AU/yr
        self.V_p6 = (2 * np.pi * r2) / (self.system.rotational_periods[6] / 365)
        self.deltav1 = self.Vperi*e_theta - self.current_velocity
        self.deltav2 = self.Vaphe * self.deltav1 / np.linalg.norm(self.deltav1)
        self.boost(self.current_fuel + 1100, deltav= self.deltav1)
        self.part_of_voyage(4250)
        self.boost(self.current_fuel + 1000, deltav=self.deltav2)

    def check_gravpull(self):
        self.find_distance_to_destination()
        planet6 = np.array((self.instance.interpolynomial[6][0](self.current_year), self.instance.interpolynomial[6][1](self.current_year)))
        dis_roc_to_plan6 = np.linalg.norm(planet6 - self.current_position)
        if dis_roc_to_plan6 < self.l:
            print("We are in orbit!!")
            return
        print(f"We are NOT in orbit. distance is {dis_roc_to_plan6}")
        return

if __name__ == "__main__":
    """
    VErdier for Hohmann transfer, altså init pos etter launch ved Hohmann kriterier
    """
    time_after_launch = 2.059 + 1.0234195e-5
    instance = Journey(time_after_launch, np.array([0.09960206, 0.07912785]), np.array([-8.97560312, 10.3206985]))
    instance.Hohmann()
