from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts
import numpy as np

mission = SpaceMission(76408)
system = SolarSystem(76408)
shortcuts = SpaceMissionShortcuts(mission, [45535])

"""
Following values are from the engine.py program
the performance of our rocket engine.

This is the shortcut given to verify the launch results using our own engines
performance
"""

Thrust = 9.518 * 10 **4
fuel_loss_per_second = 3.389
initial_fuel_mass = 1500
launch_duration = 450
x0 = 0.1280914683 + 7744.8003307 * 10 ** 3 * 6.68458712 * 10** (-12)      # [AU]
y0 = 0
t = 0
mission.set_launch_parameters(Thrust, fuel_loss_per_second, initial_fuel_mass, launch_duration, (x0, y0), t)
mission.launch_rocket()
consumed_fuel_mass, final_time, final_position, final_velocity = shortcuts.get_launch_results()
mission.verify_launch_result(final_position)
