from ast2000tools.shortcuts import SpaceMissionShortcuts
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
from imganaly import BA
from scvelocity import SCVEL
from position import Position
import Generated_pics
import numpy as np

mission = SpaceMission(76408)
system = SolarSystem(76408)

Thrust = 9.518 * 10 **4
Thrust_per_box = 6.189 * 10 ** (-9)
fuel_loss_per_box = 2.204 * 10 ** (-12)
fuel_loss_per_second = 3.389
launch_duration = 434
initial_fuel_mass = 1500
x0 = 0.1280914683 + 7744.8003307 * 10 ** 3 * 6.68458712 * 10** (-12)      # [AU]
y0 = 0
time_of_launch = 0

# Verify launch results needed to run the orientation verification
shortcuts = SpaceMissionShortcuts(mission, [45535])
mission.set_launch_parameters(Thrust, fuel_loss_per_second, initial_fuel_mass, launch_duration, (x0, y0), time_of_launch)
mission.launch_rocket()
consumed_fuel_mass, final_time, final_position, final_velocity = shortcuts.get_launch_results()
mission.verify_launch_result(final_position)

# Finds angle using the 360 refrence pictures
mission.take_picture("finalimg.png")
Imga = BA("finalimg.png")
angle = Imga.finding_angle()

# Fins shuttles velocity relative to star usign dopplershift
scvel = SCVEL(76408)
lambda_1, lambda_2 = mission.measure_star_doppler_shifts()
scvel.wavelen_to_radvel(lambda_1, lambda_2)
scvel.find_vel()
xvel = scvel.xvel
yvel = scvel.yvel

#Position
distances = mission.measure_distances()
instance = Position(final_time, distances)
xpos, ypos = instance.find_position()

#Verifying manual orientation
mission.verify_manual_orientation(np.array([xpos, ypos]), np.array([xvel[0], yvel[0]]), angle)

if __name__ == "__main__":
    pass
