from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts
mission = SpaceMission(76408)
system = SolarSystem(76408)

time_of_launch = 19.4
Thrust = 9.518 * 10 **4
Thrust_per_box = 6.189 * 10 ** (-9)
fuel_loss_per_box = 2.204 * 10 ** (-12)
fuel_loss_per_second = 3.389
launch_duration = 434
initial_fuel_mass = 1500
x0 = 0.1280914683 + 7744.8003307 * 10 ** 3 * 6.68458712 * 10** (-12)      # [AU]
y0 = 0
shortcuts = SpaceMissionShortcuts(mission, [45535])
mission.set_launch_parameters(Thrust, fuel_loss_per_second, initial_fuel_mass, launch_duration, (x0, y0), time_of_launch)
mission.launch_rocket()
consumed_fuel_mass, final_time, final_position, final_velocity = shortcuts.get_launch_results()

mission.verify_launch_result(final_position)

mission.take_picture("finalimg.png")
Imga = BA("finalimg.png")
angle = Imga.finding_angle()



mission.verify_manual_orientation(final_position, final_velocity, angle)




travel = mission.begin_interplanetary_travel()
# travel.coast_until_time(2.1090023)
# travel.boost([1.33448142, 0.11111635])
# travel.take_picture('testing.png')
# travel.coast_until_time(2.11149023)
# travel.boost([-0.04120767, 0.6119064])
# travel.coast_until_time(2.11369)
# travel.take_picture('another_one.png')
time, position, velocity = travel.orient()
