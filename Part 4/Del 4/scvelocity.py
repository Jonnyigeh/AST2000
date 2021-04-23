# EGEN KODE
import numpy as np
import ast2000tools.constants as astc
from ast2000tools.space_mission import SpaceMission

class SCVEL:
    """
    Calculates spaceshuttles velocity by look at the doppler shift
    from two distant stars that we already know the radial velocity
    and also the dopplershift from a inertia system, so we look at the
    difference in the known dopplershift from the measured one.

    We use lambda_0 as the spectralline for hydrogen
    """
    def __init__(self, seed):
        self.lambda_0 = 656.3                      # Spectralline hydrogen [nm]
        self.mission = SpaceMission(seed)
        Doppler_shifts = self.mission.star_doppler_shifts_at_sun            # Tuple of doppler shifts [nm]
        Star_angles = self.mission.star_direction_angles                    # Tuple of aximuthal angle of ref.stars relative to solarsystem x-axis [deg]
        self.ds_star1, self.ds_star2 = Doppler_shifts
        self.ang_star1, self.ang_star2 = Star_angles
        #Positions_after_launch = 0.128154526, 1.14510588e-04                # Posisjon etter launch relativ til sol
        #self.x0, self.y0 = Positions_after_launch

    def wavelen_to_radvel(self, lambda_1=0, lambda_2=0):
        """
        Transforms measured dopplershift into
        radial velocity by using the formula
        for doppler shift. Set to zero by default

        radial velocity as measured by sun ref.sys.
        Saves into attributes that can be accessed using property funcs.
        """
        def formula(delta_lambda, lambda_0):
            return delta_lambda * astc.c_AU_pr_yr / lambda_0
        self._radvel_star1 = formula(self.ds_star1, self.lambda_0)
        self._radvel_star2 = formula(self.ds_star2, self.lambda_0)
        self.radvel_1 = formula(lambda_1, self.lambda_0)
        self.radvel_2 = formula(lambda_2, self.lambda_0)

    def find_vel(self):
        """
        Find velocity along x, and y axis from
        doppler shift along the phi_1, and phi_2 axis
        from ref.star angles.
        Defines new coord.axis uhat1 and uhat2 from these
        angles using the formulas in hint 3 proj4.

        Saves x-vel and y-vel to arrays can be accessed by property.
        Galilean coord.transformation.
        """
        phi_1 = self.ang_star1 * np.pi/180
        phi_2 = self.ang_star2 * np.pi/180
        uhat_1 = np.array((np.cos(phi_1), np.sin(phi_1)))
        uhat_2 = np.array((np.cos(phi_2), np.sin(phi_2)))

        phivel_1 = (self._radvel_star1 - self.radvel_1)
        phivel_2 = (self._radvel_star2 - self.radvel_2)
        phivel = np.array(([phivel_1], [phivel_2]))
        transformation_mat = 1 / np.sin(phi_2 - phi_1) * np.array(([np.sin(phi_2), -np.sin(phi_1)], [-np.cos(phi_2), np.cos(phi_1)]))
        self._xvel, self._yvel = np.matmul(transformation_mat, phivel)

    @property
    def radvel_star1(self):
        """
        Returns stars radial velocity in m/s
        """
        return self._radvel_star1
    @property
    def radvel_star2(self):
        """
        Returns stars radial velocity in m/s
        """
        return self._radvel_star2
    @property
    def xvel(self):
        """
        Returns velocity along x-xis in AU/yr
        """
        return self._xvel
    @property
    def yvel(self):
        """
        Returns velocity along y-axis in AU/yr
        """
        return self._yvel

if __name__ == "__main__":
    instance = SCVEL(76408)
    instance.wavelen_to_radvel(0.10236087076815302, -0.0873910486582896) # Disse er verdiene er tatt fra SpaceMission.measure_star_doppler_shifts()
    instance.find_vel()
    breakpoint()
