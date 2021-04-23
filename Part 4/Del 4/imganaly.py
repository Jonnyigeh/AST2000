# EGEN KODE
import numpy as np
from ast2000tools.space_mission import SpaceMission
from imgmanip import BM
from PIL import Image

class BA(BM):
    """
    Class to analyze a picture given to decide
    what direction our camera is facing assuming it is
    in the xy-plane (theta_0 = 90 deg).
    """
    def __init__(self, bildefil, FOV=70):
        super().__init__(bildefil, FOV=70)

    def finding_angle(self, theta_0=np.pi/2):
        """
        Will find the angle $\phi$ that our picture (bildefil)
        is centered about (Can be in the interval [0,359] degrees)
        Compares the colorcode for the 16 pixels in the center of the image to the
        respective 16 pixels in the referance image and checks for equality.
        All RGB values must match!

        Returns the angle phi = i according to referance picture nr. i
        """
        pixels = np.array(self.img)
        center_pixels = pixels[241:245, 321:325]
        for i in range(360):
            ref_img = Image.open(f"Flatnr_{i}.png")
            ref_pixels = np.array(ref_img)
            center_ref_pixels = ref_pixels[241:245, 321:325]
            if np.all(center_pixels == center_ref_pixels):
                Angle_phi = i
                print(f"Angle that spaceshuttle was facing when {self.bildefil} was taken: phi = {Angle_phi} degrees")
                break
            else:
                continue

        return Angle_phi


if __name__ == "__main__":
    instance1 = BA("sample0000.png")
    angle = instance1.finding_angle()
    instance2 = BA("sample0200.png")
    angle = instance2.finding_angle()
    instance3 = BA("sample0435.png")
    angle = instance3.finding_angle()
    instance4 = BA("sample0911.png")
    angle = instance4.finding_angle()
    instance5 = BA("sample1400.png")
    angle = instance5.finding_angle()
    instance6 = BA("sample1900.png")
    angle = instance6.finding_angle()
