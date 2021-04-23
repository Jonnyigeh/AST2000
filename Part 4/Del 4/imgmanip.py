# EGEN KODE
import numpy as np
from PIL import Image
from ast2000tools.space_mission import SpaceMission
"""
Dette er kodesnippet fra oppgaveteksten

img = Image.open("example.png")
pixels = np.array(img)
width = len(pixels[0, :])
redpixs = [(255,0,0) for ii in range(width)]
pixels[500,:] = redpixs
img2 = Image.fromarray(pixels)
img2.save("Examplewithredline.png")
"""

class BM:                       # BM = Bildemanipulering
    """
    Initialiseres ved å lage instanse med filnavn i en streng "example.png"
    og field of view FOV i grader. FOV satt til default 70 grader.
    Trenger også phi_0 og theta_0 (Altså hvor bildet er sentrert rundt)
    """
    def __init__(self, bildefil, FOV=70):
        self.bildefil = bildefil
        self.img = Image.open(self.bildefil)
        self.FOV = FOV
    def find_size(self):
        """
        Returns size of matrix of pixelgrid,
        width of pixelgrid i.e
        height of pixelgrid i.e
        """
        pixels = np.array(self.img)
        self._size = np.size(pixels)
        self._width = len(pixels[0, :])
        self._height = len(pixels[:, 0])

    def stereographic_proj(self, phi_0, theta_0):
        """
        Computes the stereographic projection using
        the formulas given in TASK 3 part 2 coordinate transformation

        Saves the solutions to array theta and phi.
        """
        self.theta_0 = theta_0
        self.phi_0 = phi_0
        alpha = self.FOV * np.pi / 180
        Xmax = 2 * np.sin(alpha/2) / (1 + np.cos(alpha / 2))
        Xmin = - Xmax
        Ymax = Xmin                                             # Flippet Y-array for å unngå oppned speiling av bildet
        Ymin = Xmax
        self.x = np.linspace(Xmin, Xmax, self._width)
        self.y = np.linspace(Ymin, Ymax, self._height)
        X, Y = np.meshgrid(self.x, self.y)
        rho = np.sqrt(X ** 2 + Y ** 2)
        beta = 2 * np.arctan(rho / 2)
        with np.errstate(divide="ignore", invalid="ignore"):    # Ignores division by zero in origin
            self.theta = self.theta_0 - np.arcsin(np.cos(beta) * np.cos(self.theta_0) + Y / rho * np.sin(beta) * np.sin(self.theta_0))                  # 3a
            self.phi = self.phi_0 + np.arctan(X * np.sin(beta) / (rho * np.sin(self.theta_0) * np.cos(beta) - Y * np.cos(self.theta_0) * np.sin(beta))) # 3b
        self.theta[240][320] = self.theta_0    # Found where the denomitor in the two equations equals zero and set that index to their intial value (Should only be at origin)
        self.phi[240][320] = self.phi_0        # Used np.where(denomiator == 0) and saw that phi/theta at this index == "nan"
        self.phi[242][320] = self.phi_0

    def generate_png(self):
        """
        Converts theta and phi values into sky_image_pixels using spacemission function
        and from pixels to RGB values using himmelkule.npy array
        Lastly it will create the PNG-file with these colors
        and return this PIL.image file.
        """
        himmelkule = np.load("himmelkule.npy")
        pixel_index = np.zeros((480, 640))
        pixels = np.zeros((480, 640, 3), dtype="uint8")
        for i in range(480):
            for j in range(640):
                pixel_index[i][j] = SpaceMission.get_sky_image_pixel(self.theta[i][j], self.phi[i][j])
                pixels[i][j][0], pixels[i][j][1], pixels[i][j][2] = himmelkule[int(pixel_index[i][j])][2], himmelkule[int(pixel_index[i][j])][3], himmelkule[int(pixel_index[i][j])][4]
        skypic = Image.fromarray(pixels)
        return skypic

    """
    Property funksjoner for å returnere verdier
    knytter til png filen
    Størrelse (P_h * P_b)
    Høyde (antall piksler)
    Bredde (antall piksler)
    """
    @property
    def size(self):
        return self._size
    @property
    def width(self):
        return self._width
    @property
    def height(self):
        return self._height


if __name__ == "__main__":
    for i in range(360):
        """
        Produserer de 360 referansebildene fra 2.3
        """
        instance = BM("sample0000.png", 70)
        instance.find_size()
        instance.stereographic_proj(i*np.pi/180, np.pi/2)
        skypic = instance.generate_png()
        skypic.save(f"Flatnr_{i}.png")
