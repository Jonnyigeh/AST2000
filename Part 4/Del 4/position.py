#EGEN KODE
import numpy as np
from ast2000tools.shortcuts import SpaceMissionShortcuts
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem


class Position:
    def __init__(self, time, distances):
        """
        Konstruktør tar paramtrene
        - Time: tid ved triangulering
        - distances: Avstand fra Spacemission.measure_distance()
        """
        self.N = 100000000
        self.time = np.array([time]) #note! must give as an one dimensional array, siden snarveien tar bare arrays.
        self.mission = SpaceMission(76408)
        self.shortcuts = SpaceMissionShortcuts(self.mission, [1925])
        self.object_position = self.shortcuts.compute_planet_trajectories(self.time)
        self.x1 = np.array([self.object_position[0][0], self.object_position[1][0]]) #posisjonen til planet med index 0 i tid "time"
        self.x2 = np.array([self.object_position[0][1], self.object_position[1][1]]) #posisjonen til planet med index 11 i tid "time"
        self.dist_s = distances[-1]                                                  #Distansen fra sonden til sola
        self.dist_1 = distances[0]                                                   #Distansen fra sonden til planet med index 0
        self.dist_2 = distances[1]                                                   #Distansen fra sonden til planet med index 1


    def find_position(self):
        """
        Finner posisjon til spaceshuttle ved triangulering
        og returnerer tuppel (x, y):
        Printer også posisjon til terminal.

        definerer nye x1/2 arrays som fixed_x1/2. Grunnet til det er at når man tar r_vec - x1, får man en array med 4 elementer.
        Dette blir feil når vi finne absoluttverdi av denne vektoren.
        """
        self.fixed_x1 = np.array([np.linalg.norm(self.x1[0]), np.linalg.norm(self.x1[1])])
        self.fixed_x2 = np.array([np.linalg.norm(self.x2[0]), np.linalg.norm(self.x2[1])])

        theta = 0
        for i in np.linspace(0,2*np.pi,self.N):

            self.r_s = self.dist_s
            self.r_vec = np.array([np.cos(i), np.sin(i)]) * self.r_s
            tol = 0.000001
            if (np.linalg.norm(self.r_vec - self.fixed_x1) - self.dist_1) < tol and (np.linalg.norm(self.r_vec - self.fixed_x2) - self.dist_2) < tol:
                theta = i
                print("PING! You've got an attractive candidate for theta! ;)")
                break

        print(f"the angle is: {theta}")
        x = self.r_s * np.cos(theta)
        y = self.r_s * np.sin(theta)
        print(f"your spacecrafts position is ({x},{y}) relative to the sun")


        return x, y

if __name__ == "__main__":
    distances = mission.measure_distances()
    instance = Position(final_time, distances)
    instance.find_position()
