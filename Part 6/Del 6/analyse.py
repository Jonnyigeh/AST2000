# EGEN KODE
import matplotlib.pyplot as plt
import numpy as np


class Analyse:
    def __init__(self):
        self.c = 299792458
        self.k = 1.38064852e-23
        self.v = 10000                                          # 10 km/s
        self.data = np.load('spectrum.npy')
        self.noise = np.load("sigma_noise.npy")
        self.u = 1.661e-27
        self.table = np.zeros((12, 4))

    def f(self, lambda_i, F_min, lambda_0, T, m):
        """
        Gir Gaussfordeling gitt lambda, lambda_0, temperatur og Fmin
            ved gaussisk linje profil formel
        """
        F = 1 + (1 - F_min) * np.exp(-0.5 * ((lambda_i - lambda_0) / self.sigma(lambda_0, m, T)) ** 2)
        return F

    def sigma(self, lambda_0, m, T):
        """
        Beregner sigma (bredde til gausskurve)
            gitt temperatur, masse og bølgelengde
                FWHM / sqrt(8)ln2
        """
        sigma = (lambda_0 / self.c * np.sqrt( (self.k * T) / m ))
        return sigma

    def find_and_plot(self, data, noise, m, l_0):
        """
        Utfører chimetoden for å finne
            beste parametre til Gaussisk linjefunksjon
                Plotter målingene med støy, med tilnærmet gaussfunksjon.

        Returnerer beste parametre
        """
        self.data_s = data
        self.data_s[:,0] *= 1e-9          # omgjør bølgelengde til enhet [m]
        self.noise_s = noise
        self.noise_s[:,0] *= 1e-9
        T = np.linspace(150, 450, 30)
        Fmin = np.linspace(0.6, 1, 8)
        lambda_val = np.linspace(self.data_s[:,0][0], self.data_s[:,0][-1], 40) #len(self.data[:,0]))
        minval = 100000
        for j in range(len(Fmin)):
            for k in range(len(T)):
                for i in range(len(lambda_val)):
                    chisquare = np.sum( ((self.data_s[:,1] - self.f(self.data_s[:,0], Fmin[j], lambda_val[i], T[k], m) ) / self.noise_s[:,1] ) ** 2 )
                    if chisquare < minval:
                        minval = chisquare
                        parameter = np.array([Fmin[j], lambda_val[i] * 1e9, T[k]])
                        index = np.array([j, i, k])
        self.data_s[:,0] *= 1e9               # Omgjør tilbake til [nm]
        self.noise_s[:,0] *= 1e9
        plt.plot(self.data_s[:,0], self.data_s[:,1])
        plt.plot(self.data_s[:,0], self.f(self.data_s[:,0], parameter[0], parameter[1], parameter[0], m))
        plt.legend(["Measured data", "Gaussian line profile"])
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Flux")
        plt.title(fr"Chisquare approximation for $\lambda_0 = {l_0}$")
        if l_0 == 720 or l_0 == 1600 or l_0 == 2200:
            plt.savefig(f"absorptionline_{l_0}.jpg")
        plt.show()
        return parameter

    def O2_632(self):
        """
        Verdier for forventet spektrallinje til O2 ved 632nm
            Henter ut sliced array med interval av bølgelengder i område rundt 632 nm
                +- dopplershift som regnes ut av maxhastighet 10km/s + termisk bevegelse
        Kaller så funksjonen som utfører chimetoden på dette intervallet, og lagrer det i
            tabellen som gir oss oversikt over spektrallinjer vi finner.

        Tilsvarende for alle de andre forventede abs.linjene
        """
        m = 2*15.999*self.u                                                                               # mass of molecule
        lambda_0 = 632                                                                                    # Wavelength expected spectralline
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c                        # maximum possible dopplershift with v = 10000 m/s + v_max of molecules in gases [nm]
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 3)) < tol)          # Lower index for interval [lambda_0 - deltaLambda, lambda_0 + deltaLambda]
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 4)) < tol)          # Upper index for interval [lambda_0 - deltaLambda, lambda_0 + deltaLambda]
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]                                      # Slice the array with all wavelengths so we just have the interval
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]                                    # likewise with the noisearray
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[0,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])

    def O2_690(self):
        m = 2*15.999*self.u
        lambda_0 = 690
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 4)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 4)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[1,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])

    def O2_760(self):
        m = 2*15.999*self.u
        lambda_0 = 760
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 3)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 3)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[2,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])


    def H2O_720(self):
        m = (2*1.00784 + 15.999) * self.u
        lambda_0 = 720
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 3)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 3)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[3,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])


    def H2O_820(self):
        m = (2*1.00784 + 15.999) * self.u
        lambda_0 = 820
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 3)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 4)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[4,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])


    def H2O_940(self):
        m = (2*1.00784 + 15.999) * self.u
        lambda_0 = 940
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 4)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 3)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[5,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])


    def CO2_1400(self):
        m = (12.0107 + 2*15.999) * self.u
        lambda_0 = 1400
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 6)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 4)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[6,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])


    def CO2_1600(self):
        m = (12.0107 + 2*15.999) * self.u
        lambda_0 = 1600
        delta_lambda = ((self.v + np.sqrt(self.k * 450 / m) ) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 4)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 3)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[7,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])


    def CH4_1660(self):
        m = (12.0107 + 4*1.00784) * self.u
        lambda_0 = 1660
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 4)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 3)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[8,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])


    def CH4_2200(self):
        m = (12.0107 + 4*1.00784) * self.u
        lambda_0 = 2200
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 3)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 4)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[9,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])


    def CO_2340(self):
        m = (12.0107 + 15.999) * self.u
        lambda_0 = 2340
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 4)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 4)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[10,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])

    def N2O_2870(self):
        m = (2*14.0067 + 15.999) * self.u
        lambda_0 = 2870
        delta_lambda = ( (self.v + np.sqrt(self.k * 450 / m)) * lambda_0) / self.c
        tol = 1e-4
        lower_index = np.where(abs(self.data[:,0] - np.round(lambda_0 - delta_lambda, 4)) < tol)
        upper_index = np.where(abs(self.data[:,0] - np.round(lambda_0 + delta_lambda, 4)) < tol)
        self.data_s = self.data[lower_index[0][0]:upper_index[0][0]]
        self.noise_s = self.noise[lower_index[0][0]:upper_index[0][0]]
        Fmin, lambda_, Temp = self.find_and_plot(self.data_s, self.noise_s, m, lambda_0)
        self.table[11,:] = np.array([lambda_0, Fmin, abs(lambda_0 - lambda_), Temp])


if __name__ == '__main__':
    instance = Analyse()
    instance.O2_632()
    instance.O2_690()
    instance.O2_760()

    instance.H2O_720()
    instance.H2O_820()
    instance.H2O_940()

    instance.CO2_1400()
    instance.CO2_1600()

    instance.CH4_1660()
    instance.CH4_2200()

    instance.CO_2340()

    instance.N2O_2870()

    print("Wavelength |  Fmin  | delta_Lambda | Temperature")
    np.set_printoptions(precision=4, suppress=True)
    print(instance.table)
