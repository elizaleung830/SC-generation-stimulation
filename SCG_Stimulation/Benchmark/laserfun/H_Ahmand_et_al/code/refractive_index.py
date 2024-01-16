import csv
import math
import pandas as pd
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

"""
From refractiveIndex.info
"""


def n_Si3N4(wavelength, fit = False):
    """
    stiometric SiN, data interplot from https://opg.optica.org/ol/fulltext.cfm?uri=ol-40-21-4823&id=331311
    might not be valid for value outside 0.31-5.507um
    :param wavelength: in um
    :param fit: determine if the function use sellmeier equation or equation from fitting data
    :return: linear refractive index of Si3N4/ SiN
    """
    if fit == False:
        #if wavelength >= 0.31 and wavelength <= 5.507:
            return math.sqrt(
                (3.0249 * wavelength ** 2) / (wavelength ** 2 - 0.1353406 ** 2) + (40314 * wavelength ** 2) / (
                        wavelength ** 2 - 1239.842 ** 2) + 1)
        #else:
            #raise ValueError(f"wavelength provided is {wavelength}um, is out of the range for Si3N4")
    elif fit == True:
        if wavelength >= 0.31 and wavelength <= 5.507:
            n_x, n_y = list(np.split(pd.read_csv(
                "../reference_data/n_Si3N4.csv", dtype=np.float64
            ).values, 2, axis=1))
            n_x = np.array(list(map(lambda v: v[0], n_x)))
            n_y = np.array(list(map(lambda v: v[0], n_y)))

            y_spl = UnivariateSpline(n_x, n_y, s=0.0012, k=3)
            return y_spl(wavelength)

        else:
            raise ValueError(f"wavelength provided is {wavelength}um, is out of the range for Si3N4")




def n_MgF2(wavelength, ray="o"):
    """
    valid for  (0.2–7)um
    :param wavelength: in um
    :return: linear refractive index of MgF2
    """
    if ray == "o":
        return math.sqrt(0.48755108*wavelength**2/(wavelength**2 - 0.04338408**2) + 0.39875031 * wavelength**2 / (wavelength**2 - 0.09461442**2) + 2.3120353 * wavelength**2 /(wavelength**2 - 23.793604**2) + 1)
    elif ray == "e":
        return math.sqrt(1+ 0.41344023*wavelength**2 /(wavelength**2 - 0.03684262**2) + 0.50497499*wavelength**2/(wavelength**2-0.09076162**2)+2.4904862*wavelength**2/(wavelength**2-23.771995**2))
    elif ray == "film":
        return 1.422

def n_Air(wavvelength):
    return 1.0003

'''
n_x, n_y = list(np.split(pd.read_csv(
    "../reference_data/n_Si3N4.csv", dtype=np.float64
).values, 2, axis=1))
n_x = np.array(list(map(lambda v: v[0]*1000, n_x)))
n_y = np.array(list(map(lambda v: v[0], n_y)))

y_spl = UnivariateSpline(n_x, n_y, s = 0.0012, k =3)

x = [x for x in range(310,8000,10)]

plt.plot(x,y_spl(x), color = "red")
plt.show()

# Code used to combine the sellmeir equation and data plot
with open("n_Si3N4.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    for x in range(310,5070):
        x = round(x*1e-3,4)
        y = n_Si3N4(x)
        writer.writerow([str(x), str(y)])
'''
