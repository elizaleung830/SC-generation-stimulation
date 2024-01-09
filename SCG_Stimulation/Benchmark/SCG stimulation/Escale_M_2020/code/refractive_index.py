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
    stiometric SiN, data fitted from https://opg.optica.org/ol/fulltext.cfm?uri=ol-40-21-4823&id=331311
    might not be valid for value outside 0.31-5.507um if using fit.
    :param wavelength: in um
    :param fit: determine if the function use sellmeier equation or by using data interplot
    :return: linear refractive index of Si3N4/ SiN
    """
    if fit == False:
        if wavelength >= 0.31 and wavelength <= 5.507:
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


def n_LNOI(wavelength):
    '''
    valid for 0.4 - 5 um, n(o) of Lithium niobate
    :param wavelength:
    :return:
    '''
    if wavelength >= 0.4 and wavelength <=5:
        return math.sqrt(2.6734*wavelength**2/(wavelength**2 - 0.01764) + 1.2290*wavelength**2/(wavelength**2 - 0.05914) + 12.614*wavelength**2/(wavelength**2 - 474.6)+1)
    else:
        raise ValueError("invalid wavelength for lithium niobate, must be between 0.4-5um")
def n_MgF2(wavelength):
    """
    valid for  (0.2–7)um
    :param wavelength: in um
    :return: linear refractive index of MgF2
    """
    return math.sqrt(0.48755108*wavelength**2/(wavelength**2 - 0.04338408**2)+0.39875031*wavelength**2/(wavelength**2-0.09461442**2)+2.3120353*wavelength**2/(wavelength**2 - 23.793604**2)+1)

def n_SiO2(wavelength,type="FusedSilica"):
    if type == "FusedSilica":
        if wavelength < 0.21 or wavelength > 6.7:
            raise ValueError(f"wavelength provided is {wavelength}um, is out of the range for {type}")

        return np.sqrt( 0.6961663* wavelength**2/(wavelength**2 - 0.0684043**2)+(0.4079426*wavelength**2/(wavelength**2-0.1162414**2))+(0.8974794*wavelength**2/(wavelength**2-9.896161**2))+1)

    if type == "flim":
        return None


def n_Air(wavvelength):
    return 0.05792105/(238.0185-wavvelength**(-2))+0.00167917/(57.362-wavvelength**(-2)) + 1

"""
x = [x * 1e-3 for x in range(400, 5000, 10)]
y = [n_LNOI(i) for i in x]
plt.plot(x,y)
plt.show()
"""

