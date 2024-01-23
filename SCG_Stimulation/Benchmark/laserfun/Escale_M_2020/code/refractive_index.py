import math
import pandas as pd
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.constants import pi, c
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


def n_LNOI(wavelength, ray="o"):
    '''
    valid for 0.4 - 5 um, n(o) of Lithium niobate
    :param wavelength:
    :return:
    '''
    if ray == "o":
        if wavelength >= 0.4 and wavelength <=5:
            return math.sqrt(2.6734*wavelength**2/(wavelength**2 - 0.01764) + 1.2290*wavelength**2/(wavelength**2 - 0.05914) + 12.614*wavelength**2/(wavelength**2 - 474.6)+1)
        else:
            raise ValueError("invalid wavelength for lithium niobate, must be between 0.4-5um")
    elif ray == "e":
        if wavelength >= 0.4 and wavelength <=5:
            return math.sqrt(2.9804*wavelength**2/(wavelength**2-0.02047)+ 0.5981*wavelength**2/(wavelength**2 -0.0666) + 8.9543 * wavelength**2 / (wavelength**2 - 416.08)+1)
        else:
            raise ValueError("invalid wavelength for lithium niobate, must be between 0.4-5um")

def n_LNOI_1(wls,doped="undoped",  T=24.5):
    """
    Refractive index of congruent lithium niobate.
    :param wls: in um
    References
    ----------
    Dieter H. Jundt, "Temperature-dependent Sellmeier equation for the index of
     refraction, ne, in congruent lithium niobate," Opt. Lett. 22, 1553-1555
     (1997). https://doi.org/10.1364/OL.22.001553

    """
    if doped == "undoped":
        # Undoped
        a1 = 5.35583
        a2 = 0.100473
        a3 = 0.20692
        a4 = 100.
        a5 = 11.34927
        a6 = 1.5334e-2
        b1 = 4.629e-7
        b2 = 3.862e-8
        b3 = -0.89e-8
        b4 = 2.657e-5
    elif doped == "5%":
        # Doped 5% MgO:LN
        a1 = 5.756
        a2 = 0.0983
        a3 = 0.2020
        a4 = 189.32
        a5 = 12.52
        a6 = 1.32e-2
        b1 = 2.86e-6
        b2 = 4.7e-8
        b3 = 6.113e-8
        b4 = 1.516e-4

    f = (T-24.5)*(T+570.82)
    n2 = (a1 + b1*f + (a2 + b2*f)/(wls**2 - (a3 + b3*f)**2)
          + (a4 + b4*f)/(wls**2 - a5**2) - a6*wls**2)
    return n2**0.5

def n_SiO2(wavelength,type="FusedSilica"):
    if type == "FusedSilica":
        if wavelength < 0.21 or wavelength > 6.7:
            raise ValueError(f"wavelength provided is {wavelength}um, is out of the range for {type}")

        return np.sqrt( 0.6961663* wavelength**2/(wavelength**2 - 0.0684043**2)+(0.4079426*wavelength**2/(wavelength**2-0.1162414**2))+(0.8974794*wavelength**2/(wavelength**2-9.896161**2))+1)

    elif type == "flim":
        return 1.45


def n_Air(wavvelength):
    return 0.05792105/(238.0185-wavvelength**(-2))+0.00167917/(57.362-wavvelength**(-2)) + 1


