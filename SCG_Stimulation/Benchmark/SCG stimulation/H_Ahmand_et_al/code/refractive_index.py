import math
import pandas as pd
import scipy
import numpy as np
import matplotlib.pyplot as plt
"""
From refractiveIndex.info
"""
def n_Si3N4(wavelength):
    """
    stiometric SiN, data interplot from https://opg.optica.org/ol/fulltext.cfm?uri=ol-40-21-4823&id=331311
    :param wavelength: in um
    :return: linear refractive index of Si3N4/ SiN
    """
    if wavelength >=0.31 and wavelength <= 5.507:
        return math.sqrt((3.0249 * wavelength ** 2) / (wavelength ** 2 - 0.1353406 ** 2) + (40314 * wavelength ** 2) / (
                wavelength ** 2 - 1239.842 ** 2) + 1)
    elif wavelength <= 0.31 and wavelength >= 0.207:
        return math.sqrt(1+ ((2.8939 * wavelength**2)/(wavelength**2- 0.13967**2)))
    elif wavelength >= 5.507:
        ## Not as accurate as the above
        """
        ref_data = pd.read_csv("./n_si3n4.csv", dtype=np.float64)
        ref_x, ref_y = np.split(ref_data.values, 2, axis=1)
        ref_x = list(map(lambda v: v[0],ref_x))
        ref_y = list(map(lambda v: v[0], ref_y))
        """
        ref_data = np.load("n_Si3N4.npz")
        n_interplolate = scipy.interpolate.interp1d(ref_data['ref_x'], ref_data['ref_y'])
        return n_interplolate(wavelength)
    else:
        raise ValueError(f"wavelength provided is {wavelength}um, is out of the range for Si3N4")

def n_MgF2(wavelength):
    """
    valid for  (0.2–7)um
    :param wavelength: in um
    :return: linear refractive index of MgF2
    """
    return math.sqrt(0.48755108*wavelength**2/(wavelength**2 - 0.04338408**2)+0.39875031*wavelength**2/(wavelength**2-0.09461442**2)+2.3120353*wavelength**2/(wavelength**2 - 23.793604**2)+1)

def n_Air(wavvelength):
    return 0.05792105/(238.0185-wavvelength**(-2))+0.00167917/(57.362-wavvelength**(-2)) + 1

