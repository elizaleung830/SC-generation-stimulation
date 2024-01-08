from femwell.mesh import mesh_from_OrderedDict
from skfem.io import from_meshio
from femwell.visualization import plot_domains
from typing import OrderedDict
import shapely
import pandas as pd
import laserfun as lf
from generate_neff_and_aeff import get_neff_and_aeff
from refractive_index import n_MgF2, n_Si3N4, n_Air
from collections import OrderedDict
import scipy
import numpy as np
import matplotlib.pyplot as plt

# waveguide parameters
width = 6  # um
height = 0.8  # um
thickness = 0.8  # um

n2 = 2.5e-19  # m^2/W n2 is the nonlinear refractive index at the center
Alpha = 0.7  # loss (dB/cm)

wavelength_range = [310, 5500]
wavelegnth_step = 100  # 50nm steps

n_core = n_Si3N4
n_lower_cladding = n_MgF2
n_air = n_Air

# Construct waveguide geometry
core = shapely.geometry.box(-width / 2, 0, +width / 2, height)
lower_cladding = shapely.geometry.box(-6, -3, 6, 0)
air_cladding = shapely.geometry.box(-6, 0, 6, 6)
air = shapely.geometry.box(-8, -3, 8, 8)
polygons = OrderedDict(
    core=core,
    lower_cladding=lower_cladding,
    air_cladding=air_cladding,
    air= air
)

# Define material property and resolution of waveguide
resolutions = dict(core={"resolution": 0.04, "distance": 0.2},
                   air_cladding={"resolution": 0.15, "distance": 0.5},
                   lower_cladding={"resolution": 0.15, "distance": 0.5},
                   air={"resolution": 0.5, "distance": 1})

n_dict = {"core": n_core, "lower_cladding": n_lower_cladding, "air_cladding": n_air, "air": n_air}

print("start")
# Calculate dispersion and gamma
aeff_list, neff_list, wls = get_neff_and_aeff(polygons, n_dict, wavelength_range, wavelegnth_step, resolutions)

##save data
np.savez(f"data_h_{height}_w_{width}", wls=wls, aeff_list=aeff_list, neff_list=neff_list)

print("end")
print(aeff_list)
print(wls)
