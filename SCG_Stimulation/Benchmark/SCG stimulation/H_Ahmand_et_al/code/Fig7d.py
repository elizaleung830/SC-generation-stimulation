import shapely
from generate_neff_and_aeff import get_neff_and_aeff
from refractive_index import n_MgF2, n_Si3N4, n_Air
from collections import OrderedDict
import numpy as np

# waveguide parameters
width = 6  # um
height = 1  # um
thickness = 0.8  # um

n2 = 2.5e-19  # m^2/W n2 is the nonlinear refractive index at the center
Alpha = 0.7  # loss (dB/cm)

wavelength_range = [207, 8000]
wavelegnth_step = 100  # steps

n_core = n_Si3N4
n_lower_cladding = n_MgF2
n_air = n_Air

# Construct waveguide geometry
core = shapely.geometry.box(-width / 2, 0, +width / 2, height)
lower_cladding = shapely.geometry.box(-width, -3, width, 0)
air_cladding = shapely.geometry.box(-width, 0, width, width)
air = shapely.geometry.box(-10, -3, 10, 10)
polygons = OrderedDict(
    core=core,
    lower_cladding=lower_cladding,
    air_cladding=air_cladding,
    air=air
)

# Define material property and resolution of waveguide
resolutions = dict(core={"resolution": 0.04, "distance": 0.1},
                   air_cladding={"resolution": 0.15, "distance": 0.2},
                   lower_cladding={"resolution": 0.15, "distance": 0.2},
                   air={"resolution": 0.5, "distance": 0.5})

n_dict = {"core": n_core, "lower_cladding": n_lower_cladding, "air_cladding": n_air, "air": n_air}

print("start")

# Calculate dispersion and gamma
aeff_list, neff_list, wls = get_neff_and_aeff(polygons, n_dict, wavelength_range, wavelegnth_step, resolutions)

##plot data
np.savez(f"data_h_{height}_w_{width}", wls=wls, aeff_list=aeff_list, neff_list=neff_list)

print("end")
print(aeff_list)
print(wls)
