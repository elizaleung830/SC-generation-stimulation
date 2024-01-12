from femwell.maxwell.waveguide import compute_modes
from skfem import Basis, ElementTriP0
from tqdm import tqdm
from femwell.visualization import plot_domains
from skfem.io import from_meshio
import shapely
import matplotlib.pyplot as plt
from refractive_index import n_MgF2, n_Si3N4, n_Air
from collections import OrderedDict
from femwell.mesh import mesh_from_OrderedDict
import numpy as np

# waveguide parameters
width = 6  # um
height = 0.8  # um
thickness = 0.8  # um

n2 = 2.5e-19  # m^2/W n2 is the nonlinear refractive index at the center
Alpha = 0.7  # loss (dB/cm)

wavelength_range = [310, 5100]
wavelegnth_step = 70  #  steps

n_core = n_Si3N4
n_lower_cladding = lambda w: n_MgF2(w, ray="e") # TODO: check effect of this line
n_air = n_Air

# Construct waveguide geometry
core = shapely.geometry.box(-width / 2, 0, +width / 2, height)
lower_cladding = shapely.geometry.box(-6, -6, 6, 0)
air = shapely.geometry.box(-6, 0, 6, 6)
polygons = OrderedDict(
    core=core,
    lower_cladding=lower_cladding,
    air= air
)

# Define material property and resolution of waveguide
resolutions = dict(core={"resolution": 0.04, "distance": 0.1},
                   lower_cladding={"resolution": 0.15, "distance": 0.2},
                   air={"resolution": 0.2, "distance": 0.2})

n_dict = {"core": n_core, "lower_cladding": n_lower_cladding, "air": n_air}

# Calculate dispersion and gamma
mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions))
mesh.draw().show()
plot_domains(mesh)
plt.show()

print("start")
# Calculate dispersion and gamma
basis0 = Basis(mesh, ElementTriP0())
epsilon = basis0.zeros()
wavelength_list = np.linspace(wavelength_range[0], wavelength_range[1], wavelegnth_step)
neff_list = []
aeff_list = []
for wavelength in tqdm(wavelength_list):
    wavelength = wavelength * 1e-3
    for subdomain, n in n_dict.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n(wavelength) ** 2
    modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=3, order=1)
    modes_sorted = modes.sorted(key=lambda mode: -np.real(mode.n_eff))
    mode = modes_sorted[0]
    #mode.show(mode.E.real, direction ="x")
    neff_list.append(np.real(mode.n_eff))
    aeff_list.append(mode.calculate_effective_area())

neff_list = np.array(neff_list)
aeff_list = np.array(aeff_list)
wls = np.array(wavelength_list)

##save data
np.savez(f"data_h_{height}_w_{width}_ne", wls=wls, aeff_list=aeff_list, neff_list=neff_list)

print("end")
print(aeff_list)
print(neff_list)
print(wls)
