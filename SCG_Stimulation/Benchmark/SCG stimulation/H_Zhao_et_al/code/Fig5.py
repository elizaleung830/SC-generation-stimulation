import shapely
import scipy
from generate_neff_and_aeff import get_neff_and_aeff
from refractive_index import  n_Si3N4, n_SiO2
from collections import OrderedDict
import numpy as np
from femwell.visualization import plot_domains
import femwell
from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import compute_modes
from skfem import Basis, ElementTriP0
from skfem.io import from_meshio
import numpy as np
import matplotlib.pyplot as plt
from time import sleep

from tqdm import tqdm

def gamma_function(pump_wl):  # provide the nonlinearity at the pump to the NLSE
    aeff_interp = scipy.interpolate.interp1d(wls, aeff_list)
    gamma = 2 * np.pi * n2 / (pump_wl * 1e-9 * aeff_interp(pump_wl) * 1e-12)
    return gamma


# waveguide parameters
width = 0.88 # um
height = 0.69 # um
length = 7.5 *1e-3 # m

n2 = 2.4e-19     # m^2/W n2 is the nonlinear refractive index at the center
Alpha = 0     # loss (dB/cm)

wavelength_range = [210,2500]
wavelegnth_step = 50 # 50nm steps

n_core = n_Si3N4
n_cladding = n_SiO2
n_buried_oxide = n_SiO2
#print(n_Si3N4(6))

# Construct waveguide geometry
core = shapely.geometry.box(-width/2, 0, +width/2, height)
cladding = shapely.geometry.box(-width*2, 0, width*2, height*3)
buried_oxide = shapely.geometry.box(-width*2,-height*2,width*2,0)
polygon = OrderedDict(
    core = core,
    cladding = cladding,
    buried_oxide = buried_oxide,
)

# Define material property and resolution of waveguide
resolutions = dict(core={"resolution": 0.02, "distance": 0.2},
                   cladding={"resolution": 0.05, "distance": 0.2},
                   buried_oxide={"resolution": 0.05, "distance": 0.2} )

n_dict = {"core": n_core,"cladding":n_cladding ,"buried_oxide": n_buried_oxide}

print("start")
# Calculate dispersion and gamma
mesh = from_meshio(mesh_from_OrderedDict(polygon, resolutions, default_resolution_max=2))
mesh.draw().show()
plot_domains(mesh)
plt.show()

basis0 = Basis(mesh, ElementTriP0())
epsilon = basis0.zeros()
wavelength_list = np.linspace(wavelength_range[0], wavelength_range[1], wavelegnth_step)
neff_list = []
aeff_list = []
count = 0
"""
for wavelength in tqdm(wavelength_list):
    wavelength = wavelength * 1e-3
    for subdomain, n in n_dict.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n(wavelength) ** 2
    modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=3, order=1)
    modes_sorted = modes.sorted(key=lambda mode: -mode.calculate_power(elements="core").real)
    modes_sorted[0].show(modes_sorted[0].E.real, direction="x")
"""
wavelength = 1.055
for subdomain, n in n_dict.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n(wavelength) ** 2
modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=3, order=1)
modes_sorted = modes.sorted(key=lambda mode: -mode.calculate_power(elements="core").real)
modes_sorted[0].show(modes_sorted[0].E.real, direction="x")
aeff = modes_sorted[0].calculate_effective_area()
gamma = 2 * np.pi * n2 / (wavelength * 1e-6 * aeff * 1e-12)
print(f"Gamma is {gamma}/W/m at pump wavelength(1.055um), 3.25/W/m in the paper")

aeff_list = np.array(aeff_list)
neff_list = np.array(neff_list)
wls = np.array(wavelength_list)
##plot data
np.savez(f"data_h_{height}_w_{width}", wls=wls, aeff_list=aeff_list, neff_list=neff_list)

print("end")
print(aeff_list)
print(wls)
