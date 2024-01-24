
import math

from skfem import Basis, ElementTriP0
from tqdm import tqdm
from femwell.mesh import mesh_from_OrderedDict
from skfem.io import from_meshio
from femwell.visualization import plot_domains
import shapely
from femwell.maxwell.waveguide import compute_modes
from refractive_index import n_SiO2, n_Air, n_LNOI
from collections import OrderedDict
import scipy
import numpy as np
from matplotlib import pyplot as plt
from shapely.ops import unary_union
from shapely.geometry import Polygon

wavelength_range = [400, 1500]
wavelegnth_step = 70  # steps

n_core = lambda w: n_LNOI(w, ray="o") # TODO: Check effect of this line
n_ridge = n_core
n_buffer = lambda w: n_SiO2(w)# TODO: Check effect of this line
n_air = n_Air

# waveguide parameters
width = 0.7  # um
total_height = 0.4  # um
ridge_height = 0.25
box_height = 4

triangle_height = 0.25
triangle_width = triangle_height/ math.tan(60 * math.pi /180)


# Construct waveguide geometry
core_trapiz = Polygon([(width/2 + triangle_width,0 ),(-width/2 - triangle_width, 0),(-width/2, triangle_height) , (+width/2, triangle_height)])
core_box = shapely.geometry.box(-(width/2 + triangle_width), -0.15, (width/2 + triangle_width), 0)
core = unary_union([core_trapiz, core_box])
ridge = shapely.geometry.box(-box_height/2, -0.15,box_height/2 , 0)
core = unary_union([core_trapiz, ridge])


buffer = shapely.geometry.box(-box_height/2,-box_height/2,box_height/2,-0.15)
air = shapely.geometry.box(-box_height/2,-box_height/2,box_height/2,box_height/2)


polygon = OrderedDict(
    core = core,
    buffer = buffer,
    air= air
)

# Define material property and resolution of waveguide
resolutions = dict(core={"resolution": 0.015, "distance": 0.1},
                   buffer={"resolution": 0.1, "distance": 0.5},
                   air={"resolution": 0.1, "distance": 0.5})

n_dict = {"core": n_core, "buffer": n_buffer, "air": n_air}

mesh = from_meshio(mesh_from_OrderedDict(polygon, resolutions, filename = "mesh.msh"))
mesh.draw().show()
plot_domains(mesh)
plt.show()
print(mesh)


#----------------------FEM solver-------------------------------
print("start")
# Calculate dispersion and gamma
basis0 = Basis(mesh, ElementTriP0())
epsilon = basis0.zeros()
wavelength_list = np.linspace(wavelength_range[0], wavelength_range[1], wavelegnth_step)
neff_list_te = []
aeff_list_te = []
neff_list_tm = []
aeff_list_tm = []

for wavelength in tqdm(wavelength_list):
    wavelength = wavelength * 1e-3
    for subdomain, n in n_dict.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n(wavelength) ** 2
    modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=3, order=1)
    ## te mode
    modes_sorted = modes.sorted(key=lambda mode: -np.real(mode.te_fraction))
    mode = modes_sorted[0]
    neff_list_te.append(np.real(mode.n_eff))
    aeff_list_te.append(mode.calculate_effective_area())

    ## tm mode
    modes_sorted = modes.sorted(key=lambda mode: -np.real(mode.tm_fraction))
    if modes_sorted[0].tm_fraction < 0.7:
        print(f"at {wavelength}um, mode has highest tm_fraction of f{modes_sorted[0].tm_fraction}")
    mode = modes_sorted[0]
    neff_list_tm.append(np.real(mode.n_eff))
    aeff_list_tm.append(mode.calculate_effective_area())

neff_list_te = np.array(neff_list_te)
aeff_list_te = np.array(aeff_list_te)

neff_list_tm = np.array(neff_list_tm)
aeff_list_tm = np.array(aeff_list_tm)
wls = np.array(wavelength_list)

##save data
np.savez(f"data_h_{ridge_height}_w_{width}_no", wls=wls, aeff_list_te=aeff_list_te, neff_list_te=neff_list_te, neff_list_tm=neff_list_tm,aeff_list_tm=aeff_list_tm)

print("end")

print(wls)


