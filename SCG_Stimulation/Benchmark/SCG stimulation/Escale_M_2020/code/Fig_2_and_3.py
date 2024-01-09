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
wavelegnth_step = 50  # 50nm steps

n_core = n_LNOI
n_ridge = n_LNOI
n_buffer = n_SiO2
n_air = n_Air

# waveguide parameters
width = 0.7  # um
total_height = 0.4  # um
ridge_height = 0.25
box_height = 3

triangle_height = 0.25
triangle_width = triangle_height/ math.tan(60 * math.pi /180)
left_triangle = [(-width/2, 0), (-width/2 - triangle_width, 0), (-width/2, triangle_height)]
right_triangle = [(+width/2, 0), (width/2 + triangle_width, 0), (+width/2, triangle_height)]


# Construct waveguide geometry
core_trapiz = Polygon([(width/2 + triangle_width,0 ),(-width/2 - triangle_width, 0),(-width/2, triangle_height) , (+width/2, triangle_height)])
core_box = shapely.geometry.box(-(width/2 + triangle_width), -0.15, (width/2 + triangle_width), 0)
core = unary_union([core_trapiz, core_box])

ridge = shapely.geometry.box(-box_height, -0.15,box_height , 0)

buffer = shapely.geometry.box(-box_height,-box_height/2,box_height,-0.15)
air = shapely.geometry.box(-box_height,-box_height/2,box_height,box_height*1.5)


polygon = OrderedDict(
    core = core,
    ridge = ridge,
    buffer = buffer,
    air= air
)

# Define material property and resolution of waveguide
resolutions = dict(core={"resolution": 0.01, "distance": 0.1},
                   ridge ={"resolution": 0.01, "distance": 0.1},
                   buffer={"resolution": 0.1, "distance": 0.5},
                   air={"resolution": 0.1, "distance": 0.5})

n_dict = {"core": n_core, "ridge": n_ridge, "buffer": n_buffer, "air": n_air}

mesh = from_meshio(mesh_from_OrderedDict(polygon, resolutions))
mesh.draw().show()
plot_domains(mesh)
plt.show()

#----------------------FEM solver-------------------------------
print("start")
# Calculate dispersion and gamma
mesh = from_meshio(mesh_from_OrderedDict(polygon, resolutions))
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
    neff_list.append(np.real(mode.n_eff))
    aeff_list.append(mode.calculate_effective_area())

neff_list = np.array(aeff_list)
aeff_list = np.array(neff_list)
wls = np.array(wavelength_list)

##save data
np.savez(f"data_h_{ridge_height}_w_{width}", wls=wls, aeff_list=aeff_list, neff_list=neff_list)

print("end")
print(aeff_list)
print(wls)
