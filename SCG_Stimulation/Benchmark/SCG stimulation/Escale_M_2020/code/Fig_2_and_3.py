import math
from femwell.mesh import mesh_from_OrderedDict
from skfem.io import from_meshio
from femwell.visualization import plot_domains
import shapely
from shapely.ops import clip_by_rect
import pandas as pd
import laserfun as lf
from refractive_index import n_SiO2, n_Air, n_LNOI
from collections import OrderedDict
import scipy
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from matplotlib import pyplot as plt
import geopandas as gpd
from matplotlib import pyplot as plt
from shapely.ops import unary_union
from shapely.geometry import Polygon

n2 = 2.5e-19  # m^2/W n2 is the nonlinear refractive index at the center

wavelength_range = [310, 1500]
wavelegnth_step = 20  # 50nm steps

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
resolutions = dict(core={"resolution": 0.05, "distance": 0.1},
                   ridge_right={"resolution": 0.1, "distance": 0.1},
                   ridge_left={"resolution": 0.1, "distance": 0.1},
                   buffer={"resolution": 0.15, "distance": 0.5},
                   air={"resolution": 0.15, "distance": 0.5})

n_dict = {"core": n_core, "ridge": n_ridge, "buffer": n_buffer, "air": n_air}

mesh = from_meshio(mesh_from_OrderedDict(polygon, resolutions))
mesh.draw().show()
plot_domains(mesh)
plt.show()


print("start")
# Calculate dispersion and gamma
aeff_list, neff_list, wls = get_neff_and_aeff(polygons, n_dict, wavelength_range, wavelegnth_step, resolutions)

##save data
np.savez(f"data_h_{height}_w_{width}", wls=wls, aeff_list=aeff_list, neff_list=neff_list)

print("end")
print(aeff_list)
print(wls)
