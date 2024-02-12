import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import griddata

ref = pd.read_csv("AC_deformation.csv", dtype=np.float64)
ref = np.array(ref.values.tolist())
e = ref[:,2]
x = ref[:,0]
y = ref[:,1]
ref = [[i[0],i[1]]  for i in ref]
extent = (min(x),max(x),min(y),max(y))

grid_x, grid_y = np.mgrid[min(x):max(x):200j, min(y):max(y):200j]
grid = griddata(ref, e, (grid_x,grid_y), method='cubic')

plt.imshow(grid.T, extent=extent, cmap="jet", origin='lower')
plt.colorbar()
plt.show()
