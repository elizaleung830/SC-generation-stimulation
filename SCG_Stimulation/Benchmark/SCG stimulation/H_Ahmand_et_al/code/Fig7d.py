'''
import shapely
from generate_neff_and_aeff import get_neff_and_aeff
from refractive_index import n_MgF2, n_Si3N4, n_Air
from collections import OrderedDict
import numpy as np

# waveguide parameters
width = 6 # um
height = 1 # um
thickness = 0.8 # um

n2 = 2.5e-19     # m^2/W n2 is the nonlinear refractive index at the center
Alpha = 0.7     # loss (dB/cm)

wavelength_range = [207,8000]
wavelegnth_step = 100 # 50nm steps

n_core = n_Si3N4
n_lower_cladding = n_MgF2
n_air = n_Air
#print(n_Si3N4(6))

# Construct waveguide geometry
core = shapely.geometry.box(-width/2, 0, +width/2, height)
lower_cladding = shapely.geometry.box(-width, -3, width, 0)
air_cladding = shapely.geometry.box(-width,0,width,width)
air = shapely.geometry.box(-10,-3,10,10)
polygons = OrderedDict(
    core = core,
    lower_cladding = lower_cladding,
    air_cladding = air_cladding,
    air = air
)

# Define material property and resolution of waveguide
resolutions = dict(core={"resolution": 0.04, "distance": 0.2},
                   air_cladding={"resolution": 0.15, "distance": 0.5},
                   lower_cladding={"resolution": 0.15, "distance": 0.5},
                   air = {"resolution": 0.5, "distance": 1} )

n_dict = {"core": n_core,"lower_cladding":n_lower_cladding ,"air_cladding": n_air, "air":n_air}

print("start")
# Calculate dispersion and gamma
aeff_list, neff_list, wls = get_neff_and_aeff(polygons,n_dict,wavelength_range,wavelegnth_step,resolutions)

##plot data
np.savez(f"data_h_{height}_w_{width}", wls=wls, aeff_list=aeff_list, neff_list=neff_list)

print("end")
print(aeff_list)
print(wls)
'''
import laserfun as lf
import matplotlib.pyplot as plt
import scipy
import numpy as np
import pandas as pd

# waveguide parameters
width = 6 # um
height = 1 # um

n2 = 2.5e-19     # m^2/W n2 is the nonlinear refractive index at the center
Alpha = 0.7     # loss (dB/cm)

# pulse parameters
FWHM = 50 * 1e-3  # pulse duration (ps)
pulseWL = 1550  # pulse central wavelength (nm)

power = 5000  # W
GDD = False   # Group delay dispersion (ps^2)
TOD = False  # Third order dispersion (ps^3)
FOD = False  # Fourth order dispersion

# simulation parameters
Window = 4  # simulation window (ps)
Steps = 100  # simulation steps
Points = 2 ** 14  # simulation points
rtol = 1e-4  # relative error for NLSE integrator
atol = 1e-4  # absolute error
Raman = True  # Enable Raman effect?
Steep = True  # Enable self steepening?

# ----------- Build Waveguide -----------
n2 = 2.5e-19  # m^2/W n2 is the nonlinear refractive index at the center
Alpha = 0.7  # loss (dB/cm)

data = np.load(f"data_h_{height}_w_{width}.npz")
wls = data['wls']
neff_list = data['neff_list']
aeff_list = data['aeff_list']
print(neff_list)
print(aeff_list)
print(wls)

def disp_function(z=0):  # provide effective index to the NLSE
    return (wls, neff_list)


# Calculate Gamma
def gamma_function(pump_wl):  # provide the nonlinearity at the pump to the NLSE
    aeff_interp = scipy.interpolate.interp1d(wls, aeff_list)
    gamma = 2 * np.pi * n2 / (pump_wl * 1e-9 * aeff_interp(pump_wl) * 1e-12)
    return gamma

# create the pulse:
p = lf.Pulse(pulse_type='sech', fwhm_ps=FWHM, center_wavelength_nm=pulseWL,
             time_window_ps=Window, power_is_avg=False, wav=wls, power=power,RI=neff_list, GDD=GDD,TOD=TOD, FOD=FOD, npts=Points)

# create the waveguide
f = lf.Fiber(10 * 1e-3, center_wl_nm=pulseWL, dispersion_format='GVD',
             gamma_W_m=gamma_function(pulseWL), loss_dB_per_m=Alpha * 100)

f.set_dispersion_function(disp_function, dispersion_format='n')

# propagate the pulse using the NLSE
results = lf.NLSE(p, f, raman=Raman, shock=Steep, nsaves=Steps, rtol=rtol,
                  atol=atol, print_status=True)
#Plot adjusting
ref_fig5d = pd.read_csv(
    "../reference_data/fig7d.csv", dtype=np.float64
)

fig, axes = results.plot(wavelength=True, show=False, tlim=(-2, 2), flim=(0, 8000), wmax = 8000, wmin = 200)
ref_fig5d_x, ref_fig5d_y = np.split(ref_fig5d.values, 2, axis=1)

axes[0][0].plot(ref_fig5d_x, ref_fig5d_y, c="green", label="ref spectrum", zorder = 0)

axes[1][0].set_xlim(0,8000)
plt.show()
print(f"Gamma is {gamma_function(1550)}/W/m at pump wavelength(1.55um)")