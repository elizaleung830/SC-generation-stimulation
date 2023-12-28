# https://iopscience.iop.org/article/10.1088/1555-6611/aaf63d
# H Ahmad et al 2019 Laser Phys. 29 025301
import laserfun as lf
import scipy
import numpy as np

# pulse parameters
FWHM = 0.05       # pulse duration (ps)
pulseWL = 1550        # pulse central wavelength (nm)
peak_power = 5000 # (W)
pulse_type = "sech"

GDD = True            # Group delay dispersion (ps^2)
TOD = True           # Third order dispersion (ps^3)
FOD = True           # Fourth order dispersion

# simulation parameters
Window = 3.0    # simulation window (ps)
Steps = 100     # simulation steps
Points = 2**12  # simulation points
rtol = 1e-4     # relative error for NLSE integrator
atol = 1e-4     # absolute error
Raman = True    # Enable Raman effect?
Steep = False   # Enable self steepening?

# ----------- Build Pulse -----------
pulse = lf.Pulse(pulse_type=pulse_type, center_wavelength_nm=pulseWL,
                 fwhm_ps=FWHM, power=peak_power, power_is_avg=False)

# ----------- Build Waveguide -----------
# waveguide parameters
width = 6 # um
height = 0.8 # um
thickness = 0.8 # um

def disp_function(z=0):  # provide effective index to the NLSE
    return (wls, neff)

def gamma_function(z=0):  # provide the nonlinearity at the pump to the NLSE
    aeff_interp = scipy.interpolate.interp1d(wls, aeff)
    return 2*np.pi*n2/(pulseWL*1e-9*aeff_interp(pulseWL)*1e-12)

n2 = 2.0e-19     # m^2/W n2 is the nonlinear refractive index at the center
Alpha = 0.0      # loss (dB/cm)

wavelength_range = [500,5000]
wavelegnth_step = 50 # 50nm steps

waveguide = lf.Fiber(thickness, center_wl_nm=pulseWL, dispersion_format='GVD',
             gamma_W_m= gamma_function(), loss_dB_per_m=Alpha*100)

waveguide.set_dispersion_function(disp_function, dispersion_format='n')
# plotting
