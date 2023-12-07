"""Tools for working with laser pulses."""

import numpy as np
import sympy as sym
import scipy.fftpack as fft


# speed of light in m/s and nm/ps
c_mks = 299792458.0
c_nmps = c_mks * 1e9/1e12


class Pulse:
    """Generate a new pulse object based on a pre-defined pulse shapes.

    Customized pulses can be generated by calling this function and then
    manually setting the time- or frequency domain pulse-profile using
    ``pulse.at`` or ``pulse.aw``. 

    Note that the pulse object is intrinsically a single laser pulse, and the
    time and frequency domain electric fields refer to the single pulse, not
    the average power of the laser at some repetition rate. The `
    `power_is_avg`` and ``frep_MHz`` values are used only to set the electric
    field during the generation of the pulse and are not stored in the pulse
    object. 

    Parameters
    ----------
    pulse_type : string
        The shape of the pulse in the time domain. Options are:

        - sech, which produces a hyperbolic secant (sech) shaped pulse
          ``A(t) = sqrt(power) * sech(t/T0)``, where ``T0=fwhm/1.76``
        - gaussian, which produces a Gaussian shaped pulse
          ``A(t) = sqrt(power) * exp(-(t/T0)^2/2)``, where ``T0=fwhm/1.76``
        - sinc, which uses a sin(x)/x (sinc) functionT
          ``A(t) = sqrt(power) * sin(t/T0)/(t/T0)``, were ``T0=fwhm/3.79``

    center_wavelength_nm : float
        The center wavelength of the pulse in nm.
    fwhm_ps : float
        The full-width-at-half-maximum of the pulse in picoseconds.
    time_window_ps : float
        The time window in picoseconds. This is the full-width of the time
        window, so the times will go from ``-time_window_ps * 0.5`` to
        + ``time_window_ps * 0.5``.
    power : float
        this is either the peak power of the pulse or the average power of
        the pulse-train, depending on ``power_is_avg``. In both cases the units
        are in watts.
    epp : float or None
        the energy-per-pulse in Joules. If this is not None (the default),
        then this overrides the power argument to set the pulse energy.
    npts : int
        the number of points in the time and frequency grids. Using powers
        of 2 might be beneficial for the efficiency of the FFTs in the
        NLSE algorithm. ``2**12`` is a good starting point.
    power_is_avg : boolean
        determines if the power is peak power or average power. Note that this 
        value is only used for calculating the initial pulse amplitude and is
        not saved as an intrinsic characteristic of the pulse object.
    frep_MHz : float
        the repetition rate in MHz. Used for setting the peak power of the
        pulse to match the average power of the pulse train. Similar to 
        ``power_is_avg``, the rep rate is not saved as part of the pulse class.
    GDD : boolean
        the group-delay-dispersion (in ps^2) to apply to the pulse.
    TOD : float
        the third-order dispersion (in ps^3) to apply to the pulse.
    FOD : float
        the fourth-order dispersion (in ps^3) to apply to the pulse.
    """

    def __init__(self, pulse_type='sech',length=1, center_wavelength_nm=1550, 
                 fwhm_ps=0.2, time_window_ps=10.0, power=1, epp=None,
                 npts=2**12, power_is_avg=False, frep_MHz=100, GDD=False, TOD=0,
                 FOD=0, RI=1.5, wav=[]):

        self.npts = npts
        self.center_wavelength_nm = center_wavelength_nm
        self.time_window_ps = time_window_ps

        T0_ps = fwhm_ps/1.76

        if pulse_type == 'sech':
            # from https://www.rp-photonics.com/sech2_shaped_pulses.html
            self.at = np.sqrt(power)/np.cosh(self.t_ps/(T0_ps))

        elif pulse_type == 'gaussian':
            # from https://www.rp-photonics.com/gaussian_pulses.html
            self.at = (np.sqrt(power) *
                       np.exp(-2.77*0.5*self.t_ps**2/(T0_ps**2)))

        elif pulse_type == 'sinc':
            T0_ps = fwhm_ps/3.7909885   # previously: T0_ps = FWHM_ps/3.7909885
            # numpy.sinc is sin(pi*x)/(pi*x), so we divide by pi
            self.at = np.sqrt(power) * np.sinc(self.t_ps/(T0_ps*np.pi))

        else:
            raise ValueError('Pulse type not recognized.')
            
        # GDD response
        if GDD is True:
            length = length
            wavelength = wav*1e9
            
            index = np.searchsorted(wavelength, self.center_wavelength_nm*1e9)
            
            dn = RI[index-1] - RI[index+1]
            dfreq = 2*np.pi*c_mks/(wavelength[index-1]-wavelength[index+1])
            
            T_g = (length/c_mks)*(dn/dfreq)
            #GDD_c = -T_g/self.w0_THz
           
            GDD_c = T_g/dfreq
            print(GDD_c)
        else:
            GDD_c = 0
            
        # TOD response
        if TOD is True:          
            TOD_c = GDD_c/dfreq
            print(TOD_c)
        else:
            TOD_c = 0
            
        # FOD response
        if FOD is True:
            FOD_c = TOD_c/dfreq
            print(FOD_c)
        else:
            FOD_c = 0

        if power_is_avg:
            self.at = self.at * np.sqrt(power / (frep_MHz*1.0e6 * self.epp))

        if epp is not None:
            self.epp = epp

        if pulse_type in ['sech', 'gaussian', 'sinc']:
            self.chirp_pulse_W(GDD_c, TOD_c, FOD_c)

    # FUNDAMENTAL PROPERTIES
    # These are the 4 fundamental properties that describe the pulse.
    # Everything else is derived from these.
    # Note that underscores are used for the actual variables (self._npts),
    # while self.npts redirects to the getter/setter methods.

    # 1. npts:
    @property
    def npts(self):
        """Set/get Number of points (int) in the time (and frequency) grid."""
        return self._npts

    @npts.setter
    def npts(self, new_npts):
        self._npts = int(new_npts)

    # 2. center frequency/wavelength:
    @property
    def centerfrequency_THz(self):
        """Set/get The center frequency (float) of the pulse in THz."""
        return self._centerfrequency_THz

    @centerfrequency_THz.setter
    def centerfrequency_THz(self, new_centerfrequency):
        assert (new_centerfrequency > 1.0)
        self._centerfrequency_THz = new_centerfrequency

    # 3. time window:
    @property
    def time_window_ps(self):
        """Set/get the time window of the time grid in picoseconds (float)."""
        return self._time_window_ps

    @time_window_ps.setter
    def time_window_ps(self, time_window_ps):
        self._time_window_ps = time_window_ps

    # 4. amplitude in frequency domain:
    @property
    def aw(self):
        """Set/get the complex amplitude of the pulse in the frequency domain.

        Corresponds to the frequencies in the pulse.aw array. Note that this
        is the complex amplitude and that the intensity is the square of the
        absolute value.
        
        The units are sqrt(W) or sqrt(J)*sqrt(Hz), so abs(aw)^2 * deltaF will
        provide J/bin. (90% sure this comment about units is correct.)
        """
        return self._aw

    @aw.setter
    def aw(self, aw_new):
        if 'self._aw' not in locals():  # if aw doesn't exist, make blank array
            self._aw = np.zeros((self._npts,), dtype=np.complex128)
        self._aw[:] = aw_new

    # DERIVED PROPERTIES:

    # center angular frequency
    @property
    def w0_THz(self):
        """Get the center *angular* frequency (THz)."""
        return 2.0 * np.pi * self._centerfrequency_THz

    # center wavelength:
    @property
    def center_wavelength_nm(self):
        """Set/get the center wavelength of the grid in units of nanometers."""
        return c_nmps / self._centerfrequency_THz

    @center_wavelength_nm.setter
    def center_wavelength_nm(self, wl):
        self.centerfrequency_THz = c_nmps / wl

    # frequency grid
    @property
    def v_THz(self):
        """Get the *relative* *angular* frequency grid in THz."""
        return 2.0*np.pi*np.arange(-self.npts/2, self.npts/2)/(self.npts *
                                                               self.dt_ps)

    @property
    def w_THz(self):
        """Get the absolute *angular* frequency grid (THz)."""
        return self.v_THz + self.w0_THz

    @property
    def f_THz(self):
        """Get the absolute frequency grid in THz."""
        return (self.v_THz + self.w0_THz)/(2 * np.pi)

    # wavelength grid:
    @property
    def wavelength_nm(self):
        """Wavelength grid in nanometers."""
        return c_nmps / self.f_THz

    # time grid:
    @property
    def t_ps(self):
        """Get the temporal grid in ps."""
        return np.linspace(-self._time_window_ps / 2.0,
                           self._time_window_ps / 2.0,
                           self._npts, endpoint=False)

    # dt:
    @property
    def dt_ps(self):
        """Return time grid spacing in ps."""
        return self._time_window_ps / np.double(self._npts)

    @property
    def df_THz(self):
        """Frequency grid spacing in THz."""
        f_THz = self.f_THz
        return f_THz[1] - f_THz[0]

    # amplitude in the time domain:
    @property
    def at(self):
        """Amplitude of the time-domain electric field (complex).

        Units are ``sqrt(W)``, which can be considered ``sqrt(J)/sqrt(s)``,
        so units of energy per time bin when the absolute value is squared.
        """
        return IFFT_t(self._aw.copy())

    @at.setter
    def at(self, at_new):
        self.aw = FFT_t(at_new)

    # epp, energy per pulse:
    @property
    def epp(self):
        """Energy per pulse in Joules."""
        return (self.dt_ps * 1e-12) * np.trapz(abs(self.at)**2)

    @epp.setter
    def epp(self, desired_epp_J):
        self.at = self.at * np.sqrt(desired_epp_J / self.epp)
        
        
    # intensity in various units
    def psd(self, units='mW', rep_rate=1):
        """Return the power spectral density (PSD) in various units. Set the 
        rep_rate to 1 to get "per pulse" units. Otherwise, the rep-rate will be
        used to scale the per-pulse numbers to average power units. By default,
        the rep-rate of the pulse object is used.
        
        Unit options are:
        
        ``'mW'``, mW for each data point. (Not actually PSD units.) 
        ``'mW/THz'``, mW per THz. 
        ``'dBm/THz'``, 10*log10(mW) per THz. 
        ``'mW/nm'``, mW per nanometer.
        ``'dBm/nm'``, 10*log10(mW) per nanometer.
        
        Note that for the "per nanometer" units, the data is still delivered on
        a grid that is evenly spaced in *frequency*, not wavelength. So, if
        integrating the PSD, it is necessary to take into account the changing
        size of the wavelength bins to recover the correct value for the
        average power. The psd_wavelength function provides both the evenly
        spaced wavelength grid and the y-axis unit conversion.

        Parameters
        ----------
        units : str
            Determines the units of the intensity of power-spectral-density.
            See above for options.
        rep_rate : float
            Determines the repetition rate (in Hz) used to calculate the PSD.
            Set to 1 to get per-pulse PSD. If ``None``, then the rep-rate for
            the pulse object will be used. Note that the rep-rate provided to
            this function doesn't change the rep-rate of the pulse object, it
            is merely used for the PSD calculation.
        
        
        Returns
        -------
        psd : array
            numpy array of power spectral densities corresponding to pulse.aw.
        
        
        """
        
        if rep_rate == 1:
            print('Note: the rep-rate is set to 1; per-pulse values returned.')

        # Note: all units are calculated since these calculations are easy.
            
        # per bin units:
        f = self.f_THz
        df = (f[1]-f[0]) * 1e12  # df in Hz
        J_Hz = np.abs(self.aw)**2
        J_per_bin = J_Hz / df  # go from J*Hz/bin (native units) to J/bin
        # multiply by rep rate to get W/bin, and then mW/bin:
        mW_per_bin = J_per_bin * rep_rate * 1e3  
        
        # per THz units:
        mW_per_THz = mW_per_bin / (df * 1e-12)
        dBm_per_THz = 10 * np.log10(mW_per_THz)
            
        # per wavelength units
        wl_nm = c_nmps / f
        wl_m = wl_nm * 1e-9 
        nm_per_bin = wl_m**2 / 3e8 * df * 1e9  # Jacobian from THz to nm
        mW_per_nm = mW_per_bin / nm_per_bin # convert to mW/nm
        dBm_per_nm = 10 * np.log10(mW_per_nm)  # convert to dBm/nm
        
        # return the requested PSD:
        if   units == 'mW/bin':
            return mW_per_bin
        elif units == 'mW/THz':
            return mW_per_THz
        elif units == 'dBm/THz':
            return dBm_per_THz
        elif units == 'mW/nm':
            return mW_per_nm
        elif units == 'dBm/nm':
            return dBm_per_nm
        else:
            raise ValueError('Units not recognized.')
        
    def psd_wavelength(self, wl_min=500, wl_max=2500, wl_step=0.2):
        """-- Not yet implemented -- 
        evenly spaced wavelength and intensity in various units"""
        
        pass
       
        
    def clone_pulse(self, p):
        """Copy all parameters of pulse_instance into this one."""
        self.npts = p.npts
        self.centerfrequency_THz = p.centerfrequency_THz
        self.time_window_ps = p.time_window_ps
        self.aw = p.aw

    def create_cloned_pulse(self):
        """Create and return new pulse instance identical to this instance."""
        p = Pulse()
        p.clone_pulse(self)
        return p


    def add_noise(self, noise_type='sqrt_N_freq'):
        """Add random intensity and phase noise to a pulse.

        Parameters
        ----------
        noise_type : string
            The method used to add noise. The options are:

            - ``sqrt_N_freq`` : adds noise to each bin in the frequency domain.
              The average noise added is proportional to sqrt(N), and where N
              is the number of photons in that frequency bin.

            - ``one_photon_freq``` : which adds one photon of noise to each
              frequency bin.
        """
        # Get the number of photons/second in each frequency bin:
        size_of_bins = self.df_THz * 1e-12                 # Bin width in [Hz]
        power_per_bin = np.abs(self.aw)**2 / size_of_bins  # [J*Hz]/[Hz] = [J]

        h = 6.62607004e-34

        photon_energy = h * self.f_THz * 1e-12  # h nu [J]
        photons_per_bin = power_per_bin/photon_energy  # photons / second
        photons_per_bin[photons_per_bin < 0] = 0  # must be positive.

        # now generate some random intensity and phase arrays:
        size = np.shape(self.aw)[0]
        random_intensity = np.random.normal(size=size)
        random_phase = np.random.uniform(size=size) * 2 * np.pi

        if noise_type == 'sqrt_N_freq':
            # Gaussian noise with sqrt(photons_per_bin)
            noise = (random_intensity * np.sqrt(photons_per_bin) *
                     photon_energy * size_of_bins * np.exp(1j*random_phase))

        elif noise_type == 'one_photon_freq':  # one photon per bin
            noise = (random_intensity * photon_energy * size_of_bins *
                     np.exp(1j*random_phase))
        else:
            raise ValueError('noise_type not recognized.')

        self.aw = self.aw + noise
        
        
    def chirp_pulse_W(self, GDD=0, TOD=0, FOD=0.0, w0_THz=None):
        r"""Alter the phase of the pulse.

        Apply the dispersion coefficients :math:`\beta_2, \beta_3, \beta_4`
        expanded around frequency :math:`\omega_0`.

        Parameters
        ----------
        GDD : boolean
            Group delay dispersion (:math:`\beta_2`) [ps^2], defaults to 0.
        TOD : boolean, optional
            Third order dispersion (:math:`\beta_3`) [ps^3], defaults to 0.
        FOD : boolean, optional
            Fourth order dispersion (:math:`\beta_4`) [ps^4], defaults to 0.
        w0_THz : float, optional
            Center freq. of dispersion expansion, defaults to grid center freq.

        Notes
        -----
        The convention used for dispersion is

        .. math:: E_{new} (\omega) = \exp\left(i \left(
           \frac{1}{2} GDD\, \omega^2 + \frac{1}{6}\, TOD \omega^3 +
           \frac{1}{24} FOD\, \omega^4 \right)\right) E(\omega)
        """
        if w0_THz is None:
            self.aw = np.exp(1j * (GDD / 2.0) * self.v_THz**2 +
                             1j * (TOD / 6.0) * self.v_THz**3 +
                             1j * (FOD / 24.0) * self.v_THz**4) * self.aw
        else:
            V = self.w_THz - w0_THz
            self.aw = np.exp(1j * (GDD / 2.0) * V**2 +
                             1j * (TOD / 6.0) * V**3 +
                             1j * (FOD / 24.0) * V**4) * self.AW

    def calc_width(self, level=0.5):
        """Calculate the pulse width.
        
        For example, the full-width-at-half-maxmimum
        (FWHM) or the 1/e width. If the pulse has multiple crossings (for example,
        if it reaches the 0.5 level multiple times) the widest extent will be
        returned.

        Parameters
        ----------
        level : float
            the fraction of the peak to calculate the width. 
            Must be between 0 and 1.
            Default is 0.5, which provides the FWHM.
            1/e = 0.367879
            1/e^2 = 0.135335

        Returns
        -------
        width : float
            the width of the pulse in picoseconds.
        """
    
        def find_roots(x, y):
            s = np.abs(np.diff(np.sign(y))).astype(bool)
            return x[:-1][s] + np.diff(x)[s]/(np.abs(y[1:][s]/y[:-1][s])+1)
        
        it = np.abs(self.at)**2
        it = it/np.max(it)
        roots = find_roots(self.t_ps, it - level)
        width = np.max(roots) - np.min(roots)
        return(width)
    
    def transform_limit(self):
        """"Return a transform-limited pulse (flat spectral phase)."""
        newpulse = self.create_cloned_pulse()
        newpulse.aw = np.abs(newpulse.aw)
        return newpulse
    

def FFT_t(A, ax=0):
    """Do a FFT with fft-shifting."""
    A = A.astype('complex128')
    return fft.ifftshift(fft.ifft(fft.fftshift(A, axes=(ax,)),
                                  axis=ax), axes=(ax,))


def IFFT_t(A, ax=0):
    """Do an iFFT with fft-shifting."""
    A = A.astype('complex128')
    return fft.ifftshift(fft.fft(fft.fftshift(A, axes=(ax,)),
                                 axis=ax), axes=(ax,))