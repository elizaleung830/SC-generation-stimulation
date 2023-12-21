import femwell
from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import compute_modes
from skfem import Basis, ElementTriP0
from skfem.io import from_meshio
import numpy as np

def get_neff_and_aeff(polygon,n_dict:dict, wavelength_range:list, step, resolutions, default_resolution_max=2, aeff_equation = 2):
    """
    :param polygon:
    :param n_dict:
    :param wavelength_range: wavelength in nm
    :param step:
    :param resolutions:
    :param default_resolution_max:
    :param aeff_equation:
    :return: aeff_list, neff_list: numpy array object of the
    """
    mesh = from_meshio(mesh_from_OrderedDict(polygon, resolutions, default_resolution_max=default_resolution_max))
    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros()
    for subdomain, n in n_dict.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n ** 2

    wavelength_list = np.linspace(wavelength_range[0], wavelength_range[1], step)
    neff_list = []
    aeff_list = []
    for wavelength in wavelength_list:
        wavelength = wavelength *1e-3
        modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=1, order=1)
        for mode in modes:
            if mode.tm_fraction > 0.5:
                neff_list.append(np.real(mode.n_eff))
                aeff_list.append(mode.calculate_effective_area(equation=aeff_equation))
                break
    return np.array(aeff_list), np.array(neff_list), np.array(wavelength_list)