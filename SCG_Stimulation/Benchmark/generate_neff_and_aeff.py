from femwell.visualization import plot_domains

import femwell
from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import compute_modes
from skfem import Basis, ElementTriP0
from skfem.io import from_meshio
import numpy as np
import matplotlib.pyplot as plt
from alive_progress import alive_bar
from time import sleep
from tqdm import tqdm

def get_neff_and_aeff(polygon,n_dict:dict, wavelength_range:list, step, resolutions, default_resolution_max=2, aeff_equation = 2, plot = True):
    """
    :param polygon:
    :param n_dict:
    :param wavelength_range: wavelength in nm
    :param step:
    :param resolutions:
    :param default_resolution_max:
    :param aeff_equation:
    :param plot: plot the mesh
    :return: aeff_list, neff_list: numpy array object of the
    """
    mesh = from_meshio(mesh_from_OrderedDict(polygon, resolutions, default_resolution_max=default_resolution_max))
    if plot == True:
        mesh.draw().show()
        plot_domains(mesh)
        plt.show()

    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros()
    wavelength_list = np.linspace(wavelength_range[0], wavelength_range[1], step)
    neff_list = []
    aeff_list = []
    count = 0
    for wavelength in tqdm(wavelength_list):
        wavelength = wavelength * 1e-3
        for subdomain, n in n_dict.items():
            epsilon[basis0.get_dofs(elements=subdomain)] = n(wavelength) ** 2
        modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=3, order=1)
        for mode in modes:
            print(f"effective area:{mode.calculate_effective_area()}")
            if mode.te_fraction > 0.5:
                count += 1
                neff_list.append(np.real(mode.n_eff))
                aeff_list.append(mode.calculate_effective_area())
                break

    return np.array(aeff_list), np.array(neff_list), np.array(wavelength_list)