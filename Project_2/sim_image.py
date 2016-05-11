import numpy as np
import matplotlib as plt
from scipy.interpolate import RegularGridInterpolator


def project_fst(mol, R):
    """
    :param mol: an NxNxN array that contains rho
    :param R: the rotation matrix to use for the image projection
    :return: an image projected using the rotation matrix R and the data in mol
    """
    N = mol.shape[0]
    N_range = np.linspace(-1, 1, N)
    inter = RegularGridInterpolator((N_range, N_range, N_range), mol, method='linear', bounds_error=False, fill_value=0)
    pos = np.concatenate((N_range.reshape(N, 1, 1, 1) * np.ones(N * N).reshape(1, N, N, 1),
                          N_range.reshape(1, N, 1, 1) * np.ones(N * N).reshape(N, 1, N, 1),
                          N_range.reshape(1, 1, N, 1) * np.ones(N * N).reshape(N, N, 1, 1)),
                         axis=3)
    # pos gives a matrix of size N which at [x,y,z] returns a value which we can transform into the image plane space
    r_mol = inter(np.dot(pos, R))
    return np.real(np.fft.ifft2(np.fft.fftn(r_mol)[:, :, 0]))
