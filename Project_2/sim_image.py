import numpy as np
import matplotlib as plt
from scipy.interpolate import RegularGridInterpolator


def project_fst(mol, R):
    """
    Creates an "image" given NxNxN matrix representation of mol and rotation
    matrix R
    Args:
         mol: An NxNxN array that serves as an approximation for a function
              from R^3 --> R
         R:   rotation matrix [a,b,c] (a,b,c correspond to row vectors)

    Returns:
        NxN 2darray "Image"
    """
    N = mol.shape[0]
    N_range = np.linspace(-1, 1, N)
    inter = RegularGridInterpolator((N_range, N_range, N_range), mol, method='linear', bounds_error=False, fill_value=0)
    pos = np.concatenate((N_range.reshape(N, 1, 1, 1) * np.ones(N * N).reshape(1, N, N, 1),
                          N_range.reshape(1, N, 1, 1) * np.ones(N * N).reshape(N, 1, N, 1),
                          N_range.reshape(1, 1, N, 1) * np.ones(N * N).reshape(N, N, 1, 1)),
                         axis=3)
    r_mol = inter(np.dot(pos, R))
    return np.real(np.fft.ifft2(np.fft.fftn(r_mol)[:, :, 0]))
