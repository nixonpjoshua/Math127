import numpy as np
from scipy.interpolate import RegularGridInterpolator


def compute_h(rot, size):
    """
    Computes "h"
    Args:
        rot: a rotation matrix
        size: size of the image

    Returns:
        "h" a constant used in computing noisy image
    """
    c_vec = rot[:, 2]
    # c_vec is 3rd column of rotation matrix

    pos = np.concatenate((np.arange(size).reshape(size, 1, 1, 1)*np.ones(size*size).reshape(1, size, size, 1),
                          np.arange(size).reshape(1, size, 1, 1)*np.ones(size*size).reshape(size, 1, size, 1),
                          np.arange(size).reshape(1, 1, size, 1)*np.ones(size*size).reshape(size, size, 1, 1)), axis=3)
    # pos[x,y,z] = [x,y,z], returns its own index
    return np.sum(size*np.sinc(size * np.pi * np.dot(pos, c_vec)))

def compute_fltr(rots, size):
    """
    Computes filter for taking noise out of images
    Args:
        rots: A list of rotations

    Returns:
        filter is a scalar that gives that filters noisy image
    """
    ans = 0
    for i in xrange(len(rots)):
        ans += compute_h(rots[i], size)
    return ans


def back_project(data, use_filter = True, ):
    N = data[0][0].shape[0]
    r = np.zeros(N)
    r[N/4:(N/4)*3] = 1
    s = np.fft.fftshift(np.fft.fft(r))
    s = np.tile(s[np.newaxis, np.newaxis, :], (N, N, 1))
    B = np.zeros((N, N, N))
    if use_filter:
        H = np.zeros((N, N, N))
    for image, R in data:
        I_hat = np.fft.fftshift(np.fft.fft2(image))
        I_hat = np.tile(I_hat[..., np.newaxis], (1, 1, N))
        image = np.real(np.fft.ifftn(np.fft.ifftshift(I_hat*s)))
        # rotate and interpolate
        N_range = np.linspace(-1, 1, N)
        inter = RegularGridInterpolator((N_range, N_range, N_range), image, method='linear', bounds_error=False, fill_value=0)

        x, y, z = np.meshgrid(N_range, N_range, N_range)
        C = [x.flatten(), y.flatten(), z.flatten()]
        B += inter(np.dot(R.transpose(), C).T).reshape(N, N, N)
        if use_filter:
            H += compute_h(R, N)
    if use_filter:
        B_hat = np.fft.fftshift(np.fft.fftn(B))
        return np.real(np.fft.ifftn(np.fft.ifftshift(B_hat/H)))
    else:
        return B
