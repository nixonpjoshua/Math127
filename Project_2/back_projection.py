import numpy as np
from scipy.interpolate import RegularGridInterpolator


def compute_h(rot, D):
    """
    :param rot: the rotation matrix for a particular image
    :param D: the size of the rho
    :return: the component of h for this image's rotation matrix
    """
    c_vec = rot[:, 2]
    # c_vec is 3rd column of rotation matrix

    pos = np.concatenate((np.arange(D).reshape(D, 1, 1, 1)*np.ones(D*D).reshape(1, D, D, 1),
                          np.arange(D).reshape(1, D, 1, 1)*np.ones(D*D).reshape(D, 1, D, 1),
                          np.arange(D).reshape(1, 1, D, 1)*np.ones(D*D).reshape(D, D, 1, 1)), axis=3)
    # pos[x,y,z] = [x,y,z], returns its own index
    return np.sum(D*np.sinc(D * np.dot(pos, c_vec)))


def back_project(data, use_filter=True, D_percent=.5):
    """
    :param data: tuples containing images as np arrays and their rotations
    :param use_filter: whether or not to apply to transfer function filter
    :param D_percent: what percent of the total 'box' does the rho take up
    :return:
    """
    N = data[0][0].shape[0]
    r = np.zeros(N)
    D = int(D_percent*N)
    r[D/2:N-(D/2)] = 1
    s = np.fft.fftshift(np.fft.fft(r))
    # r is the rect function
    s = np.tile(s[np.newaxis, np.newaxis, :], (N, N, 1))
    B = np.zeros((N, N, N))
    if use_filter:
        H = 0
    for image, R in data:
        I_hat = np.fft.fftshift(np.fft.fft2(image))
        I_hat = np.tile(I_hat[..., np.newaxis], (1, 1, N))
        image = np.real(np.fft.ifftn(np.fft.ifftshift(I_hat*s)))
        # rotate and interpolate
        N_range = np.linspace(-1, 1, N)
        inter = RegularGridInterpolator((N_range, N_range, N_range), image, method='linear', bounds_error=False, fill_value=0)
        pos = np.concatenate((N_range.reshape(N, 1, 1, 1) * np.ones(N * N).reshape(1, N, N, 1),
                              N_range.reshape(1, N, 1, 1) * np.ones(N * N).reshape(N, 1, N, 1),
                              N_range.reshape(1, 1, N, 1) * np.ones(N * N).reshape(N, N, 1, 1)),
                             axis=3)
        # pos gives a matrix of size N which at [x,y,z] returns a vector which we can transform from the image plane space
        B += inter(np.dot(pos, R.transpose()))
        # we use R.transpose to transform all the vectors away from image plane space to the original
        if use_filter:
            H += compute_h(R, D)
    if use_filter:
        # H is defined in fourier space
        B_hat = np.fft.fftshift(np.fft.fftn(B))
        return np.real(np.fft.ifftn(np.fft.ifftshift(B_hat/H)))
    else:
        return B
