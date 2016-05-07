import numpy as np

# note b_i's are slightly diff from class
# they have already been cleaned up

# Introduce h_i 's for filter similar to 
# b_i for noisy image
# TODO: Conceptually understand what subparts are 
#       so we can come up with better names

def compute_b(I):
    """
    Computes "b"
    Args:
        I: the image i.e. a 2d array calculated from fst_mol

    Returns:
        "b" a 3D array used in computing noisy image
    """
    len_I = I.shape[0]
    # len_Image need    # fourier_image_ext = np.array([fourier_image for i in xrange(len_I)]) not equal D keeping it general
    # TODO ask, for our code it is ok if its always
    # equal to D
    fourier_image = np.fft.fft2(I)
    z = np.arange(len_I)
    sinc_terms = np.complex128(np.sinc(len_I*z)).reshape(len_I, 1, 1)
    # TODO note that now at z[0] it evaluates to len_I, not sure which is correct should revisit equation

    # z_list = []
    # for z in xrange(len_I):
    #     # Compute b_i
    #     if z == 0:
    #         sinc_term = 1 #TODO note that at zero we used to have it be one
    #     else:
    #         sinc_term = np.sinc(len_I*z)*len_I
    #         # sinc_term  = np.sin(len_I*np.pi*z)/np.pi*z
    #     z_list.append(np.complex128(sinc_term))
    # z_list = np.array(z_list).reshape(len_I, 1, 1)
    return sinc_terms*fourier_image


def compute_noisy(images):
    """
    Computes "B"
    Args:
        images: A list of Images i.e. 2d arrrepresentingays calculated from fst_mol

    Returns:
        "noisy" a 3D array representing a sampled function,
        unfiltered reconstruction
    """
    ans = np.zeros(images[0].shape)
    for i in xrange(len(images)):
        ans = ans + compute_b(images[i])
    return ans


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


def back_project(D, images, rotations):
    """
    Creates a 3D backprojection using images gathered from fst_mol
    Args:
        images:    A list of Images
        rotations: A list of rotation matrices
        D: size of the image i.e. mol.shape[0]

    Returns:
        DxDxD Array that represents the 3D reconstruction of the
        molecule
    """
    # is filtered but still need to perform a base change
    return np.real(np.fft.ifftn(compute_noisy(images) * compute_fltr(rotations, images[0].shape[0])))
