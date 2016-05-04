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
    # len_Image need not equal D keeping it general
    # TODO ask, for our code it is ok if its always
    # equal to D
    fourier_image = np.fft.fft2(I)
    ans = np.zeros((len_I, len_I, len_I), dtype=np.complex128)
    for x in xrange(len_I):
        for y in xrange(len_I):
            for z in xrange(len_I):
                # Compute b_i
                if z == 0:
                    sinc_term = 1 #TODO check that this is the right thing to do
                else:
                    sinc_term  = np.sin(len_I*np.pi*z)/np.pi*z

                ans[x][y][z] = fourier_image[x][y]*np.complex128(sinc_term)
    return ans


def compute_noisy(images):
    """
    Computes "B"
    Args:
        images: A list of Images i.e. 2d arrrepresentingays calculated from fst_mol

    Returns:
        "noisy" a 3D array representing a sampled function,
        unfiltered reconstruction
    """
    num_ims = len(images)
    ans = np.zeros(images[0].shape)
    for i in xrange(num_ims):
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
    ans = 0
    # TODO double check with Dynerman this is how 
    # how to compute each part of the sum for the 
    # filter
    for x in xrange(size):
        for y in xrange(size):
            for z in xrange(size):
                # x refers to input of the function
                pos = np.array([x, y, z])
                x = size*np.pi*np.dot(pos, c_vec)
                if x == 0:
                    ans += size
                else:
                    ans += size*np.sin(x)/x
    return ans 


def compute_fltr(rots, size):
    """
    Computes filter for taking noise out of images
    Args:
        rots: A list of rotations

    Returns:
        filter is a scalar that gives that filters noisy image
    """
    num_rots = len(rots)
    ans = 0
    for i in xrange(num_rots):
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
    noisy = compute_noisy(images)
    fltr = compute_fltr(rotations, images[0].shape[0])
    # is filtered but still need to perform a base change
    mol_hat = np.fft.ifftn(noisy * fltr)
    length = mol_hat.shape[0]
    ans = np.zeros((length, length, length))
    for x in xrange(length):
        for y in xrange(length):
            for z in xrange(length):
                ans[x, y, z] = np.linalg.norm(mol_hat[x, y, z])
    return ans