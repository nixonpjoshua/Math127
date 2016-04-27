import numpy as np

# note b_i's are slightly difffrom class
# they have already been cleaned up

"""
Computes "b"
Args: 
    I: the image i.e. a 2d array calculated from fst_mol
    
Returns:
    "b" a 3D array used in computing noisy image
"""
def compute_b(I):
    len_I = I.shape[0]
    # len_Image need not equal D keeping it general
    ans = np.zeros((len_I, len_I, len_I))
    for x in np.arange(len_I):
            for y in np.arange(len_I):
                for z in np.arange(len_I):
                # Compute b_i
                fourier_im = np.fft.fft2(I)[x][y]
                sinc_term  = np.div(np.sin(len_image*np.pi*z),pi*z)
                ans[x][y][z] = fourier_im*sinc_term
    return ans 

"""
Computes "B"
Args: 
    images: A list of Images i.e. a 2d arrays calculated from fst_mol
     
Returns:
    "noisy" a 3D array representing a sampled function, 
    unfiltered reconstruction
"""
def compute_noisy(images):
    num_ims = len(images)
    b_is    = np.zeros(num_ims)
    ans     = np.zeros((num_ims,num_ims))
    for i in np.arange(num_ims):
        ans = ans + compute_B(images[i])
    return ans


"""
Computes no
Args: 
    images:    A list of Images 
    rotations: A list of rotation matrices 
    D: size of the image i.e. mol.shape[0]
     
Returns:
    DxDxD Array that represents the 3D reconstruction of the 
    m
"""

def compute_fltr(rotations):
    num_rots = len(rotations)
    for i in len(rot):
    return


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


def back_project(D, images, rotations):
    noisy     = compute_noisy(images)
    fltr      = compute_fltr(rotations)
    filtered  = noisy * fltr
    # is filtered but still need to perform a base change
    mol_hat   = np.fft.ifftn()
    return mol_hat



 