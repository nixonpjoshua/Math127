import numpy as np


def l(x, y, z, D):
    """
    :param x, y, z: coordinates
    :param D: size of input
    :return:
    """
    #TODO wtf is delta?
    return delta(x,y)*rect(z, D)


def rect(z, D):
    """
    :param z: coordinate
    :param D: size of input
    :return: is rect 1 or not
    """
    #TODO how big should this matrix be, right now its size D?
    if -D/2 <= z <= D/2:
        return 1
    else:
        return 0


def b_i(x, y, z, D, R, I):
    """
    :param x, y, z: coordinates
    :param R: transform matrix for i
    :param I: is the image
    :param D: the size of input
    :return: calculates little b for a given image position
    """
    I_3D = R[:, 0]*x + R[:, 1]*y
    l_val = l(x, y, z, D)
    return realify(convolution3D(I_3D, l_val))


def realify(m):
    one = m.shape[0]
    two = m.shape[1]
    three = m.shape[2]
    out = np.zeros((one, two, three))
    for i in xrange(one):
        for j in xrange(two):
            for k in xrange(three):
                # TODO should this be norm or just the real part
                out[i, j, k] = np.linalg.norm(conv[i, j, k])
    return out

def back_project(D, images, rotations):
    """
    :param D: size of mol
    :param images: all the images
    :param rotations: parallel array to images giving the orientation
    :return:
    """
    B = np.zeros((D, D, D))
    H = np.zeros(B.shape)
    for x in xrange(D):
        for y in xrange(D):
            for z in xrange(D):
                B_val = np.zeros(images[0].shape)
                H_val = np.zeros(images[0].shape)
                for index in xrange(len(images)):
                    I = images[index]
                    R = rotations[index]
                    B_val += np.fft.fftn(b_i(x, y, z, D, R, I))
                    #TODO should there be pi in here or does sinc do that for us
                    H_val += np.sinc(D*np.pi*np.dot(np.array([x, y, z]), R[:, 2]))
                B[x,y,z] = B_val
                H[x,y,z] = H_val
    return realify(np.fft.ifftn(B/H))



def convolution3D(f, g):
    """
    this assumes that f.shape[n] = g.shape[n] and that they are both 3d
    note that due to ifft this will return an array of complex numbersg always
    :param f: function 1
    :param g: function 2
    :return: the convolution of f and g computed using the fft
    """
    f = np.fft.fftn(f)
    g = np.fft.fftn(g)
    I = f.shape[0]
    J = f.shape[1]
    K = f.shape[2]
    res = np.zeros((I, J, K), dtype=np.complex128)
    for i in xrange(I):
        for j in xrange(J):
            for k in xrange(K):
                res[i, j, k] = f[i, j, k] * g[i, j, k]
    return np.fft.ifftn(res)