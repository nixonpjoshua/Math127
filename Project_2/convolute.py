import numpy as np
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