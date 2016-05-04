import numpy
from mayavi.mlab import *

def test_contour3d():
    x, y, z = numpy.ogrid[-5:5:64j, -5:5:64j, -5:5:64j]

    scalars = x * x * 0.5 + y * y + z * z * 2.0

    obj = contour3d(scalars, contours=4, transparent=True)
    return obj

test_contour3d()