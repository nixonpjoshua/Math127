import numpy as np

#######################################################################
## got from http://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html ##
#######################################################################

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np


test = np.array([[[ 0,  1,  2],
        [ 3,  4,  5],
        [ 6,  7,  8]],
       [[ 9, 10, 11],
        [12, 13, 14],
        [15, 16, 17]],
       [[18, 19, 20],
        [21, 22, 23],
        [24, 25, 26]]])

x = test[0]
y = test[1]
z = test[2]


print("--------------------------------")
print("the bottom slice is:")
print(x)
print("--------------------------------")
print("the middle slice is:")
print(y)
print("--------------------------------")
print("the top slice is:")
print(z)
print("--------------------------------")

	



# -*- coding: utf-8 -*-
"""
slider 3D numpy array

"""

import numpy
import pylab
from matplotlib.widgets import Slider

data = numpy.random.rand(100,256,256) #3d-array with 100 frames 256x256

ax = pylab.subplot(111)
pylab.subplots_adjust(left=0.25, bottom=0.25)

frame = 0
l = pylab.imshow(data[frame,:,:]) #shows 256x256 image, i.e. 0th frame

axcolor = 'lightgoldenrodyellow'
axframe = pylab.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
sframe = Slider(axframe, 'Frame', 0, 100, valinit=0)

def update(val):
    frame = numpy.around(sframe.val)
    pylab.subplot(111)
    pylab.subplots_adjust(left=0.25, bottom=0.25)
    pylab.imshow(data[frame,:,:])

sframe.on_changed(update)

pylab.show()





































