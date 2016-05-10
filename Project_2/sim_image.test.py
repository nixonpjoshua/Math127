import numpy as np
from MRCFile import MRCFile
from sim_image import *
import matplotlib.pyplot as plt
from vanheel import  random_rotation_matrix
import time
import matplotlib.image

"""
Playing around with the data adnd the syntax
"""

zika_153 = MRCFile('zika_153.mrc')  # initializing object
zika_153.load_all_slices()
mol = zika_153.slices

# print(type(mol))

# print(mol)

# Testing Whether centering the image function changes the fast fourier transform
# 3x3x3 array

# test = np.array([[[ 0,  1,  2],
#         [ 3,  4,  5],
#         [ 6,  7,  8]],
#        [[ 9, 10, 11],
#         [12, 13, 14],
#         [15, 16, 17]],
#        [[18, 19, 20],
#         [21, 22, 23],
#         [24, 25, 26]]])

# print("the bottom slice is", test[0])
# print("-------------------------------")
# print("the middle slice is", test[1])
# print("--------------------------------")
# print("the top slice is", test[2])
#
# R = np.array([[1,0,0],
#               [0,1,0],
#               [0,0,1]])
start = time.time()
image = project_fst(mol,random_rotation_matrix())
end = time.time()
print str(end - start)
plt.imshow(image)
plt.show()






































