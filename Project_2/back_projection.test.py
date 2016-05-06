import numpy as np
from MRCFile import MRCFile
from sim_image import *
from back_projection import *
import matplotlib.pyplot as plt
import matplotlib.image
import sys

#from mayavi.mlab import * TODO talk to Dynerman see if can get this working

# Cooking up our own set of rotation matrices
# For more information see 
# https://en.wikipedia.org/wiki/Rotation_matrix


"""
Cooks up "evenly spaced" rotation matrices that rotate 
the yz plane about the x direction

Args: 
    n: number of rotation matrices 
     
Returns:
    list of rotation matrices
"""


def rotate_x(n):
	# range of angle is from 0 to pi
	step = np.pi/n
	# + 1e-5 hacky fix to make range inclusive
	thetas = np.arange(0,np.pi + 1e-5 ,step) 
	# TODO fix for efficieny purposes
	matrices = []
	for j in np.arange(len(thetas)):
		i = thetas[j]
		new_matrix = np.array([[ 1, 0, 0],
							   [ 0 , np.cos(i), -np.sin(i)],
							   [ 0 , np.sin(i), np.cos(i)]
							   ])
		matrices.append(new_matrix)
	return np.array(matrices)

"""
Cooks up "evenly spaceed" rotation matrices that rotate 
the xz plane about the y direction

Args: 
    n: number of rotation matrices 
     
Returns:
    np array of rotation matrices
"""
def rotate_y(n):
	# range of angle is from 0 to pi
	step = np.pi/n
	# + 1e-5 hacky fix to make range inclusive
	thetas = np.arange(0,np.pi + 1e-5 ,step) 
	# TODO fix for efficieny purposes
	matrices = []
	for j in np.arange(len(thetas)):
		i = thetas[j]
		new_matrix = np.array([[ np.cos(i), 0, np.sin(i)],
							   [ 0 , 1, 0],
							   [ -np.sin(i) , 0, np.cos(i)]
							   ])
		matrices.append(new_matrix)
	return np.array(matrices)

"""
Cooks up "evenly spaceed" rotation matrices that rotate 
the xz plane about the y direction

Args: 
    n: number of rotation matrices 
     
Returns:
    np array of rotation matrices
"""
def rotate_z(n):
	# range of angle is from 0 to pi
	step = np.pi/n
	# + 1e-5 hacky fix to make range inclusive
	thetas = np.arange(0,np.pi + 1e-5 ,step) 
	# TODO fix for efficieny purposes
	matrices = []
	for j in np.arange(len(thetas)):
		i = thetas[j]
		new_matrix = np.array([[ np.cos(i), -np.sin(i), 0],
							   [ np.sin(i) , np.cos(i), 0],
							   [ 0 , 0, 1]
							   ])
		matrices.append(new_matrix)
	return np.array(matrices)

zika_153 = MRCFile('zika_153.mrc')  # initializing object
zika_153.load_all_slices()
mol = zika_153.slices






##############################
#    Test Script			#
##############################






zika_153 = MRCFile('zika_153.mrc')  # initializing object
zika_153.load_all_slices()
mol = zika_153.slices

rot_x = rotate_x(1)
rot_y = rotate_y(1)
rot_z = rotate_z(1)

rots = [0]*3

for i in xrange(3):
	if i < 1:
		rots[i] = rot_x[i]
	elif i < 2:
		rots[i] = rot_y[i-1]
	else:
		rots[i] = rot_z[i-2]

images = [0]*3
for i in np.arange(len(rots)):
	images[i] = project_fst(mol,rots[i])

# b
# 153 hardcoded for zika 153
#

test   = back_project(153, [images[0]], [rots[0]])
print(test)
visual = write_file(test)





















