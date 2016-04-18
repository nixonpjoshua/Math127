import numpy as np 
from MRCFile import MRCFile

"""
Playing around with the data adnd the syntax
"""

zika_153 = MRCFile('zika_153.mrc') # initializing object
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

"""
Playing around with utility functions
"""
def center_maker(mol):
	# Assumes Cube input
	N = mol.shape[0] - 1 # size - 1
	# old list is a list of the original coordinates
	def center(old_list):
		mid  = N/2
		new  = [0,0,0]
		for i in np.arange(len(old_list)):
			old = old_list[i]
			dist = abs(old-mid)
			if old == mid:
				new[i] = 0
			elif old < mid:
				new[i] = -1*dist
			else:
				new[i] = dist
		return new
	return center


center_zika = center_maker(mol)
print(center_zika([0,0,0]))





































