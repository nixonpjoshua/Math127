import numpy as np 

"""
Goal

1.  Compute 3D DFT of molecule mol, producing F(mol) another 3D Array
2.  Restrict F(mol) to the image plane P_f thus producing a 2D Array
2a. Interpolation buisness that I don't quite understand how to code up
"""

## Alex's Approach
def alex_project_fst(mol, R):
	nice        = np.fft.fftn(mol) 
	restriction = restriction(mol,R)
	return np.fft.ifftn(np.fft.fftn(mol))

"""
Creates a 2D Array that represents a restriction of mol.
Args: 
	mol: an NxNXN matrix that contains voxels of electro potential
	R: 3x3 rotation matrix 
Returns:
	Returns a restriction of mol onto the image plane
Note:
	Interpolation is required for this step
"""
def restriction(mol,R):
	return #2Darray


"""
Creates function that transforms old coordinates into 0 centered
coordinates.
Args: 
	mol: an NxNXN matrix that contains voxels of electro potential
Returns:
	Function that given a list of 3 corrdinates transforms them into 
	coordinate system centered at 0.

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




"""

"""

## Josh's Approach
def josh_project_fst(mol, R):
	nice = np.fft.fftn(mol)
	def transform(x, y):
		# use R to convert to i, j, k from x, y
	def interpolate(i, j, k):
		# use nice to interpolate and return values
	output = #2d array
	for x xrange(len(output)):
		for y xrange(len(output)):
			output[x][y] = interpolate(transform(x, y))
	return output













































