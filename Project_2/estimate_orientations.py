import numpy as np

def estimate_3(im1,im2,im3):
	"""
	Gives Coordinate frames corresponding to three different directions
	used as a helper function for estimate_orientations

	Args: 
    	im1: A 2d array representing image 1
    	im2: A 2d array representing image 2
    	im3: A 2d array representing image 3

    Returns:
    	Three 3*3 matrices corresponding to the three
    	different viewing directions for each image
	"""
	return

#####################################################################
# Meat of argument goes here i.e. how to find rotation matrices #####
#####################################################################

def compute_common_lines(im1,im2,im3):
	"""
	Gives Coordinate frames corresponding to three different directions
	used as a helper function for estimate_orientations

	Args: 
    	im1: A 2d array representing image 1
    	im2: A 2d array representing image 2
    	im3: A 2d array representing image 3

    Returns:
    	# TODO: Figure out what data structure lines should be 
    	# represented as does it have to be values in the array?
	"""
	return

def compute_common_lines(im1,im2):
	"""
	Gives Coordinate frames corresponding to three different directions
	used as a helper function for estimate_orientations

	Args: 
    	im1: A 2d array representing image 1
    	im2: A 2d array representing image 2

    Returns:
    	# TODO: Figure out what data structure lines should be 
    	# represented as does it have to be values in the array?
	"""
	return


# Leave all of the trig to Prof. Dynerman

#
# Am confused from here on as to how to "set" the first rotaion frame at 
# F_1 = np.array([[1,0,0],
#				  [0,1,0],
#				  [0,0,1]
#				])
#

#####################################################################
# Uses estimate_3 successively to find rotation matrices and then ###
# usesback projection to try and find rho                         ###
#####################################################################

def estimate_orientations():
	"""
	"""
	return 















