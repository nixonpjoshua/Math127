import numpy as np

def estimate_3(im1,im2,im3):
	"""
	Computes Coordinate frames corresponding to three different directions
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

def get_common(im1,im2):
	"""
	Computes common line between 2 imagesas expressed by vector of values

	Args: 
    	im1: A 2d array representing image 1
    	im2: A 2d array representing image 2

    Returns:
    	A vector with array 
	"""

	return


def get_common_lines(im1,im2,im3):
	"""
	Computes common lines as expressed by array of values

	Args: 
    	im1: A 2d array representing image 1
    	im2: A 2d array representing image 2
    	im3: A 2d array representing image 3


    Returns:
    	A list with commmon lines ordered s.t.
    	[cl_(1,2) , cl_(1,3) , cl_(2,3) ]
	"""
	cl_12 = get_common(im1,im2)
	cl_13 = get_common(im1,im3)
	cl_23 = get_common(im2,im3)
	return [cl_12, cl_13, cl_23]

#####################################################################
###       Dynerman said he would give us most of the rest         ###
#####################################################################









