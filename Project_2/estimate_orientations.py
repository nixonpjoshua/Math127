import numpy as np
from scipy import interpolate

def make_interpolator(img): 
	"""
    Creates an interpolator object for a given img
    
    Args:
        img: An NxN inarray
    Returns:
        Interpolator for the corresponding image
    """
	size = img.shape[0]
	x = np.linspace(-1,1,size) 
	y = np.linspace(-1,1,size)
	im_interpolater = interpolate.RegularGridInterpolator((x,y),img,method = "linear", bounds_error = False, fill_value = 0)
	return im_interpolater

def sample_line(theta,gran):
    """
    Computes a sampled line in viewing sqaure through origin with angle counterclockwise from the horizon equal to theta, and granularity of sampling equal to n
    
    Args:
        gran: granularity level  
        theta: Angle from horizon that determines the line 
    Returns:
        Line in square through origin with angle counteclockwise from the horizon equal to theta
    """
    if np.isclose([theta],[0]) or np.isclose([theta],[np.pi/2]) or np.isclose([theta],[np.pi]) or np.isclose([theta],[3*np.pi/2]) or  np.isclose([theta],[2*np.pi]):
        rs = np.linspace(-1,1,gran) 
    else:
        max_r = min( abs(1./np.cos(theta)) , abs(1./np.sin(theta)))
        rs = np.linspace(-max_r,max_r,gran) 

    xs = np.cos(theta)*rs
    ys = np.sin(theta)*rs

    return zip(xs,ys)

def sample_lines(n,gran):
    """
    Creates n evenly spaced (w.r.t. theta) sample lines for the square
    
    Args:
        n: number of sample lines  
        gran: granularity level 
    Returns:
        n sample lines as repreented by a list of xs and ys
    """
    ans = [0]*n
    thetas = np.linspace(0,np.pi,n) 
    for i in xrange(n):
        theta  = thetas[i]
        ans[i] =  sample_line(theta,gran)
    return ans

def im1_vector(list_coords):
	"""
    Computes interpolation vector for image 1 over a specified line
    
    Args:
        list_coords: A list of coordinates represented as tuples (i.e. samplings of a line)
    Returns:
    	vector represented by an np array
    """
	ans = np.zeros(len(list_coords))
	for i in xrange(len(list_coords)):
	    coords = list_coords[i]
	    ans[i] = im1_interpolater(list(coords))[0]
	return ans

def im2_vector(list_coords):
	"""
    Computes interpolation vector for image 2 over a specified line
    
    Args:
        list_coords: A list of coordinates represented as tuples (i.e. samplings of a line)
    Returns:
    	vector represented by an np array
    """
	ans = np.zeros(len(list_coords))
	for i in xrange(len(list_coords)):
	    coords = list_coords[i]
	    ans[i] = im2_interpolater(list(coords))[0]
	return ans

def find_common_line(lines):
	"""
	Finds common line over the space of lines 

	Args:
	    lines: a list of lines where each entry contains a list of tuples that represents a line
	Returns:
			indeces of common lines
	"""
	max_IP      = -1  # will be over written because IP must always be positive
	max_index_1 = .2  # If not overwritten then error will be thrown because cannot
	                # index with decimal
	max_index_2 = .2
	# Meat of the argument go through all n**2 line combos 
	for i in xrange(len(lines)):
	    line_1 = lines[i]
	    for j in xrange(len(lines)):
	        line_2 = lines[j]
	        IP = np.vdot(im1_vector(line_1),im2_vector(line_2))
	        if  IP > max_IP:
	            max_index_1 = i
	            max_index_2 = j
	            max_IP = IP
	return [max_index_1, max_index_2]

def find_xy(theta):
	"""
	Computes unit vector xy coordinates given a line through the origin as represented by theta
	Args:
		im1: A NxN array that represents an image
		im2: A NxN array that represents
	    num_lines: a list of lines where each entry contains a list of tuples that represents a line
	    gran: number of sampling points for each line
	Returns:
		A list of x,y coordinates multiplied by necessary sign component
	"""
    # Get correct directionality
    if theta <= np.pi:
        sign = 1
    else:
        sign = -1
    # Get Unit vector i.e solve for x when r = 1
    x = sign*(np.cos(theta))
    y = sign*(np.sin(theta))
    return [x,y]    

def common_line(im1,im2, num_lines, gran):
	"""
	Computes common line over the space of lines 
	Args:
		im1: A NxN array that represents an image
		im2: A NxN array that represents
	    num_lines: a list of lines where each entry contains a list of tuples that represents a line
	    gran: number of sampling points for each line
	Returns:
		A list of x,y coordinates multiplied by necessary sign component
	"""
	#TODO abstraction barrier violation make sure it works need thetas because each line is uniquely represented by a theta
	thetas = np.linspace(0,np.pi,num_lines) 

	im1_interpolater = make_interpolator(im1)
	im2_interpolater = make_interpolator(im2)
	
	lines = sample_lines(num_lines,gran)

	com_line = find_common_line(lines)

	line_1_index = com_line[0] #common line 1 index
	line_2_index = com_line[1] #common line 2 index

	theta_1 = thetas[line_1_index] 
	theta_2 = thetas[line_2_index]

	com_1   = find_xy(theta_1)
	com_2   = find_xy(theta_2)

	return [com_1,com_2] #currently retuning common line as list of thetas corresponding to the images

def common_lines(Images, num_lines, gran):
	"""
	Finds all common lines for space of images
    
    Args:
    	Images: A list of NxN arrays that represent images
        num_lines: a list of lines where each entry contains a list of tuples that represents a line
        gran: number of sampling points for each line
    Returns:
    	A matrix "L" of unit vectors where L_ij is the unit vector for the line in image I that represents the "commonality between i and j.
	"""
	length = len(images)
	L = np.zeros((length, length, 2))
	for i in xrange(length):
		for j in xrange(length):
			L[i,j] = np.array(common_line(Images[i],Images[j],num_lines,gran)
	return L

im1 = np.zeros((9, 9))
im2 = np.zeros((9, 9))
im1[4] = 1
im2[:,4] = 1
# im1 = np.array( [ [1,2,3,4,5,6,7,8,9,1],
#                   [1,2,3,4,5,6,7,8,1,10],
#                   [1,2,3,4,5,6,7,1,9,10],
#                   [1,2,3,4,5,6,1,8,9,10],
#                   [1,2,3,4,5,1,7,8,9,10],
#                   [1,2,3,4,1,6,7,8,9,10],
#                   [1,2,3,1,5,6,7,8,9,10],
#                   [1,2,1,4,5,6,7,8,9,10],
#                   [1,1,3,4,5,6,7,8,9,10],
#                   [1,2,3,4,5,6,7,8,9,10]
#                 ])


# im2 = np.array( [[10,9,8,7,6,5,4,3,2,1],
#                   [9,8,7,6,5,4,3,2,1,1],
#                   [9,8,7,6,5,4,3,1,2,1],
#                   [9,8,7,6,5,4,1,3,2,1],
#                   [9,8,7,6,5,1,4,3,2,1],
#                   [9,8,7,6,1,5,4,3,2,1],
#                   [9,8,7,1,6,5,4,3,2,1],
#                   [9,8,1,7,6,5,4,3,2,1],
#                   [9,1,8,7,6,5,4,3,2,1],
#                   [1,2,3,4,5,6,7,8,9,10]
#                 ])



im1_interpolater = make_interpolator(im1)
im2_interpolater = make_interpolator(im2)
	
print(common_line(im1,im2, num_lines = 100, gran = 100))

	













