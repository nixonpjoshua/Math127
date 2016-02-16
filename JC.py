# Group 6: UPGMA Algorithm using Jukes Cantor Distance

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou

import matplotlib.pyplot as plt
import numpy as np
import math

############# Defined for Testing ##############

A = np.array([[5, 3, 2],
              [1, 4, 5],
              [6, 7, 8],
             ])

B = np.array([[5, 3, 9],
              [3, 4, 1],
              [6, 7, 8],
             ])

################################################
"""
Constructs a Jukes Cantor transition Matrix with a specified alpha level a
Args:
     a: alpha level for the Jukes Cantor Matrix
Returns:
     Transition Matrix corresponding to the Jukes-Cantor Algorithm
"""

def JC_matrix(a):

    """
    >>> np.trace(JC_matrix(.25))
    3.0
    """

    b = a/3
    M = np.array([[1-a, b, b, b],
                 [b, 1-a, b, b],
                 [b, b, 1-a, b],
                 [b, b, b, 1-a]])
    return M


"""
Computes proportion of differing letters from two strings of the same size
Args:
    s1: string 1
    s2: string 2 
Returns:
    Throws error if the strings are not of the same length
    Else, returns proportion (in between 0 and 1) of differing letters
"""

def prop_diff(s1,s2):
    # """ LOOK curious as to why this isn't working
    # >>> prop_diff("ATTGAC","ATGGCC") 
    # float(2)/float(6)  
    # """
    if len(s1) != len(s2):
        raise ValueError("Cannot compute compare DNA sequences of differing lenth")
    diffs = 0
    i     = 0
    while i < len(s1):
        if s1[i] != s2[i]:
            diffs += 1
        i += 1
    return float(diffs)/float(len(s1))

"""
Computes the JC distance between two sequences.
Args:
    s1: string 1
    s2: string 2 
Returns:
    Throws error if the strings are not of the same length
    Else, computes JC distance
"""

def JC_distance(s1,s2):
    prop_diff = prop_diff(s1,s2)
    return 1 - (np.log(1 - 4/3*prop_diff))

"""
Computes the index of the minimum distance in distance matrix.
Args:
    M: distance matrix
Returns:
    index of the minimum entry of the Matrix
"""
def find_min(M):
    """
    >>> find_min(np.array([[0, .45, .27, .53],[0,   0, .40, .50],[0,   0,   0, .62],[0,   0,   0,  0]]))
    2
    """
    upr   = M.shape[0]*M.shape[1]
    small = 99 # Dummy value. Distances must always be < 1
    index = 0
    i     = 0
    while i < upr:
        if (M.item(i) < small) and (M.item(i) != 0):
            small = M.item(i) 
            index = i
        i += 1
    return index

"""
Computes the column number of the nearest neighbor of the matrix  
Args:
    nearest:  index of the smallest distance in the matrix
    num_cols: number of columns in the array
Returns:
    subsequent matrix using UPGMA
"""

def get_row(M):
	"""
	>>> get_row(A)
	1.0
	"""
	nearest  = find_min(M)
	row_size = M.shape[1]
	ans      = math.floor(nearest/row_size)
	return ans

"""
Computes the column number of the nearest neighbor of the matrix  
Args:
    nearest:  index of the smallest distance in the matrix
    num_cols: number of columns in the array
Returns:
    subsequent matrix using UPGMA
"""

def get_col(M):
	"""
	>>> get_col(A)
	2
	"""
	num_cols = M.shape[1]
	nearest  = find_min(M)
	col_num  = nearest % num_cols
	return  col_num

"""
Computes the subsequent Distance Matrix of the UPGMA algorithm 
Args:
    M: the transition Matrix for the Jukes Cantor Algorithm 
Returns:
    subsequent matrix using UPGMA
"""

def new_dist(M):
    """
    >>> new_dist(np.array([[0, .45, .27, .53],[0,   0, .40, .50],[0,   0,   0, .62],[0,   0,   0,  0]]))
    np.array([[0, .425, .575],[0,    0,  .50],[0,    0,    0]])
    """
    min_row   = get_row(M)       # equiv. to  species 1
    min_col   = get_col(M)       # equiv. to  species 2
    num_rows  = M.shape[0] 
    num_cols  = M.shape[1]
    new_spec  = min(min_row,min_col)
    ans       = np.zeros((num_rows - 1, num_cols - 1))
    
    # will use new spec as the row of my "new species" in the new matrix

    for i in xrange(num_rows):
    	for j in xrange(num_cols):
    		if i >= j: #only analyzing upper triangle of entries
    			pass
    		elif i == min_row and j == min_col:
    			pass
    		else:
    			ans[i][j] = (M[][] + M[][])/2
    		else:
    			ans[i][j] = M[i][j]
  
    return ans

"""
Work-Horse function of UPGMA Algorithm, that given an array of sequences
returns a phylogenetic tree of distances represented as a list of lists
"""

def JC_UPGMA(M):
    return "Need to implement this function"





















