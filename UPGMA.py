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


def closest_neighbors(M):
    size = len(M)
    min = 100000
    coordinates = (0,0)
    for i in xrange(size-1):
        for j in xrange(i + 1, size):
            if M[i,j] <= min:
                min = M[i,j]
                coordinates = (i,j)
    return coordinates


"""
Computes the subsequent Distance Matrix of the UPGMA algorithm
Args:
    M: the transition Matrix for the Jukes Cantor Algorithm
Returns:
    subsequent matrix using UPGMA
"""

def UPGMA(M):
    """
    >>> UPGMA(np.array([[0, .45, .27, .53],[0,   0, .40, .50],[0,   0,   0, .62],[0,   0,   0,  0]]))
    np.array([[0, .425, .575],[0,    0,  .50],[0,    0,    0]])
    """
    taxa1, taxa2 = closest_neighbors(M)
    size  = len(M)
    new_size = size - 1
    ans       = np.zeros((new_size, new_size))

    # will use the 0th row and column for the new species in the matrix
    #
    # copies over the vales from the old matrix
    computed_col = 1
    new_row = 1
    for i in xrange(size - 1):
        if i != taxa1 and i != taxa2:
            new_col = new_row + 1
    	    for j in xrange(i + 1, size):
               if j != taxa1 and j != taxa2:
                   ans[new_row,new_col] = M[i,j]
                   new_col += 1
            new_row +=1
    # compute the first row of entries
    new_col = 1
    for j in xrange(size):
        if j != taxa1 and j != taxa2:
            #exploits the fact that we have an upper triangular so col > row always
            ans[0,new_col] = (M[min(taxa1, j), max(taxa1, j)] + M[min(taxa2, j), max(taxa2, j)]) / 2
            new_col += 1
    return ans

"""
Work-Horse function of UPGMA Algorithm, that given an array of sequences
returns a phylogenetic tree of distances represented as a list of lists
"""

def JC_UPGMA(M):
    return "Need to implement this function"





















