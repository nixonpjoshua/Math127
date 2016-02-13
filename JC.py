# Group 6: UPGMA Algorithm using Jukes Cantor Distance

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou


import matplotlib.pyplot as plt
import numpy as np

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
Constructs a random probability vector, represented as an array of size n
Args:
     n: number of entries in the probability vector
Returns:
     random probability vector
"""

def rand_vector(n):
    a = np.random.rand(n)
    return a/sum(a)

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
    # """
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
    return "Need to implement this function"

"""
Work-Horse function of UPGMA Algorithm, that given an array of sequences
returns a phylogenetic tree of distances represented as a list of lists
"""

def JC_UPGMA(M):
    return "Need to implement this function"





















