# Group 6: Neighbhor Joining Algorithm 

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
    2
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
    """
    >>> prop_diff("ATTGAC",ATGGCC) 
    float(2)/float(6)
    """
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
Finds the equilibrium point of a transition Matrix

Args:
    M: the transition Matrix for the Jukes Cantor Algorithm 

Returns:
    The equilibrium point of a transition Matrix
"""

def find_eq(M):
    D, V = np.linalg.eig(M)
    for x in xrange(D.size):
        if abs(D[x] - 1) < .0001:
            return V[:, x]/sum(V[:, x])

"""
Counts the number of time steps until the Jukes-Counter Algorithm
raches a steady state.

Args:
          M: the transition Matrix for the Jukes Cantor Algorithm 
        p_t: the initial probability vector specified by 4.4.3
    epsilon: the acceptable error bound on the equilibrium value

Returns:
    Number of iterations required for the model to converge 
    to within epsilon of the equilibrium value
"""
def counter(epsilon, p_t, M):
    p_eq = find_eq(M)
    def is_within_epsilon(p_t):
        t = True
        for i in xrange(p_t.size):
            t = t and abs(p_eq[i] - p_t[i]) < epsilon
        return t
    count = 0
    while not is_within_epsilon(p_t):
        p_t = M.dot(p_t)
        count += 1
    return count


























