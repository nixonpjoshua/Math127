# Group 6: mutate function 

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou

import numpy as np
from numpy import linalg as LA

"""
Constructs a Jukes Cantor transition Matrix with a specified alpha level a
Args:
     a: alpha level for the Jukes Cantor Matrix
Returns:
     Transition Matrix corresponding to the Jukes-Cantor Algorithm
"""

def JC_matrix(a):
    b = a/3
    M = np.array([[1-a, b, b, b],
                 [b, 1-a, b, b],
                 [b, b, 1-a, b],
                 [b, b, b, 1-a]])
    return M

"""
    Given a character 'A' 'G' 'C' or 'T' returns number corresponding 
    to the index in whcih we need to index in our vector
    Args:
        nuc: character 'A' 'G' 'C' or 'T'
    Returns:
        number 0,1,2,3 indicatin what index to set to 1 
"""

def DNA_to_position(nuc):
    if nuc == 'A':
        return 0
    if nuc == 'G':
        return 1
    if nuc == 'C':
        return 2
    if nuc == 'T':
        return 3

"""
    Given a number 0,1,2,3 returns character 'A' 'G' 'C' or 'T' 
    corresponding to the nucleic acid represented in the model
    Args:
        number 0,1,2,3 
    Returns:
        character 'A' 'G' 'C' or 'T' orresponding to the 
        nucleic acid represented in the model
"""

def position_to_DNA(pos):
    if pos == 0:
        return 'A'
    if pos == 1:
        return 'G'
    if pos == 2:
        return 'C'
    if pos == 3:
        return 'T'

"""
    Mutates the DNA string seq, int t times according to the 
    Jukes-Cantor model specified by the alpha value, a
    Args:
        a:   alpha level for the Jukes Cantor Matrix
        t:   number of time steps we plan to simulate
        seq: Sequence of DNA represented as a string
    Returns:
        simulated descendant sequence after t time steps 
        consitently represented as a string
"""

def mut(a, t, seq):
    M = JC_matrix(a)
    M = LA.matrix_power(M, t)
    def mut_helper(M,seq):
        if len(seq) == 0:
            return ''
        
        nuc                     = seq[0]
        seq                     = seq[1:]
        p                       = np.zeros(4)
        p[DNA_to_position(nuc)] = 1
        p                       = np.dot(M,p)
        rand                    = np.random.rand()

        if rand <= p[0]:
            val = 'A'
        elif rand <= p[0]+p[1]:
            val = 'G'
        elif rand <= p[0]+p[1]+p[2]:
            val = 'C'
        else:
            val = 'T'
        return val + mut_helper(M,seq)
    return mut_helper(M,seq)

print(mut(.3,4,'GATTACA'))















