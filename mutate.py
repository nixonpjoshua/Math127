import numpy as np

def transition_matrix(a):
    b = a/3
    M = np.array([[1-a, b, b, b],
                 [b, 1-a, b, b],
                 [b, b, 1-a, b],
                 [b, b, b, 1-a]])
    return M
<<<<<<< Updated upstream
def position_to_DNA(nuc):
    if nuc == 'A':
        return 0
    if nuc == 'G':
        return 1
    if nuc == 'C':
        return 2
    if nuc == 'T':
        return 3
def DNA_to_position(pos):
    if pos == 0:
        return 'A'
    if pos == 1:
        return 'G'
    if pos == 2:
        return 'C'
    if pos == 3:
=======

def position_to_DNA(nuc):
    if  nuc == 'A':
        return 0
    if  nuc == 'G':
        return 1
    if  nuc == 'C':
        return 2
    if  nuc == 'T':
        return 3

def DNA_to_position(nuc):
    if  nuc == 0:
        return 'A'
    if  nuc == 1:
        return 'G'
    if  nuc == 2:
        return 'C'
    if  nuc == 3:
>>>>>>> Stashed changes
        return 'T'
"""
    mutates the DNA string seq, int t times
    according to the Jukes-Cantor model specified by the alpha value, A
"""
def mutate(a, t, seq):
    A = transition_matrix(a)
    for i in xrange(len(seq)):
        nuc = seq[i]
        p = np.zeros(4)
        p[DNA_to_position(nuc)] = 1
        p = np.transpose(p)
        print(p)
        for i in xrange(t):
            p = A*p
        seed = np.random.rand()
        if seed <= p[0]:
            seq[i] = DNA_to_position(0)
        elif seed <= p[0]+p[1]:
            seq[i] = DNA_to_position(1)
        elif seed <= p[0]+p[1]+p[2]:
            seq[i] = DNA_to_position(2)
        else:
            seq[i] = DNA_to_position(3)
    return seq

print(mutate(.3,4,'GATTACA'))





