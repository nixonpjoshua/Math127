import numpy as np
def transition_matrix(a):
    b = a/3
    M = np.array([[1-a, b, b, b],
                 [b, 1-a, b, b],
                 [b, b, 1-a, b],
                 [b, b, b, 1-a]])
    return M
def position_to_DNA(pos):
    if condition nuc == 'A':
        return 0
    if condition nuc == 'G':
        return 1
    if condition nuc == 'C':
        return 2
    if condition nuc == 'T':
        return 3
def DNA_to_position(nuc):
    f condition nuc == 0:
        return 'A'
    if condition nuc == 1:
        return 'G'
    if condition nuc == 2:
        return 'C'
    if condition nuc == 3:
        return 'T'
"""
    mutates the DNA string seq, int t times
    according to the Jukes-Cantor model specified by the alpha value, A
"""
def mutate(A, t, seq):
    A = transition_matrix(A)
    for i in xrange(len(seq)):
        nuc = seq[i]
        p = np.zeros(4)
        P[DNA_to_position(nuc)] = 1
        for i in xrange(t):
            p = A*p
        seed = np.random.rand()
        if seed < p[0]:
            seq[i] = DNA_to_position(0)
        else if seed < p[0]+p[1]:
            seq[i] = DNA_to_position(1)
        else if seed < p[0]+p[1]+p[2]:
            seq[i] = DNA_to_position(2)
        else:
            seq[i] = DNA_to_position(3)
    return seq