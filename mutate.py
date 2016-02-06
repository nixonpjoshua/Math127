import numpy as np

def transition_matrix(a):
    b = a/3
    M = np.array([[1-a, b, b, b],
                 [b, 1-a, b, b],
                 [b, b, 1-a, b],
                 [b, b, b, 1-a]])
    return M

def DNA_to_position(nuc):
    if nuc == 'A':
        return 0
    if nuc == 'G':
        return 1
    if nuc == 'C':
        return 2
    if nuc == 'T':
        return 3

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
    mutates the DNA string seq, int t times
    according to the Jukes-Cantor model specified by the alpha value, A
"""
def mutate(a, t, seq):
    seq = list(seq)
    A = transition_matrix(a)
    for i in xrange(len(seq)):
        nuc = seq[i]
        p = np.zeros(4)
        p[DNA_to_position(nuc)] = 1
        for i in xrange(t):
            p = A.dot(p)
        seed = np.random.rand()
        if seed <= p[0]:
            seq[i] = position_to_DNA(0)
        elif seed <= p[0]+p[1]:
            seq[i] = position_to_DNA(1)
        elif seed <= p[0]+p[1]+p[2]:
            seq[i] = position_to_DNA(2)
        else:
            seq[i] = position_to_DNA(3)
    return ''.join(seq)

print(mutate(.3,4,'GATTACA'))





