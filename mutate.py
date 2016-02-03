import numpy as np
def mutate(A, t, seq):
    if t == 0:
        return seq
    A.vector = A.matrix*A.vector
    return mutate(A, t - 1, seq)