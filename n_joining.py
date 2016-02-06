# Group 6: Neighbhor Joining Algorithm (and corresponding problems)

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou


import numpy as np

def find_eq(M):

# Finds the equilibrium point of a transition Matrix

# Args:
#     M: the transition Matrix for the Jukes Cantor Algorithm 

# Returns:
#     The equilibrium point of a transition Matrix

    D, V = np.linalg.eig(M)
    for x in xrange(D.size):
        if abs(D[x] - 1) < .0001:
            return V[:, x]/sum(V[:, x])


def counter(epsilon, p_t, M):

# Counts the number of time steps until the Jukes-Counter Algorithm
# raches a steady state.

# Args:
#           M: the transition Matrix for the Jukes Cantor Algorithm 
#         p_t: the initial probability vector specified by 4.4.3
#     epsilon: the acceptable error bound on the equilibrium value

# Returns:
#     Number of iterations required for the model to converge 
#     to within epsilon of the equilibrium value

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

def problem_443_ex():
    a   = .03 # alpha level specified in the example case of problem 4.43
    b   = a/3
    M   = np.array([[1-a, b, b, b], 
                    [b, 1-a, b, b], 
                    [b, b, 1-a, b], 
                    [b, b, b, 1-a]
                   ])
    p = np.array([.2, .3, .4, .1])
    small = 0
    large = 0
    iter  = 1

    for x in xrange(iter):
        large += counter(0.05, p, M)
        small += counter(0.01, p, M)
    print('------------------------------------------------------------')
    print('average number of iterations to get with epsilon = .05 is:')
    print(large/iter)
    print('-------------------------------------------------------------')
    print('average number of iterations to get with epsilon = .01 is:')
    print(small/iter)
    print('-------------------------------------------------------------')

problem_443_ex()


