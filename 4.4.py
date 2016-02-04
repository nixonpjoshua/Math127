# Group 6: Neighbhor Joining Algorithm (and corresponding problems)

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou


import numpy as np

"""
Constructs a transition Matrix with a specified alpha level a
Args:
     a: alpha level for the Jukes Cantor Matrix
Returns:
     Transition Matrix corresponding to the Jukes-Cantor Algorithm
"""

def transition_matrix(a):
    b = a/3
    M = np.array([[1-a, b, b, b],
                 [b, 1-a, b, b],
                 [b, b, 1-a, b],
                 [b, b, b, 1-a]])
    return M

def rand_vector(n):
    a = np.random.rand(n)
    return a/sum(a)

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

"""
Performs simulations for all of the specific examples in problem 443

Args:
          p: probability vector

Returns:
    Number of iterations required for the model to converge 
    to within epsilon of the equilibrium value
"""
def problem_443_specific(p):
    M     = transition_matrix(.3)
    large = counter(0.05, p, M)
    small = counter(0.01, p, M)
    print('------------------------------------------------------------')
    print('number of iterations to get within epsilon = .05 is:')
    print(large)
    print('-------------------------------------------------------------')
    print('number of iterations to get with epsilon = .01 is:')
    print(small)
    print('-------------------------------------------------------------')


"""
Performs simulations for a set of random probability vectors

Returns:
    Average Number of iterations required for the model to converge 
    to within epsilon of the equilibrium value
"""
def problem_443_rand():
    small = 0
    large = 0
    iter  = 100
    M     = transition_matrix(.3)

    for x in xrange(iter):
        p      = rand_vector(4)
        large += counter(0.05, p, M)
        small += counter(0.01, p, M)
    print('------------------------------------------------------------')
    print('average number of iterations to get within epsilon = .05 is:')
    print(large/iter)
    print('-------------------------------------------------------------')
    print('average number of iterations to get within epsilon = .01 is:')
    print(small/iter)
    print('-------------------------------------------------------------')

def problem_443():
    # Solves Problem 443 by calling necessary functions
    p_a = np.array([.2, .3, .4, .1])
    p_c = np.array([.25, .25, .25, .25])
    p_d = np.array([0, 1, 0, 0])
    print('Problem 443a using probability vector')
    print(p_a)
    print('we get')
    problem_443_specific(p_a)
    print('443b using random probability vectors to obatain')
    print('average number of iterations to get within epsilon')
    problem_443_rand()
    print('Problem 443c using probability vector')
    print(p_c)
    print('we get')
    problem_443_specific(p_c)
    print('Problem 443c using probability vector')
    print(p_d)
    print('we get')
    problem_443_specific(p_d)


problem_443()





