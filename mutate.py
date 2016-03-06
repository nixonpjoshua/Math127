# Group 6: mutate function

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou

import numpy as np
from ete3 import Tree
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

def mutate(a, t, seq):
    M = JC_matrix(a)
    M = LA.matrix_power(M, t)
    tokens = list(seq)
    for i in xrange(len(tokens)):
        p = np.zeros(4)
        p[DNA_to_position(tokens[1])] = 1
        p = np.dot(M,p)

        rand= np.random.rand()
        if rand < p[0]:
            val = 'A'
        elif rand < p[0]+p[1]:
            val = 'G'
        elif rand < p[0]+p[1]+p[2]:
            val = 'C'
        else:
            val = 'T'
        tokens[i] = val
    return ''.join(tokens)

print(mutate(.3,4,'GATTACA'))

"""
Simulates Evolution of a DNA sequence with 
    a                : alpha level (mutation rate)
    sim_time         : total time simulating
    timestep         : timestep will dictate distance between father and son node
    seq              : original input seqeuence
    selectionFn      : Decides who survives from population
    paramater        : list contating parameters to be fed into the selection fn
    expansion_factor : multiplication factor that determines number of seqs in the population. pop_size_n = 2*pop_size_n-1 - those killed by selectionFn
"""

def evolution_simulator(a, sim_time, timestep, seq,selectionFn, paramater, expansion_factor):
    t = Tree(name= seq)
    def tree_helper(curr_time, t):
        if curr_time >= sim_time - timestep:
            return
        pop = [t.name]*expansion_factor
        pop = map(lambda x: mutate(a, timestep, x), pop) # Mutates everything in the population
        pop = selectionFn(pop) # selects the survivors
        for taxa in pop:
            new_t = t.add_child(name = taxa)
            new_t.dist = timestep
            tree_helper(curr_time+timestep, new_t)
    tree_helper(0, t)
    return t

"""
Kills on average "proportion" of the population. (Uniform Distribution)

	pop       : population of sequnces that Uniform Killing decides whether or not to kill
	parameter : parameter of uniform distribution (i.e. on average how many will survive)
"""
def uniform_killing(pop, parameter):
	# filter keeps seqeunces less than proportion
	proportion = parameter[0]
    return filter(lambda y: np.random.rand() < proportion, pop)


print('tree sim')
parameter = [.4]
print(evolution_simulator(.1, 10000, 10, 'GATTACA', uniform_killing, parameter, 2))





