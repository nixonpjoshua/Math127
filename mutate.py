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
The origonal idea for evolution of a sequence, creates a bifrucating tree which is fully bushy and has many nodes...
"""
def tree_simulator(a, t, seq):
    def tree_helper(elapsed, seq):
        t = Tree(name = seq)
        first = ''
        for e in xrange(elapsed, t):
            m = mutate(a, 1, seq)
            if m != seq:
                a = t.add_child(tree_helper(e, m))
                a.dist = abs(elapsed - e)
                first = m
                break
        for e in xrange(elapsed, t):
            m = mutate(a, 1, seq)
            if m != seq and m != first:
                b = t.add_child(tree_helper(e, m))
                b.dist = abs(elapsed - e)
                break
        return t
    return tree_helper(0, seq)

"""
Simulates Evolution of a DNA sequence with 
    a                : alpha level (mutation rate)
    sim_time         : total time simulating
    timestep         : timestep
    seq              : original input seqeuence
    selectionFn      : Decides who survives from population
    expansion_factor : multiplication factor that determines number of seqs in the population. pop_size_n = 2*pop_size_n-1
"""

def evolution_simulator(a, sim_time, timestep, seq, selectionFn, expansion_factor):
    t = Tree(name= seq)
    def tree_helper(curr_time, node):
        if curr_time >= sim_time - timestep:
            return
        pop = [node.name]*expansion_factor
        pop = map(lambda x: mutate(a, timestep, x), pop) 
        pop = selectionFn(pop) #the survivors
        for taxa in pop:
            new_node = node.add_child(name = taxa)
            new_node.dist = timestep
            tree_helper(curr_time+timestep, new_node)
    tree_helper(0, t)
    return t

"""
Kills 60 percent of the population.
"""
def simple_killing(pop):
    return filter(lambda y: np.random.rand() > .4, pop)


print('tree sim')
print(evolution_simulator(.1, 100, 10, 'GATTACA', simple_killing, 2))





