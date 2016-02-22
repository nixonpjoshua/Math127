# Group 6: UPGMA Algorithm using Jukes Cantor Distance

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou

import numpy as np
from ete3 import Tree

def closest_neighbors(M):
    size = len(M)
    min = 100000
    coordinates = (0,0)
    for i in xrange(size-1):
        for j in xrange(i + 1, size):
            if M[i,j] <= min:
                min = M[i,j]
                coordinates = (i,j)
    return coordinates

def update_UPGMA_matrix(M, taxa1, taxa2):
        size  = len(M)
        new_size = size - 1
        ans       = np.zeros((new_size, new_size))

        # will use the 0th row and column for the new species in the matrix
        #
        # copies over the vales from the old matrix
        computed_col = 1
        new_row = 1
        for i in xrange(size - 1):
            if i != taxa1 and i != taxa2:
                new_col = new_row + 1
                for j in xrange(i + 1, size):
                   if j != taxa1 and j != taxa2:
                       ans[new_row,new_col] = M[i,j]
                       new_col += 1
                new_row +=1
        # compute the first row of entries
        new_col = 1
        for j in xrange(size):
            if j != taxa1 and j != taxa2:
                #exploits the fact that we have an upper triangular so col > row always
                ans[0,new_col] = (M[min(taxa1, j), max(taxa1, j)] + M[min(taxa2, j), max(taxa2, j)]) / 2
                new_col += 1
        return ans
"""
Computes the subsequent Distance Matrix of the UPGMA algorithm
Args:
    M: the transition Matrix for the Jukes Cantor Algorithm
Returns:
    subsequent matrix using UPGMA
"""

def UPGMA(M, names):
    """
    >>> UPGMA(np.array([[0, .45, .27, .53],[0,   0, .40, .50],[0,   0,   0, .62],[0,   0,   0,  0]]))
    np.array([[0, .425, .575],[0,    0,  .50],[0,    0,    0]])
    """
    def search_nodes(trees ,name):
        for tree in trees:
            if tree.name == name:
                return tree
    trees = []
    while True:
        taxa1, taxa2 = closest_neighbors(M)
        if taxa1 > taxa2:
            tmp = taxa1
            taxa1 = taxa2
            taxa1 = tmp
        #define a new parent for the join
        t = Tree()
        #search for the children in trees and add them
        A = search_nodes(trees, names[taxa1])
        if A == None:
            A = t.add_child(name = names[taxa1])
        else:
            t.add_child(A)
            trees.remove(A)
        B = search_nodes(trees, names[taxa2])
        if B == None:
            B = t.add_child(name = names[taxa2])
        else:
            t.add_child(B)
            trees.remove(B)
        #delete old taxa names and update the new name
        new_names = [names[taxa1] + names[taxa2]]
        del names[taxa2]
        del names[taxa1]
        [new_names.append(name) for name in names]
        names = new_names
        #create the distance between children and parent
        avg_dist = M[taxa1, taxa2]/2
        A.dist = avg_dist
        B.dist = avg_dist
        #name the parent
        t.name = names[0]
        #add the new subtree
        trees.append(t)

        if len(M) <= 2:
            break
        M = update_UPGMA_matrix(M, taxa1, taxa2)
    return trees[0]

"""
Work-Horse function of UPGMA Algorithm, that given an array of sequences
returns a phylogenetic tree of distances represented as a list of lists
"""

def JC_UPGMA(M):
    return "Need to implement this function"





















