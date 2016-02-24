# Group 6: UPGMA Algorithm using Jukes Cantor Distance

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou

import numpy as np
from ete3 import Tree

def sums_others(M):
    size = len(M)
    sums = np.zeros(size)
    for i in xrange(size):
        s = 0
        for other in xrange(size):
            s += M[min(other, i), max(other, i)]
        sums[i] = s
    return sums

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

"""
Makes Q matrix that decides what will be joined 
Args:
    M: symetric matrix
Returns:
    Q matrix
"""
def make_Q_matrix(M, sums):
    N  = M.shape[0]        # number of taxa
    Q  = np.zeros(M.shape) # matrix to be returned 
    for i in xrange(N):
        for j in xrange(N):
            if i < j:
                Q[i][j] = (N-2)*M[i][j] - sums[i] - sums[j]
    return Q

"""
Computes the subsequent Distance Matrix of the UPGMA algorithm
Args:
    M: the transition Matrix for the Jukes Cantor Algorithm
Returns:
    subsequent matrix using UPGMA
"""
def UPGMA_new_dist(M, taxa1, taxa2, j):
    return (M[min(taxa1, j), max(taxa1, j)] + M[min(taxa2, j), max(taxa2, j)]) / 2

def neighbor_joining_new_dist(M, taxa1, taxa2, j):
    return (M[min(taxa1, j), max(taxa1, j)] + M[min(taxa2, j), max(taxa2, j)] - M[taxa1, taxa2])/2.0

def split_dist(M, taxa1, taxa2):
    avg_dist = M[taxa1, taxa2]/2
    return (avg_dist, avg_dist)

# TODO this means we compute sums others twice per loop
def neighbor_joining_parent_dist(M, taxa1, taxa2):
    sums = sums_others(M)
    taxa1_dist = M[taxa1, taxa2]/2.0 + (sums[taxa1] - sums[taxa2])/(2*(len(M)-2))
    return (taxa1_dist, M[taxa1, taxa2] - taxa1_dist)
    

def update_matrix(M, taxa1, taxa2, new_dist_fn):
        size  = len(M)
        new_size = size - 1
        ans       = np.zeros((new_size, new_size))

        # will use the 0th row and column for the new species in the matrix
        #
        # copies over the vales from the old matrix
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
                ans[0,new_col] = UPGMA_new_dist(M, taxa1, taxa2, j)
                new_col += 1
        return ans

def neighbor_based_method(M, names, closest_neighbors_fn new_dist_fn, parent_dist_fn):
    def search_nodes(trees ,name):
        for tree in trees:
            if tree.name == name:
                return tree
    trees = []
    while True:
        taxa1, taxa2 = closest_neighbors_fn(M)
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
        A.dist, B.dist = parent_dist_fn(M, taxa1, taxa2)
        #name the parent
        t.name = names[0]
        #add the new subtree
        trees.append(t)

        if len(M) <= 2:
            break
        M = update_matrix(M, taxa1, taxa2, new_dist_fn)
    return trees[0]

def UPGMA(M, names):
    return neighbor_based_method(M, names, closest_neighbors UPGMA_new_dist, split_dist)

def neighbor_joining(M, names):
    return neighbor_based_method(M, names, lambda m : closest_neighbors(make_Q_matrix(m)), neighbor_joining_new_dist, neighbor_joining_parent_dist)

















