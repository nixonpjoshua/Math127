# Group 6: UPGMA Algorithm using Jukes Cantor Distance

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou

import matplotlib.pyplot as plt
import numpy as np
import math
from ete3 import Tree

"""
Makes Q matrix that decides what will be joined 
Args:
	M: symetric matrix
Returns:
	Q matrix
"""

def make_Q_matrix(M):
	N  = M.shape[0]        # number of taxa
	Q  = np.zeros(M.shape) # matrix to be returned 
	for i in xrange(N):
		for j in xrange(N):
			if i < j:
				Q[i][j] = (N-2)*M[i][j]
				for k in xrange(N):
						Q[i][j] = Q[i][j] - M[i][k] - M[j][k]
			else:
				pass
