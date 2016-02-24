# Group 6: UPGMA Algorithm using Jukes Cantor Distance

# Authors: Josh Nixon
#          Alex Pearson
#          Bidit Acharya
#          Tracy Lou

import numpy as np
from neighbor_joining import *

"""
Makes Q matrix that decides what will be joined 
Args:
	M: symetric matrix
Returns:
	Q matrix
"""

a = np.array([[ 0, 5 , 9 , 9, 8],
			   [ 0, 0, 10, 10, 9],
			   [ 0, 0, 0,  8, 7],
			   [ 0, 0, 0,  0, 3],
			   [ 0, 0, 0,  3, 0]
			   ])

print a 
print make_Q_matrix(a)


