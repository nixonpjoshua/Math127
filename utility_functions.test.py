from distance_based_trees import *
import numpy as np

a = np.array([[0, .45, .27, .53],
				[0,   0, .40, .50],
				[0,   0,   0, .62],
				[0,   0,   0,  0]])

print sums_others(a)