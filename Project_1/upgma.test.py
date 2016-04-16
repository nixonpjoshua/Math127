from UPGMA import *
from ete3 import Tree
import numpy as np
a = UPGMA(np.array([[0, .45, .27, .53],
					[0,   0, .40, .50],
					[0,   0,   0, .62],
					[0,   0,   0,  0]]), 
					['a', 'b', 'c', 'd'])
print(a)

# http://www.southampton.ac.uk/~re1u06/teaching/upgma/upgma15.PNG

b = UPGMA(np.array([[0, 19.0,  27.0,  8.0, 33.0, 18.0 ,13.0],
					[0,    0,  31.0, 18.0, 36.0,  1.0, 13.0],
					[0,    0,   0,   26.0, 41.0, 32.0, 29.0],
					[0,    0,   0,      0, 31.0, 17.0, 14.0],
					[0,    0,   0,      0,    0, 35.0, 28.0],
					[0,    0,   0,      0,    0,    0, 12.0]]), 
					['a', 'b', 'c', 'd','e','f','g'])

print(b)