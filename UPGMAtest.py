from UPGMA import *
import numpy as np
a = UPGMA(np.array([[0, .45, .27, .53],[0,   0, .40, .50],[0,   0,   0, .62],[0,   0,   0,  0]]), ['a', 'b', 'c', 'd'])
print(str(a) + '\n\n')