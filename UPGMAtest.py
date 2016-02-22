from UPGMA import *
import numpy as np
a = UPGMA(np.array([[0, .45, .27, .53],[0,   0, .40, .50],[0,   0,   0, .62],[0,   0,   0,  0]]))
b = UPGMA(a)
c = UPGMA(b)
print(str(a) + '\n\n')
print(str(b) + '\n\n')
print(str(c) + '\n\n')