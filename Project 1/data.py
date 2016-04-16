from scipy.io import loadmat
from jukes_cantor import *
from distance_based_trees import *
from mutate import *


hiv_data = loadmat('flhivdata.mat')

print(hiv_data.keys())

dnt    = hiv_data["dnt"][0]

ctrl_1 = hiv_data["lc1"][0]
ctrl_5 = hiv_data["lc5"][0]

ptb    = hiv_data["ptb"][0]
ptc	   = hiv_data["ptc"][0]
ptd    = hiv_data["ptd"][0]

min_len = min(len(dnt),len(ctrl_1),len(ctrl_5),len(ptb),len(ptc),len(ptd))
print(min_len)

def chop(seq):
	ans    = min_len*["o"]
	tokens = list(seq)
	for i in xrange(min_len):
		ans[i] = tokens[i]
	return ''.join(ans)

dnt    = chop(dnt)
ctrl_1 = chop(ctrl_1)
ctrl_5 = chop(ctrl_5)

ptb    = chop(ptb)
ptc	   = chop(ptc)
ptd    = chop(ptd)

names = ["dnt", "ptb", "ptc", "ptd", "ctrl1", "ctrl5"]
seqs  = [dnt, ptb, ptc, ptd, ctrl_1, ctrl_5]

M = JC_matrix_maker(seqs)

UPGMA(M, names)

names = ["dnt", "ptb", "ptc", "ptd", "ctrl1", "ctrl5"]
seqs  = [dnt, ptb, ptc, ptd, ctrl_1, ctrl_5]
M = JC_matrix_maker(seqs)
tree = neighbor_joining(M, names)
names = ["dnt", "ptb", "ptc", "ptd", "ctrl1", "ctrl5"]

M = np.array([[ 0        ,  0.04615385,  0.04923077,  0.12923077,  0.40307692, 0.11692308],
       		  [ 0        ,  0         ,  0.06461538,  0.13230769,  0.39384615, 0.12      ],
              [ 0        ,  0         ,  0         ,  0.13538462,  0.40923077, 0.13538462],
       		  [ 0        ,  0         ,  0         ,  0         ,  0.4       , 0.11384615],
       		  [ 0        ,  0         ,  0         ,  0         ,  0         , 0.38461538],
              [ 0        ,  0         ,  0         ,  0         ,  0         , 0        ]])

tree = neighbor_joining(M, names)
tree.show()
