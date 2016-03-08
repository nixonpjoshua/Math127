import numpy as np


"""
Construct a Kimura-2 Paramter distance matrix with specified alpha, a, and beta, b.

Args:
	Given a, alpha, the rate of transitiion, and b, beta, the rate of transversions

Return:
	Transition matrix corresponding to Kimura-2 algorithm
"""

def K2_matrix(a, b):
	dg = 1-a-2*b
	M = np.array([
				 [dg, b, b, a],
				 [b, dg, a, b],
				 [b, a, dg, b],
				 [a, b, b, dg]
				 ])
	return M

"""
Compute proportion of transitions, a, and transversions, b, from two strings of the
of the same length. 

Transversions are (A<->C, A<->T, C<->G, G<->T).
Transitions are  (A<->G, C<->T)

Args: 
	s1: string 1
	s2: string 2

Return:
	Error if the strings are not of same length
	Else, return a and b of different letters

"""
def proportion_trv(s1,s2):
	num_trv = 0
	transversions = [('A','C'),('C','A'), ('A','T'), ('T','A'), ('C','G'), ('G','C'), ('G','T'), ('T','G')]
	s1     = list(s1)
	s2     = list(s2)
	length = float(len(s1))
	zipped = zip(s1,s2)
	for i in xrange(len(zipped)):
		if zipped[i] in transversions:
			num_trv +=  1
	return num_trv/length

def proportion_trs(s1,s2):
	num_trs     = 0
	transitions = [('A','G'),('G','A'), ('C','T'), ('T','C')]
	s1          = list(s1)
	s2          = list(s2)
	length      = float(len(s1))
	zipped      = zip(s1,s2)
	for i in xrange(len(zipped)):
		if zipped[i] in transitions:
			num_trs = num_trs + 1
	return num_trs/len(s1)


"""
Computes the K2 distance between the two sequences.

Args:
	s1: string 1
	s2: string 2

Returns:
	Error if the sequences are not of the same length
	Else, computes the K2 distance
"""

def K2_distance(s1, s2):
 	prop_trv = proportion_trv(s1,s2)
 	prop_trs = proportion_trs(s1,s2)
	return -0.5*(np.log(1-2*prop_trs - 2*prop_trv)) - 0.25*(np.log(1-2*prop_trv))


def kimura_matrix_maker(seqs):
  M = np.zeros((len(seqs),len(seqs)))
  for i in xrange(len(seqs) - 1):
    s1 = seqs[i]
    for j in xrange(i, len(seqs)):
      s2 = seqs[j]
      M[i][j] = K2_distance(s1,s2)
  return M



###############################
##### FOR TESTING PURPOSES ####
###############################

a = "ACT"
b = "CAT"


print(proportion_trv(a,b))


s1 = "AACTC"
s2 = "AAGTC"
s3 = "TAGTT"

seqs = [s1, s2, s3]
print(kimura_matrix_maker(seqs))

























