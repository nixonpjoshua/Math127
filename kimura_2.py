import numpy as np

###########################
"""
Construct a Kimura-2 Paramter distance matrix with specified alpha, a, and beta, b.

Args:
	Given a, alpha, the rate of transitiion, and b, beta, the rate of transversions

Return:
	Transition matrix corresponding to Kimura-2 algorithm
"""

def K2_matrix(a, b):
		
	M = np.array([1-a-2b, b, b, a],
				[b, 1-a-2b, a, b]
				[b, a, 1-a-2b, b]
				[a, b, b, 1-a-2b])
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
s1 = []
s2 = []

for i, j in zip(s1, s2):
	if len(s1) != len(s2):
        raise ValueError("Cannot compute compare DNA sequences of differing length")

count_trv = 0
count_trs = 0

transversions = ['AC', 'CA', 'AT', 'TA', 'CG', 'GC', 'GT', 'TG']
transitions = ['AG', 'GA', 'CT', 'TC']

for i, j in s1 and s2:
	if i + j in transversions:
		count_trv += 1
	if i + j in transitions:
		count_trs += 1

prop_diff_a = float(count_trs) / float(len(s1))
prop_diff_b = float(count_trv) / float(len(s1))

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
	return 0.5(np.log(1/(1-2*prop_diff_a - prop_diff_b))) + 0.25(np.log(1/(2*prop_diff_b)))




































