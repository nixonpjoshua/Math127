import numpy as np

###########################
"""
Construct a Kimura-3 Paramter distance matrix with specified alpha, a, and beta, b, and g, gamma.
Unlike K2, this assumes the two transversions occur at different rates

Args:
	Given a, alpha, the rate of transitiion, and b, beta, and g, gamma the rates of transversions

Return:
	Transition matrix corresponding to Kimura-3 algorithm
"""

def K3_matrix(a, b):
		
	M = np.array([1-a-2b, a, b, g],
				[a, 1-a-2b, g, b]
				[b, g, 1-a-2b, a]
				[g, b, a, 1-a-2b])
	return M

"""
Compute proportion of transitions, a, and transversion1, b, and transversion2, g, from two strings of the
of the same length. 

Transversion_b is (A<->C G<->T).
Transversion_g is (A<->T, C<->G)
Transitions are  (A<->G, C<->T)

Args: 
	s1: string 1
	s2: string 2

Return:
	Error if the strings are not of same length
	Else, return a, b, and g of different letters

"""
s1 = []
s2 = []

for i, j in zip(s1, s2):
	if len(s1) != len(s2):
        raise ValueError("Cannot compute compare DNA sequences of differing length")

count_trv1 = 0
count_trv2 = 0
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
Computes the K3 distance between the two sequences.

Args:
	s1: string 1
	s2: string 2

Returns:
	Error if the sequences are not of the same length
	Else, computes the K2 distance
"""

def K3_distance(s1, s2):