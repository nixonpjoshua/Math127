import numpy as np 

"""
How to write doctests
"""

def multiply(a,b):
	"""
	>>> multiply(4,3)
	12
	"""
	return a*b

def JC_matrix(a):

    """
    >>> np.trace(JC_matrix(.25))
    3.0
    """

    b = a/3
    M = np.array([[1-a, b, b, b],
                 [b, 1-a, b, b],
                 [b, b, 1-a, b],
                 [b, b, b, 1-a]
                 ])
    return M

foo = JC_matrix(.2)

bar = np.array([[-1,3],
				[2,-5]
				])

# ans = np.dot(bar,bar)
# print(ans)

# minimum = np.argmin(bar)
# print(minimum)

# print(bar[0,:])
"""
Computes number of times n can be divided by k
"""
def num_div(n,k):

	"""
	>>> num_div(15,2)
	7
	"""
	i = 1
	while True:
		if n < k*i:
			return i-1
		i += 1







