import numpy as np 

"""
1.  Compute 3D DFT of molecule mol, producing F(mol) another 3D Array
2.  Restrict F(mol) to the image plane P_f thus producing a 2D Array
2a. Interpolation buisness that I don't quite understand how to code up


"""


def project_fst(mol, R):
	nice = np.fft.fftn(mol)
	def transform(x, y):
		# use R to convert to i, j, k from x, y
	def interpolate(i, j, k):
		# use nice to interpolate and return values
	output = #2d array
	for x xrange(len(output)):
		for y xrange(len(output)):
			output[x][y] = interpolate(transform(x, y))
	return output


"""
Questions:

What type of Discrete Fourier Transform of mol are we computing?
Since we are given a 3D array mol I am pretty sure we are computing
3 dimenstional fourier transform so I looked at:

http://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftn.html#numpy.fft.fftn

"""

def compute_3D_DFT(mol):
	return










































