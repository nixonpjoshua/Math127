import numpy as np

A = np.array([[5, 3, 2],
              [1, 4, 5],
              [6, 7, 8],
             ])

B = np.array([[5, 3, 9],
              [3, 4, 1],
              [6, 7, 8],
             ])

################################################
"""
Constructs a Jukes Cantor transition Matrix with a specified alpha level a
Args:
     a: alpha level for the Jukes Cantor Matrix
Returns:
     Transition Matrix corresponding to the Jukes-Cantor Algorithm
"""

def JC_matrix(a):

    """
    >>> np.trace(JC_matrix(.25))
    3.0
    """

    b = a/3
    M = np.array([[1-a, b, b, b],
                 [b, 1-a, b, b],
                 [b, b, 1-a, b],
                 [b, b, b, 1-a]])
    return M


"""
Computes proportion of differing letters from two strings of the same size
Args:
    s1: string 1
    s2: string 2 
Returns:
    Throws error if the strings are not of the same length
    Else, returns proportion (in between 0 and 1) of differing letters
"""

def prop_diff(s1,s2):
    # """ LOOK curious as to why this isn't working
    # >>> prop_diff("ATTGAC","ATGGCC") 
    # float(2)/float(6)  
    # """
    if len(s1) != len(s2):
        raise ValueError("Cannot compute compare DNA sequences of differing length")
    diffs = 0
    i     = 0
    while i < len(s1):
        if s1[i] != s2[i]:
            diffs += 1
        i += 1
    return float(diffs)/float(len(s1))

"""
Computes the JC distance between two sequences.
Args:
    s1: string 1
    s2: string 2 
Returns:
    Throws error if the strings are not of the same length
    Else, computes JC distance
"""

def JC_distance(s1,s2):
    prop_diffs = prop_diff(s1,s2)
    return 1 - (np.log(1 - 4/3*prop_diffs))

"""
Returns JC Matrix give sequences
"""

def JC_matrix_maker(seqs):
  M = np.zeros((len(seqs),len(seqs)))
  for i in xrange(len(seqs) - 1):
    s1 = seqs[i]
    for j in xrange(i, len(seqs)):
      s2 = seqs[j]
      M[i][j] = JC_distance(s1,s2)
  return M


# human = "aactc"
# chimp = "aagtc"
# orang = "tagtt"
# seqs = [human, chimp, orang]
# print(JC_matrix_maker(seqs))

# Works according to:
# http://homes.cs.washington.edu/~ruzzo/courses/gs559/09wi/lectures/7A_distance.pdf











