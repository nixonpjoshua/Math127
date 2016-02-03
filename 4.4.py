import numpy as np
a = .03, b = a/3
M = np.array([[1-a, b, b, b], [b, 1-a, b, b], [b, b, 1-a, b]. [b, b, b, 1-a]])
p = np.array([.2, .3, .4, .1])


def find_eq(M):
    D, V = np.linalg.eig(M)
    for x in xrange(D.size):
        if D[x] == 1:
            return V[:, x]


def problem_443_a(epsilon, p_t, M):
    p_eq = find_eq(M)
    def is_within_epsilon(p_t):
        t = True
        for i in xrange(p_0.size):
            t = t and abs(p_eq - p_t[i]) < epsilon
        return t
    count = 0
    while is_within_epsilon(p_t):
        p_t = M*p_t
        count += 1
    return count
small = 0
large = 0
iter = 100

for x in xrange(iter):
    large += problem_443_a(0.05, p, p_eq, M)
    small += problem_443_a(0.01, p, p_eq, M)
print('average for .05' + large/iter)
print('average for .01' + small/iter)
