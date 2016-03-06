"""
The original idea for evolution of a sequence, creates a bifrucating tree which is fully bushy and has many nodes...
"""
def tree_simulator(a, t, seq):
    def tree_helper(elapsed, seq):
        t = Tree(name = seq)
        first = ''
        for e in xrange(elapsed, t):
            m = mutate(a, 1, seq)
            if m != seq:
                a      = t.add_child(tree_helper(e, m))
                a.dist = abs(elapsed - e)
                first  = m
                break
        for e in xrange(elapsed, t):
            m = mutate(a, 1, seq)
            if m != seq and m != first:
                b      = t.add_child(tree_helper(e, m))
                b.dist = abs(elapsed - e)
                break
        return t
    return tree_helper(0, seq)