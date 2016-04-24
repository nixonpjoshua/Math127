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


###################################################################
#### Wrong Idea for getting the data ##############################
###################################################################

    #UNUSED Data
ptb4 ='ataagacaagcacattgtaacattagtggagcaaaatggaataatactatagaacaggta'
ptb5 ='aagacaaaattaagagaacaatttgggaataaaacaataatctttaatcactcctcagga'
ptb6 ='ggggacccagaaattgtaatgcacagttttaattgtggagggg'


## Truncated data from new Matlab File  Prof. Dynerman Sent us

d0   ='ctagcagaagaagaggtagtaattagatctgccaatttcacagacaatgctaaaatcata'
d1   ='cccaacaacaatacaagaaaaggtatacatataggaccagggagagcattttatgcaaca'
d2   ='agtagagaaaaatggaataatactttaaaccaggtagttacagaattaagggaacaattt'

ptb1 ='ctagcagaagaagagatagtaattagatctgccaatttcacagacaatgctaaaatcata'
ptb2 ='atagtacagctgaatgcatctgtagaaattaattgtacaagacccgacaacaatacaaga'
ptb3 ='aaaggtatacatataggaccagggagggcattttatgcaacaggagaaataataggagat'

ptc1 ='ctagcagaagaagaggtagtaattagatctgccgatttcacagacaatgctaaaatcata'
ptc2 ='aagaaaaggtatacatataggaccagggagagcagtttatgcaacagacagaataatagg'
ptc3 ='tagttacaaaattaagagaacaatttgtgaataaaacaataatctttactcacccctcag'

ptd1 ='ctagcagaagaagaggtagtaattagatctgcaaatttctcggacaatgctaaaaccata'
ptd2 ='caataatacaagacaaagtatacctataggaccagggaaagcagtttatgcaacaggaca'
ptd3 ='taaaagaacaatttaagaataaaacaatagtcttcaatcaatcctcaggaggggacccag'

lc1a ='ctagcagaagaagaagtagtaattagatctgaaaatttcacgaataatgctaaaatcata'
lc1b ='agtatacctatgggaccagggaaagcattttatacaacagaaataataggaaatataaga'
lc1c ='gagaacaatttaagaataaaacaatagtcttcaatcactcctcaggaggggacccagaaa'

## Patients 

ptb = [ptb1, ptb2, ptb3] 
ptc = [ptc1, ptc2, ptc3]
ptd = [ptd1, ptd2, ptd3]


## Dentist

dnt  = [d0, d1,  d2]

## Control 

ctl = [lc1a, lc1b, lc1c]







