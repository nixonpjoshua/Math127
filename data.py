from scipy.io import loadmat

from mutate import *


hiv_data = loadmat('flhivdata.mat')

#UNUSED Data
ptb4 ='ataagacaagcacattgtaacattagtggagcaaaatggaataatactatagaacaggta'
ptb5 ='aagacaaaattaagagaacaatttgggaataaaacaataatctttaatcactcctcagga'
ptb6 ='ggggacccagaaattgtaatgcacagttttaattgtggagggg'


print(hiv_data.keys())

## Truncated data from new Matlab File  Prof. Dynerman Sentus

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







