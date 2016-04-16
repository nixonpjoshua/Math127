from ete3 import Tree 
import numpy as np

t = Tree(name = "Lisa")
t.add_child(name = "Alex")
t.add_child(name = "Dilly")

t = Tree(name = "GATTACA") 

def uniform_killing(pop, proportion):
    return filter(lambda y: np.random.rand() < proportion, pop) # filter keeps values less than proportion

pop = [t.name]*10

for i in xrange(0,5):
	print(uniform_killing(pop,.9))
























































