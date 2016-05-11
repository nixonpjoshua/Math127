from MRCFile import MRCFile
from sim_image import project_fst
import matplotlib.pyplot as plt
from vanheel import random_rotation_matrix
zika_153 = MRCFile('zika_153.mrc')
zika_153.load_all_slices()
mol = zika_153.slices
# image = project_fst(mol, np.eye(3))
image = project_fst(mol, random_rotation_matrix())
plt.imshow(image)
plt.show()






































