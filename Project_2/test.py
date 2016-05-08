from orientations_from_images import reconstruct_orientations
from back_projection import back_project
from sim_image import project_fst
from MRCFile import MRCFile
import numpy as np
import sys
def rotate_x(n):
    thetas = np.linspace(0, np.pi, n)
    matrices = []
    for j in np.arange(len(thetas)):
        i = thetas[j]
        new_matrix = np.array([[ 1, 0, 0],
                               [ 0 , np.cos(i), -np.sin(i)],
                               [ 0 , np.sin(i), np.cos(i)]
                               ])
        matrices.append(new_matrix)
    return np.array(matrices)
if len(sys.argv) != 2:
    print "please call with an mrc file as argument"
else:
    f = MRCFile(sys.argv[1])
    f.load_all_slices()
    images = []
    Rs = rotate_x(4)
    # TODO here we can make many images later
    for R in Rs:
        images.append(project_fst(f.slices, R))
    # orientations = reconstruct_orientations(images, num_lines=1, granularity=1)
    # print "Rs\n"
    # print Rs
    # print "\norientations\n"
    # print orientations
    f.slices = back_project(images, Rs)
    f.write_file('true_R.mrc')
    # f.slices = back_project(images, orientations)
    # f.write_file('estimated_R.mrc')
