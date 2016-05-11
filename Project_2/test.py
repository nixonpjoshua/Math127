from back_projection import back_project
from sim_image import project_fst
from MRCFile import MRCFile
from vanheel import random_rotation_matrix
import numpy as np
import sys

if len(sys.argv) != 3:
    print "please call with an mrc file as argument 1, and number of rotations as argument 2"
else:
    f = MRCFile(sys.argv[1])
    f.load_all_slices()
    images = []
    images_noisy = []
    Rs = [random_rotation_matrix() for i in xrange(int(sys.argv[2]))]
    # TODO here we can make many images later
    for R in Rs:
        image = project_fst(f.slices, R)
        images.append(image)
        noise = np.random.normal(0.0, 20*np.mean(image), image.shape)
        images_noisy.append(image + noise)
    data = zip(images, Rs)
    noisy_data = zip(images_noisy, Rs)
    f.slices = back_project(data, use_filter=False)
    f.write_file('#nofilter.mrc', overwrite=True)
    f.slices = back_project(data)
    f.write_file('filter.mrc', overwrite=True)
    f.slices = back_project(noisy_data, use_filter=False)
    f.write_file('#nofilter_noise.mrc', overwrite=True)
    f.slices = back_project(noisy_data)
    f.write_file('filter_noise.mrc', overwrite=True)
