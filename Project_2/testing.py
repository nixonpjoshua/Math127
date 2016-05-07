from MRCFile import MRCFile
import numpy as np
f = MRCFile('zika_153.mrc')
f.load_all_slices()
f.write_file('zika_update1.mrc')
f.slices = np.zeros(f.slices.shape)
f.write_file('zika_update2.mrc')