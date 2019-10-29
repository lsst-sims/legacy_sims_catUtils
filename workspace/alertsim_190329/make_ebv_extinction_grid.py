import numpy as np
import os
import h5py
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed
from lsst.sims.utils import defaultSpecMap

sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']
sed_name = 'sed_flat.txt'

lsst_bp = BandpassDict.loadTotalBandpassesFromFiles()

ebv_grid_1 = np.arange(0.01, 8.0, 0.01)
ebv_grid_2 = np.arange(9.0, 120.0, 1.0)
ebv_grid = np.concatenate([ebv_grid_1, ebv_grid_2])
ext_grid = np.zeros((6,len(ebv_grid)), dtype=float)

unextincted_sed = Sed()
unextincted_sed.readSED_flambda(os.path.join(sed_dir, defaultSpecMap[sed_name]))
unextincted_mags = lsst_bp.magListForSed(unextincted_sed)

a_x, b_x = unextincted_sed.setupCCM_ab()

for i_ebv, ebv in enumerate(ebv_grid):
    extincted_sed = Sed(wavelen=unextincted_sed.wavelen,
                        flambda=unextincted_sed.flambda)

    extincted_sed.addDust(a_x, b_x, R_v=3.1, ebv=ebv)
    mags = lsst_bp.magListForSed(extincted_sed)
    for i_bp, bp in enumerate('ugrizy'):
        ext_grid[i_bp][i_ebv] = mags[i_bp]-unextincted_mags[i_bp]

assert ext_grid.min()>0.0

with h5py.File('data/ebv_grid.h5', 'w') as out_file:
    out_file.create_dataset('ebv_grid', data=ebv_grid)
    out_file.create_dataset('extinction_grid', data=ext_grid)
