import numpy as np
import sncosmo
import h5py
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed
from lsst.sims.photUtils import CosmologyObject

import time

bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
cosmo = CosmologyObject()

grid_name = 'sne_interp_models.h5'
rng = np.random.RandomState(44)

c0_range = 0.6
x1_range = 6.0
z_range = 1.2

abs_mag_mean = -19.3
abs_mag_std = 0.3

n_samples = 100000


with h5py.File(grid_name, 'r') as in_file:
    param_mins = in_file['param_mins'].value
    d_params = in_file['d_params'].value
    t_grid = in_file['t_grid'].value

    x1_vals = param_mins[0]+x1_range*rng.random_sample(n_samples)
    c0_vals = param_mins[1]+c0_range*rng.random_sample(n_samples)
    z_vals = param_mins[2]+z_range*rng.random_sample(n_samples)
    abs_mag_vals = -20.2+rng.random_sample(n_samples)*1.8
    t_vals = -20.0+rng.random_sample(n_samples)*40.0

    mag_truth = np.zeros((6,n_samples), dtype=float)
    mag_interp = np.zeros((6,n_samples), dtype=float)

    t_start = time.time()
    for ii in range(n_samples):
        x1 = x1_vals[ii]
        c0 = c0_vals[ii]
        z = z_vals[ii]
        abs_mag = abs_mag_vals[ii]
        t = t_vals[ii]

        sn = sncosmo.Model(source='salt2-extended')
        sn.set(x1=x1,c=c0,z=z)
        sn.source.set_peakmag(abs_mag+cosmo.distanceModulus(z),
                              band='bessellb', magsys='ab')

        flambda = 10.0*sn.flux(time=t, wave=10.0*bp_dict['g'].wavelen)
        ss = Sed(flambda=flambda, wavelen=bp_dict['g'].wavelen)
        mag_truth[:,ii] = bp_dict.magListForSed(ss)

        i_x1 = np.round((x1-param_mins[0])/d_params[0]).astype(int)
        i_c0 = np.round((c0-param_mins[1])/d_params[1]).astype(int)
        i_z = np.round((z-param_mins[2])/d_params[2]).astype(int)
        d_mag = abs_mag-param_mins[3]

        tag = i_x1+i_c0*100+i_z*10000
        mag_grid = in_file['%d' % tag].value
        interp_mags = np.zeros(6, dtype=float)
        for i_bp in range(6):
            mm = np.interp(t, t_grid, mag_grid[i_bp])+d_mag
            mag_interp[i_bp,ii] = mm

        if ii>0 and ii%100 == 0:
            duration = (time.time()-t_start)/3600.0
            pred = n_samples*duration/ii
            print('%d in %e hrs; predict %e' % (ii,duration,pred))

with h5py.File('sne_interp_test_data.h5', 'w') as out_file:
    out_file.create_dataset('t', data=t_vals)
    out_file.create_dataset('c0', data=c0_vals)
    out_file.create_dataset('x1', data=x1_vals)
    out_file.create_dataset('z', data=z_vals)
    out_file.create_dataset('mag_true', data=mag_truth)
    out_file.create_dataset('mag_interp', data=mag_interp)
