import numpy as np
import sncosmo
import h5py
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed
from lsst.sims.photUtils import CosmologyObject

bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
cosmo = CosmologyObject()

grid_name = 'sne_interp_models.h5'
rng = np.random.RandomState(44)

c0_range = 0.6
x1_range = 6.0
z_range = 1.2

abs_mag_mean = -19.3
abs_mag_std = 0.1

with h5py.File(grid_name, 'r') as in_file:
    param_mins = in_file['param_mins'].value
    d_params = in_file['d_params'].value
    t_grid = in_file['t_grid'].value

    for ii in range(10):
        x1 = param_mins[0]+rng.random_sample()*x1_range
        c0 = param_mins[1]+rng.random_sample()*c0_range
        z = param_mins[2]+rng.random_sample()*z_range
        abs_mag = -20.2+rng.random_sample()*1.8
        t = rng.uniform(0.0, 5.0)
        print('\nparams %e %e %e %e -- t %e' % (x1,c0,z,abs_mag,t))
        
        sn = sncosmo.Model(source='salt2-extended')
        sn.set(x1=x1,c=c0,z=z)
        sn.source.set_peakmag(abs_mag+cosmo.distanceModulus(z),
                              band='bessellb', magsys='ab')

        flambda = 10.0*sn.flux(time=t, wave=10.0*bp_dict['g'].wavelen)
        ss = Sed(flambda=flambda, wavelen=bp_dict['g'].wavelen)
        mags_true = bp_dict.magListForSed(ss)

        i_x1 = np.round((x1-param_mins[0])/d_params[0]).astype(int)
        i_c0 = np.round((c0-param_mins[1])/d_params[1]).astype(int)
        i_z = np.round((z-param_mins[2])/d_params[2]).astype(int)
        d_mag = abs_mag-param_mins[3]
        
        tag = i_x1+i_c0*100+i_z*10000
        mag_grid = in_file['%d' % tag].value
        interp_mags = np.zeros(6, dtype=float)
        for i_bp in range(6):
            mm = np.interp(t, t_grid, mag_grid[i_bp])+d_mag
            print('interped %e actual %e -- %e' %
                  (mm, mags_true[i_bp],mm-mags_true[i_bp]))
