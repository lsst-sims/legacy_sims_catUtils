#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

import h5py

import sncosmo
import numpy as np
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed
from lsst.sims.photUtils import CosmologyObject
import time

rng = np.random.RandomState(654)
cosmo = CosmologyObject()

mag_min = -20.2
c0_min = -0.3
x1_min = -3.0
z_min = 0.01

mag_range = 1.8
c0_range = 0.6
x1_range = 6.0
z_range = 1.2

d_c0 = 0.025
d_x1 = 1.0
d_z=0.01

bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
wav_grid = bp_dict['g'].wavelen

t0_grid = np.arange(-100.0, -21.0, 10.0)
t1_grid = np.arange(-20.0, 41.0, 2.0)
t2_grid = np.arange(45.0, 106.0, 10.0)

t_grid = np.sort(np.concatenate([t0_grid,t1_grid,t2_grid]))
assert np.diff(t_grid).min()>0.9

z_grid = np.arange(z_min, z_min+z_range+0.1*d_z, d_z)
c0_grid = np.arange(c0_min, c0_min+c0_range+0.1*d_c0, d_c0)
x1_grid = np.arange(x1_min, x1_min+x1_range+0.1*d_x1, d_x1)

ct_models = len(z_grid)*len(c0_grid)*len(x1_grid)
print('generating %e models' % ct_models)

abs_mag = -19.0

t_start = time.time()
models_generated = 0
tags_used = set()
max_tag = 2**32
with h5py.File('data/sne_interp_models.h5', 'w') as out_file:
    out_file.create_dataset('t_grid', data=t_grid)
    out_file.create_dataset('param_mins',
                            data=np.array([x1_min, c0_min, z_min, abs_mag]))
    out_file.create_dataset('d_params',
                            data=np.array([d_x1, d_c0, d_z]))
    out_file.create_dataset('param_names',
            data=np.array(['x1','c0', 'z', 'abs_mag']).astype(np.string_))

    for z in z_grid:
        i_z = np.round((z-z_min)/d_z).astype(int)
        dm = cosmo.distanceModulus(z)
        assert i_z>=0
        for c0 in c0_grid:
            i_c0 = np.round((c0-c0_min)/d_c0).astype(int)
            assert i_c0>=0
            for x1 in x1_grid:
                i_x1 = np.round((x1-x1_min)/d_x1).astype(int)
                assert i_x1>=0
                tag = i_x1+i_c0*100+i_z*10000
                assert tag not in tags_used
                assert tag<max_tag
                tags_used.add(tag)
                mag_grid = np.zeros((6, len(t_grid)), dtype=float)

                sn = sncosmo.Model(source='salt2-extended')
                sn.set(x1=x1,c=c0,z=z)
                sn.source.set_peakmag(abs_mag+dm, band='bessellb', magsys='ab')

                for i_t, t in enumerate(t_grid):
                    flambda = 10.0*sn.flux(time=t, wave=10.0*wav_grid)
                    ss = Sed(flambda=flambda, wavelen=wav_grid)

                    mags = bp_dict.magListForSed(ss)
                    for i_bp in range(6):
                        mag_grid[i_bp][i_t] = mags[i_bp]

                mag_grid = np.where(np.isnan(mag_grid),
                                    99.0, mag_grid)

                out_file.create_dataset('%d' % tag, data=mag_grid)

                models_generated += 1
                if models_generated % 100 == 0:
                    duration = (time.time()-t_start)/3600.0
                    print('%d of %d in %e hrs; proj %e hrs' %
                            (models_generated, ct_models, duration,
                             ct_models*duration/models_generated))

exit()

duration = time.time()-t_start
print('one sn took %e' % duration)
print('all will take %e' % (duration*ct_models))

for bp in 'ugrizy':
    mag_grid[bp] = np.where(np.isnan(mag_grid[bp]),99.0,mag_grid[bp])

plus_c0 = c0+0.5*d_c0 #(rng.random_sample()-0.5)*d_c0
plus_x1 = x1+0.5*d_x1 #(rng.random_sample()-0.5)*d_x1
plus_absmag = absmag+0.5*d_mag #(rng.random_sample()-0.5)*d_mag
plus_z = z+0.5*d_z #(rng.random_sample()-0.5)*d_z
plus_dm = cosmo.distanceModulus(plus_z)

plus_sn = sncosmo.Model(source='salt2-extended')
plus_sn.set(x1=plus_x1,c=plus_c0,z=plus_z)
plus_sn.source.set_peakmag(plus_absmag+plus_dm, band='bessellb', magsys='ab')

minus_c0 = c0-0.5*d_c0 #(rng.random_sample()-0.5)*d_c0
minus_x1 = x1-0.5*d_x1 #(rng.random_sample()-0.5)*d_x1
minus_absmag = absmag-0.5*d_mag #(rng.random_sample()-0.5)*d_mag
minus_z = z-0.5*d_z #(rng.random_sample()-0.5)*d_z
minus_dm = cosmo.distanceModulus(minus_z)

minus_sn = sncosmo.Model(source='salt2-extended')
minus_sn.set(x1=minus_x1,c=minus_c0,z=minus_z)
minus_sn.source.set_peakmag(minus_absmag+minus_dm, band='bessellb', magsys='ab')


plus_mag_grid = {}
minus_mag_grid = {}
for bp in 'ugrizy':
    plus_mag_grid[bp] = np.zeros(len(t_grid), dtype=float)
    minus_mag_grid[bp] = np.zeros(len(t_grid), dtype=float)

for i_t, t in enumerate(t_grid):
    flambda = 10.0*plus_sn.flux(time=t, wave=10.0*wav_grid)
    ss = Sed(flambda=flambda, wavelen=wav_grid)

    mags = bp_dict.magListForSed(ss)
    for i_bp, bp in enumerate('ugrizy'):
        plus_mag_grid[bp][i_t] = mags[i_bp]


    flambda = 10.0*minus_sn.flux(time=t, wave=10.0*wav_grid)
    ss = Sed(flambda=flambda, wavelen=wav_grid)

    mags = bp_dict.magListForSed(ss)
    for i_bp, bp in enumerate('ugrizy'):
        minus_mag_grid[bp][i_t] = mags[i_bp]


for bp in 'ugrizy':
    plus_mag_grid[bp] = np.where(np.isnan(plus_mag_grid[bp]),
                                99.0,plus_mag_grid[bp])

    minus_mag_grid[bp] = np.where(np.isnan(minus_mag_grid[bp]),
                                99.0,minus_mag_grid[bp])


plt.figure(figsize=(20,20))
for i_bp, bp in enumerate('ugrizy'):
    plt.subplot(3,2,i_bp+1)
    plt.plot(t_grid,mag_grid[bp], linestyle='-')
    plt.plot(t_grid,plus_mag_grid[bp], linestyle='--')
    plt.plot(t_grid,minus_mag_grid[bp], linestyle='--')
    plt.title(bp,fontsize=20)

plt.tight_layout()
plt.savefig('interped_mags.png')

print('c0 ',c0,plus_c0,minus_c0)
print('x1 ',x1,plus_x1,minus_x1)
print('z ',z,plus_z,minus_z)
print('absmag ',absmag,plus_absmag,minus_absmag)
print('ct models %e' % ct_models)
print('floats %e' % (ct_models*len(t_grid)))

