import numpy as np
import sncosmo
import h5py
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed
from lsst.sims.photUtils import CosmologyObject

import os
import multiprocessing

import argparse

import time

def test_sne(grid_name, x1_vals, c0_vals,
             z_vals, abs_mag_vals, t_vals,
             dict_key, out_dict):

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    cosmo = CosmologyObject()

    n_samples = len(z_vals)

    mag_truth = np.zeros((6,n_samples), dtype=float)
    mag_interp = np.zeros((6,n_samples), dtype=float)

    mag_grid_dict = {}
    with h5py.File(grid_name, 'r') as in_file:
        param_mins = in_file['param_mins'].value
        d_params = in_file['d_params'].value
        t_grid = in_file['t_grid'].value
        for name in in_file.keys():
            if name == 'param_mins':
                continue
            if name == 'd_params':
                continue
            if name == 't_grid':
                continue
            mag_grid_dict[name] = in_file[name].value

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
        mag_grid = mag_grid_dict['%d' % tag]
        interp_mags = np.zeros(6, dtype=float)
        for i_bp in range(6):
            mm = np.interp(t, t_grid, mag_grid[i_bp])+d_mag
            mag_interp[i_bp,ii] = mm


        if ii>0 and ii%100 == 0:
            duration = (time.time()-t_start)/3600.0
            pred = n_samples*duration/ii
            print('%d in %e hrs; predict %e' % (ii,duration,pred))

    out_dict[dict_key] = (mag_interp, mag_truth)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_samples', type=int, default=1000000)
    parser.add_argument('--n_processes', type=int, default=30)
    args = parser.parse_args()


    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
    cosmo = CosmologyObject()

    grid_name = 'sne_interp_models.h5'
    rng = np.random.RandomState(44)

    c0_range = 0.6
    x1_range = 6.0
    z_range = 1.2

    abs_mag_mean = -19.3
    abs_mag_std = 0.3

    d_samples = args.n_samples//args.n_processes

    with h5py.File(grid_name, 'r') as in_file:
        param_mins = in_file['param_mins'].value
        d_params = in_file['d_params'].value
        t_grid = in_file['t_grid'].value

        x1_vals = param_mins[0]+x1_range*rng.random_sample(args.n_samples)
        c0_vals = param_mins[1]+c0_range*rng.random_sample(args.n_samples)
        z_vals = param_mins[2]+z_range*rng.random_sample(args.n_samples)
        abs_mag_vals = -20.2+rng.random_sample(args.n_samples)*1.8
        t_vals = -20.0+rng.random_sample(args.n_samples)*40.0


    out_dict = multiprocessing.Manager().dict()
    t_start = time.time()
    p_list = []
    dict_key_list = []
    dict_key = 0
    for i_min in range(0,args.n_samples,d_samples):
        p = multiprocessing.Process(target=test_sne,
                                    args=(grid_name,
                                          x1_vals[i_min:i_min+d_samples],
                                          c0_vals[i_min:i_min+d_samples],
                                          z_vals[i_min:i_min+d_samples],
                                          abs_mag_vals[i_min:i_min+d_samples],
                                          t_vals[i_min:i_min+d_samples],
                                          dict_key, out_dict))
        dict_key_list.append(dict_key)
        dict_key += 1
        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    mag_truth = np.concatenate([out_dict[kk][1] for kk in dict_key_list],
                               axis=1)
    mag_interp = np.concatenate([out_dict[kk][0] for kk in dict_key_list],
                                axis=1)

    assert mag_truth.shape[1] == len(t_vals)
    assert mag_interp.shape[1]== len(t_vals)

    with h5py.File('sne_interp_test_data.h5', 'w') as out_file:
        out_file.create_dataset('t', data=t_vals)
        out_file.create_dataset('c0', data=c0_vals)
        out_file.create_dataset('x1', data=x1_vals)
        out_file.create_dataset('z', data=z_vals)
        out_file.create_dataset('mag_true', data=mag_truth)
        out_file.create_dataset('mag_interp', data=mag_interp)
