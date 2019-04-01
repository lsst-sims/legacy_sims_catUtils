import h5py
import pickle
import numpy as np
import hashlib
import time

import argparse

from lsst.sims.utils import htmModule as htm
from lsst.sims.utils import angularSeparation

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fov_radius', type=float, default=2.1,
                        help='Field of view radius in degrees '
                        '(default=2.1)')
    parser.add_argument('--obs_params', type=str, default=None,
                        help='hdf5 file containing pointing parameters')
    parser.add_argument('--out_file', type=str, default=None,
                        help='file to write pickled dict mapping '
                        'htmid to obsHistID')


    args = parser.parse_args()
    if args.obs_params is None:
        raise RuntimeError("must specify obs_params")
    if args.out_file is None:
        raise RuntimeError("must specify out_file")

    trixel_level = 6
    trixel_dict = htm.getAllTrixels(trixel_level)
    trixel_arr = []
    trixel_ra = []
    trixel_dec = []
    trixel_radius = []
    for htmid in trixel_dict:
        if htm.levelFromHtmid(htmid) == trixel_level:
            tx = trixel_dict[htmid]
            trixel_arr.append(tx)
            ra, dec = tx.get_center()
            trixel_ra.append(ra)
            trixel_dec.append(dec)
            trixel_radius.append(tx.get_radius())
    del trixel_dict

    trixel_arr = np.array(trixel_arr)
    trixel_ra = np.array(trixel_ra)
    trixel_dec = np.array(trixel_dec)
    trixel_radius = np.array(trixel_radius)

    htmid_to_obs = {}

    md5_hasher = hashlib.md5()
    with open(args.obs_params, 'rb') as in_file:
        for chunk in iter(lambda: in_file.read(4096), b""):
            md5_hasher.update(chunk)
    md5_sum = md5_hasher.hexdigest()

    htmid_to_obs['md5_sum'] = md5_sum

    t_start = time.time()
    with h5py.File(args.obs_params, 'r') as obs_params:
        obs_ra = obs_params['ra'].value
        obs_dec = obs_params['dec'].value
        obs_id = obs_params['obsHistID'].value
        if obs_ra.max()<10.0:
            raise RuntimeError("I think RA is in radians")
        if obs_dec.min()>-10.0:
            raise RuntimeError("I think Dec is in radians")

        for i_obs, (ra, dec, obsid) in enumerate(zip(obs_ra, obs_dec, obs_id)):
            if i_obs>0 and i_obs%1000 == 0:
                duration = (time.time()-t_start)/3600.0
                print('ran %d in %e hrs' % (i_obs, duration))

            dd = angularSeparation(ra, dec, trixel_ra, trixel_dec)
            valid = np.where(dd<=args.fov_radius+trixel_radius)
            candidate_trixels = trixel_arr[valid]
            hs = htm.halfSpaceFromRaDec(ra, dec, args.fov_radius)
            for tx in candidate_trixels:
                if hs.contains_trixel(tx) != 'outside':
                    htmid = tx.htmid
                    if htmid not in htmid_to_obs:
                        htmid_to_obs[htmid] = []
                    htmid_to_obs[htmid].append(obsid)

    with open(args.out_file, 'wb') as out_file:
        pickle.dump(htmid_to_obs, out_file)
