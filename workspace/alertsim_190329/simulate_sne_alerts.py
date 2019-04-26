import h5py
import pickle
import numpy as np
import os

import time

from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed
from lsst.sims.photUtils import SignalToNoise as SNR
from lsst.sims.utils import htmModule as htm
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import ModifiedJulianDate
from lsst.sims.catUtils.baseCatalogModels.LocalGalaxyModels import LocalGalaxyTileObj
from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels

from lsst.sims.photUtils import CosmologyObject
from lsst.sims.catUtils.supernovae import SNUniverse
from lsst.sims.catUtils.dust import EBVbase

from alert_focal_plane import apply_focal_plane

import multiprocessing

import numbers
import argparse

_ct_sne = 0


class LocalSNeTileObj(LocalGalaxyTileObj):

    # place holder SNe parameters
    columns = [('t0', '0.0', float),
               ('c0', '0.0', float),
               ('x1', '0.0', float),
               ('abs_mag', '0.0', float),
               ('A_u', '0.0', float),
               ('A_g', '0.0', float),
               ('A_r', '0.0', float),
               ('A_i', '0.0', float),
               ('A_z', '0.0', float),
               ('A_y', '0.0', float)]


def process_sne_chunk(chunk, filter_obs, mjd_obs, m5_obs,
                      coadd_m5, obs_md_list, proper_chip,
                      invisible_set, out_data):

    sne_interp_file = 'data/sne_interp_models.h5'

    global _ct_sne

    n_t = len(filter_obs)
    n_obj = len(chunk)

    with h5py.File('data/ebv_grid.h5', 'r') as in_file:
       ebv_grid = in_file['ebv_grid'].value
       extinction_grid = in_file['extinction_grid'].value

    dust_model = EBVbase()
    ebv = dust_model.calculateEbv(equatorialCoordinates=np.array([np.radians(chunk['ra']),
                                                                  np.radians(chunk['dec'])]),
                                  interp=True)
    for i_bp, bp in enumerate('ugrizy'):
        vv = np.interp(ebv, ebv_grid, extinction_grid[i_bp])
        chunk['A_%s' % bp] = vv

    #print('processing %d' % len(chunk))
    ct_first = 0
    ct_at_all = 0
    ct_tot = 0
    coadd_visits = {}
    coadd_visits['u'] = 6
    coadd_visits['g'] = 8
    coadd_visits['r'] = 18
    coadd_visits['i'] = 18
    coadd_visits['z'] = 16
    coadd_visits['y'] = 16

    # from the overview paper
    # table 2; take m5 row and add Delta m5 row
    # to get down to airmass 1.2
    m5_single = {}
    m5_single['u'] = 23.57
    m5_single['g'] = 24.65
    m5_single['r'] = 24.21
    m5_single['i'] = 23.79
    m5_single['z'] = 23.21
    m5_single['y'] = 22.31

    gamma_coadd = {}
    for bp in 'ugrizy':
        gamma_coadd[bp] = None

    gamma_single = {}
    for bp in 'ugrizy':
       gamma_single[bp] = [None]*n_t

    n_t_per_filter = {}
    t_obs_arr = {}
    i_obs_per_filter = {}
    for i_bp, bp in enumerate('ugrizy'):
        valid = np.where(filter_obs==i_bp)
        n_t_per_filter[bp] = len(valid[0])
        i_obs_per_filter[bp] = valid[0]
        if n_t_per_filter[bp] == 0:
            continue

        t_obs_arr[bp] = mjd_obs[valid]

    # first just need to interpolate some stuff
    ct_invis = 0
    d_mag = np.zeros((n_obj, n_t), dtype=float)
    photometry_mask = np.zeros((n_obj, n_t), dtype=bool)
    photometry_mask_1d = np.zeros(n_obj, dtype=bool)
    with h5py.File(sne_interp_file, 'r') as in_file:
        param_mins = in_file['param_mins'].value
        d_params = in_file['d_params'].value
        t_grid = in_file['t_grid'].value
        abs_mag_0 = param_mins[3]
        i_x_arr = np.round((chunk['x1']-param_mins[0])/d_params[0]).astype(int)
        i_c_arr = np.round((chunk['c0']-param_mins[1])/d_params[1]).astype(int)
        i_z_arr = np.round((chunk['redshift']-param_mins[2])/d_params[2]).astype(int)
        model_tag = i_x_arr+100*i_c_arr+10000*i_z_arr

        unq_tag = np.unique(model_tag)
        for i_tag in unq_tag:
            valid_obj = np.where(model_tag == i_tag)
            if i_tag in invisible_set:
                ct_invis += len(valid_obj[0])
                continue
            d_abs_mag = chunk['abs_mag'][valid_obj]-abs_mag_0

            mag_grid = in_file['%d' % i_tag].value
            for i_bp, bp in enumerate('ugrizy'):
                if n_t_per_filter[bp] == 0:
                    continue
                valid_obs = i_obs_per_filter[bp]
                assert len(valid_obs) == len(t_obs_arr[bp])

                t_matrix = t_obs_arr[bp]-chunk['t0'][valid_obj,None]
                t_arr = t_matrix.flatten()

                sne_mag = np.interp(t_arr,
                                    t_grid, mag_grid[i_bp]).reshape((len(valid_obj[0]),
                                                                     n_t_per_filter[bp]))

                for ii in range(len(valid_obj[0])):
                    obj_dex = valid_obj[0][ii]
                    sne_mag[ii] += d_abs_mag[ii]+chunk['A_%s' % bp][ii]
                    d_mag[obj_dex, valid_obs] = sne_mag[ii]-m5_single[bp]
                    photometry_mask[obj_dex, valid_obs] = sne_mag[ii]<m5_single[bp]

    for i_obj in range(n_obj):
        if photometry_mask[i_obj,:].any():
            photometry_mask_1d[i_obj] = True


    n_unq = len(unq_tag)
    ct_detected = 0

    t_before_chip = time.time()
    chip_mask = apply_focal_plane(chunk['ra'], chunk['dec'],
                                  photometry_mask_1d, obs_md_list,
                                  filter_obs, proper_chip)
    duration = (time.time()-t_before_chip)/3600.0

    unq_arr = -1*np.ones(n_obj, dtype=int)
    mjd_arr = -1.0*np.ones(n_obj, dtype=float)
    redshift_arr = -1.0*np.ones(n_obj, dtype=float)
    snr_arr = -1.0*np.ones(n_obj, dtype=float)

    for i_obj in range(n_obj):
        if photometry_mask_1d[i_obj]:
            detected = photometry_mask[i_obj,:] & chip_mask[i_obj,:]
            if detected.any():
                unq_arr[i_obj] = chunk['galtileid'][i_obj]
                first_dex = np.where(detected)[0].min()
                snr_arr[i_obj] = 5.0*10**(-0.4*(d_mag[i_obj,first_dex]))
                mjd_arr[i_obj] = mjd_obs[first_dex]
                redshift_arr[i_obj] = chunk['redshift'][i_obj]

    valid = np.where(unq_arr>=0)
    pid = os.getpid()
    out_data[pid] = (unq_arr[valid], mjd_arr[valid],
                     snr_arr[valid], redshift_arr[valid])


if __name__ == "__main__":

    t_start = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--out_name', type=str, default=None)
    parser.add_argument('--circular_fov', default=False,
                        action='store_true')
    parser.add_argument('--fast_t0', default=False,
                        action='store_true',
                        help='use an int stored on the database to '
                        'quickly query for only a random 10% of galaxies')
    parser.add_argument('--q_chunk_size', type=int, default=50000,
                        help='number of galaxies to query from '
                             'database at once (default 5*10**4)')
    parser.add_argument('--p_chunk_size', type=int, default=10000,
                        help='number of galaxies to process at once '
                             '(default 10**4)')
    parser.add_argument('--n_threads', type=int, default=30)
    parser.add_argument('--htmid', type=int, nargs='+', default=None)
    parser.add_argument('--htmid_file', type=str, default=None)

    args = parser.parse_args()
    proper_chip = not args.circular_fov
    assert args.out_name is not None
    assert ((args.htmid is not None and args.htmid_file is None) or
            (args.htmid is None and args.htmid_file is not None))

    if args.htmid is not None:
        if isinstance(args.htmid, numbers.Number):
            htmid_list = [args.htmid]
        else:
            htmid_list = args.htmid
    else:
        htmid_list = []
        with open(args.htmid_file, 'r') as in_file:
            for line in in_file:
                if line.startswith('#'):
                    continue
                params = line.strip().split()
                htmid_list.append(int(params[0]))


    if args.fast_t0:
        fast_t0_rng = np.random.RandomState(max(htmid_list))

    invisible_file = 'data/invisible_sn_tags.txt'
    invisible_tags = set()
    with open(invisible_file, 'r') as in_file:
        for line in in_file:
            invisible_tags.add(int(line.strip()))

    coadd_m5_name = 'data/coadd_m5.txt'
    coadd_m5 = {}
    with open(coadd_m5_name, 'r') as in_file:
        for line in in_file:
            if line.startswith('#'):
                continue
            p = line.strip().split()
            coadd_m5[p[0]] = float(p[1])

    htmid_map_name = 'data/htmid_to_obs_map.pickle'
    assert os.path.isfile(htmid_map_name)
    with open(htmid_map_name, 'rb') as in_file:
        htmid_to_obs = pickle.load(in_file)

    print('%d htmid' % len(htmid_to_obs))


    try:
        gal_db = LocalSNeTileObj(database='LSST',
                                 host='epyc.astro.washington.edu',
                                 port=1433,
                                 driver='mssql+pymssql')
    except:
        gal_db = LocalSNeTileObj(database='LSST',
                                 host='localhost',
                                 port=51432,
                                 driver='mssql+pymssql')

    obs_param_name = 'data/obs_params.h5'

    rng = np.random.RandomState(htmid_list[0])

    if args.n_threads>1:
        mgr = multiprocessing.Manager()
        out_data = mgr.dict()
    else:
        out_data = {}

    n_tot = 0
    n_processed = 0

    with h5py.File(obs_param_name, 'r') as obs_params:
        for htmid_query in htmid_list:
            if htmid_query not in htmid_to_obs:
                continue
            print(htmid_query)
            query_level = htm.levelFromHtmid(htmid_query)
            trixel_query = htm.trixelFromHtmid(htmid_query)
            ra_query, dec_query = trixel_query.get_center()
            radius_query = trixel_query.get_radius()
            print(ra_query, dec_query, radius_query)

            obs_query = ObservationMetaData(pointingRA=ra_query,
                                            pointingDec=dec_query,
                                            boundType='circle',
                                            boundLength=radius_query)

            col_names = ['id', 'redshift',
                         'ra', 'dec',
                         'u_ab', 'g_ab', 'r_ab',
                         'i_ab', 'z_ab', 'y_ab',
                         't0', 'c0', 'x1', 'abs_mag',
                         'A_u', 'A_g', 'A_r', 'A_i',
                         'A_z', 'A_y']


            assert np.diff(obs_params['obsHistID']).min()>0

            obsid_query = np.array(htmid_to_obs[htmid_query])
            obs_dex = np.searchsorted(obs_params['obsHistID'].value, obsid_query)
            np.testing.assert_array_equal(obs_params['obsHistID'].value[obs_dex],
                                          obsid_query)

            ra_obs = obs_params['ra'].value[obs_dex]
            dec_obs = obs_params['dec'].value[obs_dex]
            mjd_obs = obs_params['mjd'].value[obs_dex]
            rotsky_obs = obs_params['rotSkyPos'].value[obs_dex]
            filter_obs = obs_params['filter'].value[obs_dex]
            m5_obs = obs_params['m5'].value[obs_dex]

            mjd_obj_list = ModifiedJulianDate.get_list(TAI=mjd_obs)
            obs_md_list = []
            for ii in range(len(ra_obs)):
                obs = ObservationMetaData(pointingRA=ra_obs[ii],
                                          pointingDec=dec_obs[ii],
                                          mjd=mjd_obj_list[ii],
                                          rotSkyPos=rotsky_obs[ii],
                                          bandpassName='ugrizy'[filter_obs[ii]])
                obs_md_list.append(obs)

            print('%d time steps' % len(filter_obs))

            constraint = 'redshift<=1.2 '
            if args.fast_t0:
                fast_t0_flag = fast_t0_rng.randint(0,10)
                constraint += 'AND sne_alert_flag=%d' % fast_t0_flag

            sn_frequency = 1.0/(100.0*365.0)
            midSurveyTime = 59580.0+5.0*365.25

            data_iter = gal_db.query_columns(col_names, obs_metadata=obs_query,
                                             chunk_size=args.q_chunk_size,
                                             constraint=constraint)
            p_list = []
            i_chunk = 0
            to_concatenate = []

            n_sne = 0
            n_gal = 0
            tot_unq = 0
            tot_det = 0

            t_lsst_0 = 59580.0-34.0
            t_lsst_1 = t_lsst_0+3652.5+100.0

            for chunk in data_iter:
                htmid_found = htm.findHtmid(chunk['ra'],
                                            chunk['dec'],
                                            query_level)

                valid = np.where(htmid_found==htmid_query)
                if len(valid[0]) == 0:
                    continue

                chunk = chunk[valid]

                n_gal += len(chunk)

                if args.fast_t0:
                    t0_arr = rng.uniform(59580.0, 59580.0+3652.5,
                                         size=len(chunk))
                else:
                    t0_arr = rng.uniform(midSurveyTime-0.5/sn_frequency,
                                         midSurveyTime+0.5/sn_frequency,
                                         size=len(chunk))

                is_sne = np.where(np.logical_and(t0_arr>=t_lsst_0,
                                                 t0_arr<=t_lsst_1))

                chunk = chunk[is_sne]
                t0_arr = t0_arr[is_sne]

                n_sne += len(is_sne[0])
                dt_matrix = mjd_obs-t0_arr[:,None]

                valid = np.where((dt_matrix>-34.0).any(axis=1) & (dt_matrix<100.0).any(axis=1))

                chunk['t0'] = t0_arr

                chunk = chunk[valid]

                c0_arr = np.clip(rng.normal(0.0, 0.1, size=len(chunk)), -0.3, 0.3)
                x1_arr = np.clip(rng.normal(0.0, 1.0, size=len(chunk)), -3.0, 3.0)
                abs_mag_arr = rng.normal(-19.3, 0.3, size=len(chunk))

                chunk['c0'] = c0_arr
                chunk['x1'] = x1_arr
                chunk['abs_mag'] = abs_mag_arr

                n_tot += len(chunk)
                if args.n_threads == 1:
                    process_sne_chunk(chunk, filter_obs,mjd_obs,
                                      m5_obs, coadd_m5, obs_md_list,
                                      proper_chip, invisible_tags, out_data)
                    continue

                # multiprocessing code
                if len(chunk)<args.p_chunk_size:
                    to_concatenate.append(chunk)
                    tot_sub = 0
                    for sub_chunk in to_concatenate:
                        tot_sub += len(sub_chunk)

                    if n_processed+tot_sub != n_tot:
                        raise RuntimeError('n_proc+tot %d n_tot %d'
                                           % (n_processed+tot_sub, n_tot))
                    if tot_sub<args.p_chunk_size:
                        continue
                    else:
                        chunk = np.concatenate(to_concatenate)
                        assert len(chunk)==tot_sub
                        to_concatenate = []

                for i_min in range(0, len(chunk)+1, args.p_chunk_size):
                    sub_chunk = chunk[i_min:i_min+args.p_chunk_size]
                    if len(sub_chunk)<args.p_chunk_size:
                        to_concatenate.append(sub_chunk)
                        continue

                    n_processed += len(sub_chunk)
                    assert len(sub_chunk)>=args.p_chunk_size
                    assert args.n_threads>1
                    p = multiprocessing.Process(target=process_sne_chunk,
                                                args=(sub_chunk, filter_obs, mjd_obs,
                                                      m5_obs, coadd_m5, obs_md_list,
                                                      proper_chip, invisible_tags,
                                                      out_data))
                    p.start()
                    p_list.append(p)
                    while len(p_list)>=args.n_threads:
                        exit_code_list = []
                        for p in p_list:
                            exit_code_list.append(p.exitcode)
                        for i_p in range(len(exit_code_list)-1, -1, -1):
                            if exit_code_list[i_p] is not None:
                                p_list.pop(i_p)

                tot_sub = 0
                for sub_chunk in to_concatenate:
                    tot_sub += len(sub_chunk)
                if n_processed+tot_sub!=n_tot:
                    raise RuntimeError("sums failed after processing %d %d -- %d"
                    % (n_processed+tot_sub,n_tot,tot_sub))

            if len(to_concatenate)>0:
                chunk = np.concatenate(to_concatenate)
                for i_min in range(0,len(chunk),args.p_chunk_size):
                    sub_chunk = chunk[i_min:i_min+args.p_chunk_size]
                    n_processed += len(sub_chunk)

                    assert args.n_threads>1
                    p = multiprocessing.Process(target=process_sne_chunk,
                                                args=(sub_chunk, filter_obs, mjd_obs,
                                                      m5_obs, coadd_m5, obs_md_list,
                                                      proper_chip, invisible_tags,
                                                      out_data))

                    p.start()
                    p_list.append(p)
                    while len(p_list)>=args.n_threads:
                        exit_code_list = []
                        for p in p_list:
                            exit_code_list.append(p.exitcode)
                        for i_p in range(len(exit_code_list)-1, -1, -1):
                            if exit_code_list[i_p] is not None:
                                p_list.pop(i_p)


            for p in p_list:
                p.join()

    if args.n_threads>1:
        out_data_final = {}
        for name in out_data.keys():
            out_data_final[name] = out_data[name]
    else:
        out_data_final = out_data

    print('n_lc %d' % len(out_data_final))
    print('n_sne %e' % n_sne)
    print('n_gal %e' % n_gal)
    unq_arr = []
    mjd_arr = []
    snr_arr = []
    redshift_arr = []
    for name in out_data_final:
        unq_arr.append(out_data_final[name][0])
        mjd_arr.append(out_data_final[name][1])
        snr_arr.append(out_data_final[name][2])
        redshift_arr.append(out_data_final[name][3])
    unq_arr = np.concatenate(unq_arr)
    mjd_arr = np.concatenate(mjd_arr)
    snr_arr = np.concatenate(snr_arr)
    redshift_arr = np.concatenate(redshift_arr)

    with h5py.File(out_name, 'w') as out_file:
        out_file.create_dataset('unq', data=unq_arr)
        out_file.create_dataset('mjd', data=mjd_arr)
        out_file.create_dataset('snr', data=snr_arr)
        out_file.create_dataset('redshift', data=redshift_arr)

    print('that took %e hrs' % ((time.time()-t_start)/3600.0))
    print('shld %d processed %d' % (n_tot, n_processed))
