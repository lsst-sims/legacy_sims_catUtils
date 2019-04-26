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

from lsst.sims.catUtils.dust import EBVbase

from alert_focal_plane import apply_focal_plane

import multiprocessing

import argparse
import numbers

def process_agn_chunk(chunk, filter_obs, mjd_obs, m5_obs,
                      coadd_m5, m5_single,
                      obs_md_list, proper_chip, out_data):
    t_start_chunk = time.time()
    #print('processing %d' % len(chunk))
    ct_first = 0
    ct_at_all = 0
    ct_tot = 0

    n_t = len(filter_obs)
    n_obj = len(chunk)

    agn_model = ExtraGalacticVariabilityModels()
    dust_model = EBVbase()

    with h5py.File('data/ebv_grid.h5', 'r') as in_file:
       ebv_grid = in_file['ebv_grid'].value
       extinction_grid = in_file['extinction_grid'].value

    coadd_visits = {}
    coadd_visits['u'] = 6
    coadd_visits['g'] = 8
    coadd_visits['r'] = 18
    coadd_visits['i'] = 18
    coadd_visits['z'] = 16
    coadd_visits['y'] = 16

    gamma_coadd = {}
    for bp in 'ugrizy':
        gamma_coadd[bp] = None

    gamma_single = {}
    for bp in 'ugrizy':
       gamma_single[bp] = [None]*n_t

    params = {}
    params['agn_sfu'] = chunk['agn_sfu']
    params['agn_sfg'] = chunk['agn_sfg']
    params['agn_sfr'] = chunk['agn_sfr']
    params['agn_sfi'] = chunk['agn_sfi']
    params['agn_sfz'] = chunk['agn_sfz']
    params['agn_sfy'] = chunk['agn_sfy']
    params['agn_tau'] = chunk['agn_tau']
    params['seed'] = chunk['id']+1

    ebv = dust_model.calculateEbv(equatorialCoordinates=np.array([np.radians(chunk['ra']),
                                                                  np.radians(chunk['dec'])]),
                                  interp=True)

    for i_bp, bp in enumerate('ugrizy'):
        extinction_values = np.interp(ebv, ebv_grid, extinction_grid[i_bp])
        chunk['%s_ab' % bp] += extinction_values
        chunk['AGNLSST%s' % bp] += extinction_values

    dmag = agn_model.applyAgn(np.where(np.array([True]*len(chunk))),
                              params, mjd_obs,
                              redshift=chunk['redshift'])

    dmag_mean = np.mean(dmag, axis=2)
    assert dmag_mean.shape == (6,n_obj)

    dummy_sed = Sed()
    lsst_bp = BandpassDict.loadTotalBandpassesFromFiles()
    flux_gal = np.zeros((6,n_obj), dtype=float)
    flux_agn_q = np.zeros((6,n_obj), dtype=float)
    flux_coadd = np.zeros((6,n_obj), dtype=float)
    mag_coadd = np.zeros((6,n_obj), dtype=float)
    snr_coadd = np.zeros((6,n_obj), dtype=float)
    snr_single = {}
    snr_single_mag_grid = np.arange(14.0, 30.0, 0.05)

    phot_params_single = PhotometricParameters(nexp=1,
                                               exptime=30.0)

    t_start_snr = time.time()

    for i_bp, bp in enumerate('ugrizy'):
        phot_params_coadd = PhotometricParameters(nexp=1,
                                                  exptime=30.0*coadd_visits[bp])

        flux_gal[i_bp] = dummy_sed.fluxFromMag(chunk['%s_ab' % bp])
        flux_agn_q[i_bp] = dummy_sed.fluxFromMag(chunk['AGNLSST%s' % bp] +
                                                 dmag_mean[i_bp,:])
        flux_coadd[i_bp] = flux_gal[i_bp]+flux_agn_q[i_bp]
        mag_coadd[i_bp] = dummy_sed.magFromFlux(flux_coadd[i_bp])

        (snr_coadd[i_bp],
         gamma) = SNR.calcSNR_m5(mag_coadd[i_bp],
                                 lsst_bp[bp],
                                 coadd_m5[bp],
                                 phot_params_coadd)


        (snr_single[bp],
         gamma) = SNR.calcSNR_m5(snr_single_mag_grid,
                                 lsst_bp[bp],
                                 m5_single[bp],
                                 phot_params_single)

    #print('got all snr in %e' % (time.time()-t_start_snr))


    t_start_obj = time.time()
    photometry_mask = np.zeros((n_obj, n_t), dtype=bool)
    photometry_mask_1d = np.zeros(n_obj, dtype=bool)
    snr_arr = np.zeros((n_obj, n_t), dtype=float)
    for i_bp, bp in enumerate('ugrizy'):
        valid_obs = np.where(filter_obs==i_bp)
        n_bp = len(valid_obs[0])
        if n_bp == 0:
            continue
        mag0_arr = chunk['AGNLSST%s' % bp]
        dmag_bp = dmag[i_bp][:,valid_obs[0]]
        assert dmag_bp.shape == (n_obj, n_bp)
        agn_flux_tot = dummy_sed.fluxFromMag(mag0_arr[:,None]+dmag_bp)
        q_flux = flux_agn_q[i_bp]
        agn_dflux = np.abs(agn_flux_tot-q_flux[:,None])
        flux_tot = flux_gal[i_bp][:, None] + agn_flux_tot
        assert flux_tot.shape == (n_obj, n_bp)
        mag_tot = dummy_sed.magFromFlux(flux_tot)
        snr_single_val = np.interp(mag_tot,
                                   snr_single_mag_grid,
                                   snr_single[bp])

        noise_coadd = flux_coadd[i_bp]/snr_coadd[i_bp]
        noise_single = flux_tot/snr_single_val
        noise = np.sqrt(noise_coadd[:,None]**2 + noise_single**2)
        dflux_thresh = 5.0*noise
        detected = (agn_dflux>=dflux_thresh)
        assert detected.shape == (n_obj, n_bp)
        snr_arr[:,valid_obs[0]] = agn_dflux/noise
        for i_obj in range(n_obj):
            if detected[i_obj].any():
                photometry_mask_1d[i_obj] = True
                photometry_mask[i_obj, valid_obs[0]] = detected[i_obj]

    t_before_chip = time.time()
    chip_mask = apply_focal_plane(chunk['ra'], chunk['dec'],
                                  photometry_mask_1d, obs_md_list,
                                  filter_obs, proper_chip)
    duration = (time.time()-t_before_chip)/3600.0

    unq_out = -1*np.ones(n_obj, dtype=int)
    mjd_out = -1.0*np.ones(n_obj, dtype=float)
    snr_out = -1.0*np.ones(n_obj, dtype=float)
    for i_obj in range(n_obj):
        if photometry_mask_1d[i_obj]:
            detected = photometry_mask[i_obj,:] & chip_mask[i_obj,:]

            if detected.any():
                unq_out[i_obj] = chunk['galtileid'][i_obj]
                first_dex = np.where(detected)[0].min()
                mjd_out[i_obj] = mjd_obs[first_dex]
                snr_out[i_obj] = snr_arr[i_obj, first_dex]

    valid = np.where(unq_out>=0)
    pid = os.getpid()
    out_data[pid] = (unq_out[valid], mjd_out[valid],
                     snr_out[valid])


if __name__ == "__main__":

    t_start = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--out_name', type=str, default=None)
    parser.add_argument('--circular_fov', default=False,
                        action='store_true')
    parser.add_argument('--m5_single', type=str, default='data/single_m5.txt',
                        help='File containing single visit m5 values')
    parser.add_argument('--coadd_m5', type=str, default='data/coadd_m5.txt')
    parser.add_argument('--q_chunk_size', type=int, default=10000,
                        help='number of galaxies to query from '
                             'database at once (default 10**4)')
    parser.add_argument('--p_chunk_size', type=int, default=1000,
                        help='number of galaxies to process at once '
                             '(default 10**3)')
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

    m5_single = {}
    with open(args.m5_single, 'r') as in_file:
        for line in in_file:
            if line.startswith('#'):
                continue
            p = line.strip().split()
            m5_single[p[0]] = float(p[1])

    coadd_m5 = {}
    with open(args.coadd_m5, 'r') as in_file:
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
        gal_db = LocalGalaxyTileObj(database='LSST',
                                    host='epyc.astro.washington.edu',
                                    port=1433,
                                    driver='mssql+pymssql')
    except:
        gal_db = LocalGalaxyTileObj(database='LSST',
                                    host='localhost',
                                    port=51432,
                                    driver='mssql+pymssql')


    mgr = multiprocessing.Manager()
    out_data = mgr.dict()
    obs_param_name = 'data/obs_params.h5'
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

            col_names = ['agn_sfu', 'agn_sfg', 'agn_sfr',
                         'agn_sfi', 'agn_sfz', 'agn_sfy',
                         'agn_tau', 't0_agn', 'id', 'redshift',
                         'ra', 'dec',
                         'u_ab', 'g_ab', 'r_ab',
                         'i_ab', 'z_ab', 'y_ab',
                         'AGNLSSTu', 'AGNLSSTg', 'AGNLSSTr',
                         'AGNLSSTi', 'AGNLSSTz', 'AGNLSSTy']

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

            constraint = 'isagn=1 '

            data_iter = gal_db.query_columns(col_names, obs_metadata=obs_query,
                                             chunk_size=args.q_chunk_size,
                                             constraint=constraint)

            p_list = []
            i_chunk = 0
            to_concatenate = []
            n_tot = 0
            n_processed = 0
            for chunk in data_iter:
                htmid_found = htm.findHtmid(chunk['ra'],
                                            chunk['dec'],
                                            query_level)

                valid = np.where(htmid_found==htmid_query)
                if len(valid[0]) == 0:
                    continue

                chunk = chunk[valid]
                n_tot += len(chunk)

                #process_agn_chunk(chunk, filter_obs, mjd_obs, m5_obs, coadd_m5,
                #                  out_data)

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
                    p = multiprocessing.Process(target=process_agn_chunk,
                                                args=(sub_chunk, filter_obs, mjd_obs,
                                                      m5_obs, coadd_m5, 
                                                      m5_single, obs_md_list,
                                                      proper_chip, out_data))
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
                    p = multiprocessing.Process(target=process_agn_chunk,
                                                args=(sub_chunk,
                                                      filter_obs, mjd_obs,
                                                      m5_obs, coadd_m5,
                                                      m5_single, obs_md_list,
                                                      proper_chip, out_data))
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

    unq_arr = []
    mjd_arr = []
    snr_arr = []
    for name in out_data.keys():
        unq_arr.append(out_data[name][0])
        mjd_arr.append(out_data[name][1])
        snr_arr.append(out_data[name][2])
    unq_arr = np.concatenate(unq_arr)
    mjd_arr = np.concatenate(mjd_arr)
    snr_arr = np.concatenate(snr_arr)

    print('n_lc %d' % len(unq_arr))

    with h5py.File(args.out_name, 'w') as out_file:
        out_file.create_dataset('unq', data=unq_arr)
        out_file.create_dataset('mjd', data=mjd_arr)
        out_file.create_dataset('snr', data=snr_arr)

    print('that took %e hrs' % ((time.time()-t_start)/3600.0))
    print('shld %d processed %d' % (n_tot, n_processed))
    obs_params.close()
