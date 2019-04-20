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
from lsst.sims.catUtils.baseCatalogModels.LocalStarModels import LocalStarCatalogObj

from lsst.sims.catUtils.mixins import create_variability_cache
from lsst.sims.catUtils.mixins import ParametrizedLightCurveMixin
from lsst.sims.catUtils.mixins import StellarVariabilityModels
from lsst.sims.catUtils.mixins import MLTflaringMixin

from alert_focal_plane import apply_focal_plane

import multiprocessing

import numbers
import argparse


def dflux_for_mlt(chunk, filter_obs, mjd_obs, variability_cache, dflux_out):

    valid_obj = np.where(chunk['var_type'] == 2)
    n_mlt = len(valid_obj[0])
    if n_mlt == 0:
        return

    lc_id_map = {}
    for ii in range(4):
        lc_id_map[10+ii] = 'early_inactive_%d' % ii
        lc_id_map[20+ii] = 'early_active_%d' % ii
        lc_id_map[30+ii] = 'mid_inactive_%d' % ii
        lc_id_map[40+ii] = 'mid_active_%d' % ii
        lc_id_map[50+ii] = 'late_active_%d' % ii

    mlt_model = MLTflaringMixin()
    mlt_model.photParams = PhotometricParameters(nexp=1, exptime=30.0)
    mlt_model.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
    mlt_model._actually_calculated_columns = []
    for bp in 'ugrizy':
        mlt_model._actually_calculated_columns.append('lsst_%s' % bp)

    params = {}
    params['lc'] = np.array([lc_id_map[ii] for ii in chunk['lc_id'][valid_obj]])
    params['t0'] = chunk['t0'][valid_obj]
    q_mags = {}
    for bp in 'ugrizy':
        q_mags[bp] = chunk['%smag' % bp][valid_obj]

    for i_bp, bp in enumerate('ugrizy'):
        valid_bp = np.where(filter_obs==i_bp)
        if len(valid_bp[0]) == 0:
            continue

        results = mlt_model.applyMLTflaring(np.array([True]*n_mlt),
                                            params, mjd_obs[valid_bp],
                                            parallax=chunk['parallax'][valid_obj],
                                            ebv=chunk['ebv'][valid_obj],
                                            quiescent_mags=q_mags,
                                            variability_cache=variability_cache,
                                            do_mags=False,
                                            mag_name_tuple=(bp,))

        for i_local, i_global in enumerate(valid_obj[0]):
            dflux_out[i_global][valid_bp] = results[0][i_local]


def dflux_for_kepler(chunk, filter_obs, mjd_obs, v_cache, dflux_out):
    valid_obj = np.where(chunk['var_type']==1)
    n_obj = len(valid_obj[0])
    if n_obj==0:
        return

    params = {}
    params['lc'] = chunk['lc_id'][valid_obj]
    params['t0'] = chunk['t0'][valid_obj]

    plc_model = ParametrizedLightCurveMixin()
    dmag = plc_model.singleBandParametrizedLightCurve(np.array([True]*n_obj),
                                                      params, mjd_obs,
                                                      variability_cache=v_cache)

    dummy_sed = Sed()
    for i_bp, bp in enumerate('ugrizy'):
        valid_obs = np.where(filter_obs==i_bp)
        n_t = len(valid_obs[0])
        if n_t==0:
            continue
        flux0 = dummy_sed.fluxFromMag(chunk['%smag' % bp][valid_obj])
        dmag_obs = dmag[:,valid_obs[0]].transpose()
        assert dmag_obs.shape == (n_t, n_obj)
        flux1 = dummy_sed.fluxFromMag(chunk['%smag' % bp][valid_obj]
                                      + dmag_obs)

        dflux_local = (flux1-flux0).transpose()
        for i_local, i_global in enumerate(valid_obj[0]):
            dflux_out[i_global][valid_obs] = dflux_local[i_local]

def dflux_for_rrly(chunk, filter_obs, mjd_obs, v_cache, dflux_out):
    valid_obj = np.where(chunk['var_type']==3)
    n_obj = len(valid_obj[0])
    if n_obj == 0:
        return

    print('running RRLy')
    params = {}
    params['filename'] = np.array([v_cache['rrly_map'][ii]
                                   for ii in chunk['lc_id'][valid_obj]])

    params['tStartMjd'] = chunk['t0'][valid_obj]

    stellar_model = StellarVariabilityModels()
    stellar_model.initializeVariability(doCache=True)
    dmag = stellar_model.applyRRly(np.where(np.array([True]*n_obj)),
                                   params, mjd_obs,
                                   variability_cache=v_cache).transpose(0,2,1)

    dummy_sed = Sed()
    for i_bp, bp in enumerate('ugrizy'):
        valid_obs = np.where(filter_obs==i_bp)
        if len(valid_obs[0]) == 0:
            continue
        dmag_this_band = dmag[i_bp]
        dmag_subset = dmag_this_band[valid_obs,:][0]
        flux_0 = dummy_sed.fluxFromMag(chunk['%smag' % bp][valid_obj])
        flux_1 = dummy_sed.fluxFromMag(flux_0+dmag_subset).transpose()
        assert flux_1.shape == (n_obj, len(valid_obs[0]))
        for i_local, i_global in enumerate(valid_obj[0]):
            dflux = flux_1[i_local]-flux_0[i_local]
            dflux_out[i_global][valid_obs] = dflux


def process_stellar_chunk(chunk, filter_obs, mjd_obs, m5_obs,
                          coadd_m5, obs_md_list, proper_chip,
                          variability_cache, out_data):

    t_start_chunk = time.time()
    #print('processing %d' % len(chunk))
    ct_first = 0
    ct_at_all = 0
    ct_tot = 0

    n_t = len(filter_obs)
    n_obj = len(chunk)

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

    dflux = np.zeros((n_obj,n_t), dtype=float)
    dflux_for_mlt(chunk, filter_obs, mjd_obs, variability_cache, dflux)
    dflux_for_kepler(chunk, filter_obs, mjd_obs, variability_cache, dflux)
    dflux_for_rrly(chunk, filter_obs, mjd_obs, variability_cache, dflux)

    dummy_sed = Sed()
    lsst_bp = BandpassDict.loadTotalBandpassesFromFiles()
    flux_q = np.zeros((6,n_obj), dtype=float)
    flux_coadd = np.zeros((6,n_obj), dtype=float)
    mag_coadd = np.zeros((6,n_obj), dtype=float)
    snr_coadd = np.zeros((6,n_obj), dtype=float)
    snr_single = {}
    snr_single_mag_grid = np.arange(14.0, 30.0, 0.05)

    phot_params_single = PhotometricParameters(nexp=1,
                                               exptime=30.0)

    t_start_snr = time.time()
    photometry_mask = np.zeros((n_obj, n_t), dtype=bool)
    photometry_mask_1d = np.zeros(n_obj, dtype=bool)

    for i_bp, bp in enumerate('ugrizy'):
        valid_obs = np.where(filter_obs == i_bp)
        phot_params_coadd = PhotometricParameters(nexp=1,
                                                  exptime=30.0*coadd_visits[bp])

        flux_q[i_bp] = dummy_sed.fluxFromMag(chunk['%smag' % bp])
        dflux_sub = dflux[:,valid_obs[0]]
        assert dflux_sub.shape == (n_obj, len(valid_obs[0]))
        dflux_mean = np.mean(dflux_sub, axis=1)
        assert dflux_mean.shape==(n_obj,)
        flux_coadd[i_bp] = flux_q[i_bp]+dflux_mean
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
    snr_arr = np.zeros((n_obj, n_t), dtype=float)

    for i_bp, bp in enumerate('ugrizy'):
        valid_obs = np.where(filter_obs==i_bp)
        if len(valid_obs[0])==0:
            continue
        n_bp = len(valid_obs[0])
        dflux_bp = dflux[:,valid_obs[0]]
        flux0_arr = flux_q[i_bp]
        flux_tot = flux0_arr[:,None] + dflux_bp
        mag_tot = dummy_sed.magFromFlux(flux_tot)
        snr_single_val = np.interp(mag_tot,
                                  snr_single_mag_grid,
                                  snr_single[bp])

        noise_coadd = flux_coadd[i_bp]/snr_coadd[i_bp]
        noise_single = flux_tot/snr_single_val
        noise = np.sqrt(noise_coadd[:,None]**2+noise_single**2)
        dflux_thresh = 5.0*noise
        dflux_bp = np.abs(dflux_bp)
        detected = (dflux_bp>=dflux_thresh)
        snr_arr[:,valid_obs[0]] = dflux_bp/noise
        for i_obj in range(n_obj):
            if detected[i_obj].any():
                photometry_mask_1d[i_obj] = True
                photometry_mask[i_obj,valid_obs[0]] = detected[i_obj]

    print('doing first pass of photometry took %e hrs'
    % ((time.time()-t_start_obj)/3600.0))

    t_before_chip = time.time()
    chip_mask = apply_focal_plane(chunk['ra'], chunk['decl'],
                                  photometry_mask_1d, obs_md_list,
                                  filter_obs, proper_chip)
    duration = (time.time()-t_before_chip)/3600.0

    for i_obj in range(n_obj):
        if photometry_mask_1d[i_obj]:
            detected = photometry_mask[i_obj,:] & chip_mask[i_obj,:]
            if detected.any():
                unq = chunk['simobjid'][i_obj]
                first_dex = np.where(detected)[0].min()
                out_data[unq] = (mjd_obs[first_dex],
                                 snr_arr[i_obj, first_dex])
                if detected[0]:
                    ct_first += 1
                else:
                    ct_at_all += 1

    #print('%d tot %d first %d at all %d ' %
    #(os.getpid(),ct_tot, ct_first, ct_at_all))

if __name__ == "__main__":

    t_start = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--out_name', type=str, default=None)
    parser.add_argument('--circular_fov', default=False,
                        action='store_true')
    parser.add_argument('--q_chunk_size', type=int, default=20000,
                        help='number of galaxies to query from '
                             'database at once (default 2*10**4)')
    parser.add_argument('--p_chunk_size', type=int, default=10000,
                        help='number of galaxies to process at once '
                             '(default 10**4)')
    parser.add_argument('--n_threads', type=int, default=30)
    parser.add_argument('--htmid', type=int, nargs='+', default=None)

    args = parser.parse_args()
    proper_chip = not args.circular_fov
    assert args.out_name is not None
    assert args.htmid is not None
    if isinstance(args.htmid, numbers.Number):
        htmid_list = [args.htmid]
    else:
        htmid_list = args.htmid


    variability_cache = create_variability_cache()
    plc = ParametrizedLightCurveMixin()
    plc.load_parametrized_light_curves(variability_cache=variability_cache)

    variability_cache['rrly_map'] = {}
    lc_dir = os.path.join(os.environ['SIMS_SED_LIBRARY_DIR'], 'rrly_lc')
    assert os.path.isdir(lc_dir)
    with open('data/rrly_lc_map.txt', 'r') as in_file:
        for line in in_file:
            p = line.strip().split(';')
            variability_cache['rrly_map'][int(p[0])] = os.path.join(lc_dir, p[1])

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

    # get a crowded field (galactic center)
    #from lsst.sims.utils import equatorialFromGalactic
    #ra, dec = equatorialFromGalactic(0.0, 0.0)
    #print('ra dec ',ra,dec)
    #htmid_query = htm.findHtmid(ra, dec, 6)

    try:
        star_db = LocalStarCatalogObj(database='LSST',
                                      host='epyc.astro.washington.edu',
                                      port=1433,
                                      driver='mssql+pymssql')
    except:
        star_db = LocalStarCatalogObj(database='LSST',
                                      host='localhost',
                                      port=51432,
                                      driver='mssql+pymssql')


    obs_param_name = 'data/obs_params.h5'
    obs_params = h5py.File(obs_param_name, 'r')

    mgr = multiprocessing.Manager()
    out_data = mgr.dict()

    for htmid_query in htmid_list:
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

        col_names = ['ra', 'decl',
                     'umag', 'gmag', 'rmag',
                     'imag', 'zmag', 'ymag',
                     'lc_id', 't0', 'var_type',
                     'ebv', 'parallax',
                     'simobjid']

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

        constraint = 'isvar=1 '
        htmid_21_min = htmid_query<<2*(21-query_level)
        htmid_21_max = (htmid_query+1)<<2*(21-query_level)
        constraint += 'AND htmid>=%d AND htmid<%d' % (htmid_21_min, htmid_21_max)

        data_iter = star_db.query_columns(col_names, obs_metadata=obs_query,
                                          chunk_size=args.q_chunk_size,
                                          constraint=constraint)

        p_list = []
        i_chunk = 0
        to_concatenate = []
        n_tot = 0
        n_processed = 0
        for chunk in data_iter:
            htmid_found = htm.findHtmid(chunk['ra'],
                                        chunk['decl'],
                                        query_level)

            valid = np.where(htmid_found==htmid_query)
            if len(valid[0]) == 0:
                continue

            chunk = chunk[valid]
            n_tot += len(chunk)
            #print('n_tot %e ' % n_tot)

            #i_chunk += 1
            #process_stellar_chunk(chunk, filter_obs, mjd_obs, m5_obs, coadd_m5,
            #                      obs_md_list, proper_chip, variability_cache,
            #                      out_data)

            #if i_chunk>3:
            #    break
            #continue

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

                p = multiprocessing.Process(target=process_stellar_chunk,
                                            args=(sub_chunk, filter_obs, mjd_obs,
                                                  m5_obs, coadd_m5, obs_md_list,
                                                  proper_chip, variability_cache,
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
                p = multiprocessing.Process(target=process_stellar_chunk,
                                            args=(sub_chunk, filter_obs, mjd_obs,
                                                  m5_obs, coadd_m5, obs_md_list,
                                                  proper_chip, variability_cache,
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

    out_data_final = {}
    for name in out_data.keys():
        out_data_final[name] = out_data[name]

    print('n_lc %d' % len(out_data_final))
    with open(args.out_name, 'wb') as out_file:
        pickle.dump(out_data_final, out_file)


    #with h5py.File(out_name, 'w') as out_file:
    #    print('n_lc %d' % len(out_data))
    #    for name in out_data.keys():
    #        out_file.create_dataset('%d' % name, data=out_data[name])

    print('that took %e hrs' % ((time.time()-t_start)/3600.0))
    print('shld %d processed %d' % (n_tot, n_processed))
    obs_params.close()
