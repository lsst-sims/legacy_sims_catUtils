import h5py
import pickle
import numpy as np
import os

import time

from lsst.sims.utils import angularSeparation
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed
from lsst.sims.photUtils import SignalToNoise as SNR
from lsst.sims.utils import htmModule as htm
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels.LocalGalaxyModels import LocalGalaxyTileObj
from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels

import multiprocessing

def process_agn_chunk(chunk, filter_obs, mjd_obs, m5_obs,
                      coadd_m5, out_data):

    ct_first = 0
    ct_at_all = 0
    ct_tot = 0

    agn_model = ExtraGalacticVariabilityModels()

    coadd_visits = {}
    coadd_visits['u'] = 6
    coadd_visits['g'] = 8
    coadd_visits['r'] = 18
    coadd_visits['i'] = 18
    coadd_visits['z'] = 16
    coadd_visits['y'] = 16

    n_t = len(filter_obs)
    n_obj = len(chunk)

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

    dmag = agn_model.applyAgn(np.where(np.array([True]*len(chunk))),
                              params, mjd_obs,
                              redshift=chunk['redshift'])

    dmag_mean = np.mean(dmag, axis=2)

    dummy_sed = Sed()
    lsst_bp = BandpassDict.loadTotalBandpassesFromFiles()
    flux_gal = np.zeros((6,n_obj), dtype=float)
    flux_agn_q = np.zeros((6,n_obj), dtype=float)
    flux_coadd = np.zeros((6,n_obj), dtype=float)
    mag_coadd = np.zeros((6,n_obj), dtype=float)
    snr_coadd = np.zeros((6,n_obj), dtype=float)
    snr_single = {}
    snr_single_mag_grid = {}

    phot_params_single = PhotometricParameters(nexp=1,
                                               exptime=30.0)

    t_start_snr = time.time()
    for i_bp, bp in enumerate('ugrizy'):
        phot_params_coadd = PhotometricParameters(nexp=1,
                                                  exptime=30.0*coadd_visits[bp])

        flux_gal[i_bp] = dummy_sed.fluxFromMag(chunk['%s_ab' % bp])
        flux_agn_q[i_bp] = dummy_sed.fluxFromMag(chunk['%s_ab' % bp] +
                                               dmag_mean[i_bp,:])
        flux_coadd[i_bp] = flux_gal[i_bp]+flux_agn_q[i_bp]
        mag_coadd[i_bp] = dummy_sed.magFromFlux(flux_coadd[i_bp])

        (snr_coadd[i_bp],
         gamma) = SNR.calcSNR_m5(mag_coadd[i_bp],
                                 lsst_bp[bp],
                                 coadd_m5[bp],
                                 phot_params_coadd)

        snr_single_mag_grid[bp] = np.arange(mag_coadd[i_bp].min()-2.0,
                                            mag_coadd[i_bp].max()+2.0,
                                            0.05)
        snr_single[bp] = {}
        for i_t in range(n_t):
            (snr_single[bp][i_t],
             gamma) = SNR.calcSNR_m5(snr_single_mag_grid[bp],
                                     lsst_bp[bp],
                                     m5_obs[i_t],
                                     phot_params_single)

    print('got all snr in %e' % (time.time()-t_start_snr))


    snr_arr = []
    t_start_obj = time.time()
    for i_obj in range(n_obj):
        if i_obj>0 and i_obj%100==0:
            duration = (time.time()-t_start_obj)/3600.0
            print('    %d in %e hrs' % (i_obj,duration))
        ct_tot += 1
        unq = chunk['galtileid'][i_obj]
        first_detection = None

        bp_arr = list(['ugrizy'[filter_obs[i_t]] for i_t in range(n_t)])
        mag0_arr = np.array([chunk['AGNLSST%s' % bp][i_obj] for bp in bp_arr])
        dmag_arr = np.array([dmag[filter_obs[i_t]][i_obj][i_t]
                             for i_t in range(n_t)])

        agn_flux_tot = dummy_sed.fluxFromMag(mag0_arr+dmag_arr)
        q_flux = np.array([flux_agn_q[ii][i_obj] for ii in filter_obs])
        agn_dflux = np.abs(agn_flux_tot-q_flux)
        flux_tot = np.array([flux_gal[ii][i_obj] for ii in filter_obs])
        flux_tot += agn_flux_tot
        mag_tot = dummy_sed.magFromFlux(flux_tot)

        for i_t, i_bp in zip(range(n_t), filter_obs):
            bp = 'ugrizy'[i_bp]

            snr_single_val = np.interp(mag_tot[i_t],
                                       snr_single_mag_grid[bp],
                                       snr_single[bp][i_t])

            noise_coadd = flux_coadd[i_bp][i_obj]/snr_coadd[i_bp][i_obj]
            noise_single = flux_tot[i_t]/snr_single_val
            noise = np.sqrt(noise_coadd**2+noise_single**2)
            dflux_thresh = 5.0*noise
            snr_arr.append(agn_dflux[i_t]/dflux_thresh)
            if agn_dflux[i_t]>=dflux_thresh:
                out_data[unq] = mjd_obs[i_t]
                if i_t==0:
                    ct_first += 1
                else:
                    ct_at_all += 1
                break

    snr_arr = np.array(snr_arr)
    print('%d tot %d first %d at all %d -- %e %e' %
    (os.getpid(),ct_tot, ct_first, ct_at_all,snr_arr.min(),np.mean(snr_arr)))

if __name__ == "__main__":

    fov_radius = 1.75

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

    threshold = 5000
    for kk in htmid_to_obs:
        n_obs = len(htmid_to_obs[kk])
        if n_obs>threshold and n_obs<2*threshold:
            htmid_query = kk
            break

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

    obs_param_name = 'data/obs_params.h5'
    obs_params = h5py.File(obs_param_name, 'r')

    assert np.diff(obs_params['obsHistID']).min()>0

    gal_db = LocalGalaxyTileObj(database='LSST',
                                host='epyc.astro.washington.edu',
                                port=1433,
                                driver='mssql+pymssql')

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

    print('%d time steps' % len(filter_obs))

    chunk_size = 10000
    data_iter = gal_db.query_columns(col_names, obs_metadata=obs_query,
                                     chunk_size=chunk_size,
                                     constraint='isagn=1')

    t_start = time.time()
    mgr = multiprocessing.Manager()
    out_data = mgr.dict()
    p_list = []
    out_name = '/astro/store/pogo4/danielsf/dummy_agn_lc.h5'
    i_chunk = 0
    for chunk in data_iter:
        htmid_found = htm.findHtmid(chunk['ra'],
                                    chunk['dec'],
                                    query_level)

        valid = np.where(htmid_found==htmid_query)
        if len(valid[0]) == 0:
            continue

        chunk = chunk[valid]

        print('valid %d -- %e %e' % (len(valid[0]), chunk['ra'].min(),
                                     chunk['ra'].max()))

        process_agn_chunk(chunk, filter_obs, mjd_obs, m5_obs, coadd_m5,
                          out_data)

        i_chunk += 1
        if i_chunk >= 5:
            break

        # multiprocessing code
        # p = multiprocessing.Process(target=process_agn_chunk,
        #                            args=(chunk, filter_obs, mjd_obs,
        #                                  m5_obs, coadd_m5, out_data))
        # p.start()
        # p_list.append(p)
        if len(p_list)>30:
            for p in p_list:
                p.join()
            p_list = []


    for p in p_list:
        p.join()

    with h5py.File(out_name, 'w') as out_file:
        print('n_lc %d' % len(out_data))
        for name in out_data.keys():
            out_file.create_dataset('%d' % name, data=out_data[name])

    print('that took %e hrs' % ((time.time()-t_start)/3600.0))
    obs_params.close()
