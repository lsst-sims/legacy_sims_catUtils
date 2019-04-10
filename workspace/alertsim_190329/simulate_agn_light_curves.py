import h5py
import pickle
import numpy as np
import os

import time

from lsst.sims.utils import angularSeparation
from lsst.sims.utils import htmModule as htm
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels.LocalGalaxyModels import LocalGalaxyTileObj
from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels

import multiprocessing

def process_agn_chunk(chunk, filter_obs):
    agn_model = ExtraGalacticVariabilityModels()

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
                              redshift=chunk['redshift']).transpose(1,2,0)

    for ii in range(len(chunk)):
        dmag_obs = np.array([dmag[ii][tt][filter_obs[tt]]
                             for tt in range(len(filter_obs))])

if __name__ == "__main__":

    fov_radius = 1.75

    htmid_map_name = 'data/htmid_to_obs_map.pickle'
    assert os.path.isfile(htmid_map_name)
    with open(htmid_map_name, 'rb') as in_file:
        htmid_to_obs = pickle.load(in_file)

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
                 'ra', 'dec']

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

    print('%d time steps' % len(filter_obs))

    chunk_size = 10000
    data_iter = gal_db.query_columns(col_names, obs_metadata=obs_query,
                                     chunk_size=chunk_size,
                                     constraint='isagn=1')

    t_start = time.time()
    p_list = []
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

        p = multiprocessing.Process(target=process_agn_chunk,
                                    args=(chunk, filter_obs))
        p.start()
        p_list.append(p)
        if len(p_list)>30:
            for p in p_list:
                p.join()
            p_list = []


    for p in p_list:
        p.join()

    print('that took %e hrs' % ((time.time()-t_start)/3600.0))
    obs_params.close()
