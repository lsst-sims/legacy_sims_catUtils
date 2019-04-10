import h5py
import pickle
import numpy as np
import os

import time

from lsst.sims.utils import htmModule as htm
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels.LocalGalaxyModels import LocalGalaxyTileObj
from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels

if __name__ == "__main__":

    htmid_map_name = 'data/htmid_to_obs_map.pickle'
    assert os.path.isfile(htmid_map_name)
    with open(htmid_map_name, 'rb') as in_file:
        htmid_to_obs = pickle.load(in_file)

    threshold = 1000
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

    gal_db = LocalGalaxyTileObj(database='LSST',
                                host='epyc.astro.washington.edu',
                                port=1433,
                                driver='mssql+pymssql')


    chunk_size = 10000
    data_iter = gal_db.query_columns(col_names, obs_metadata=obs_query,
                                     chunk_size=chunk_size,
                                     constraint='isagn=1')

    for chunk in data_iter:
        htmid_found = htm.findHtmid(chunk['ra'],
                                    chunk['dec'],
                                    query_level)

        valid = np.where(htmid_found==htmid_query)
        print('valid %d -- %e %e' % (len(valid[0]), chunk['ra'].min(),
                                     chunk['ra'].max()))

    obs_params.close()
