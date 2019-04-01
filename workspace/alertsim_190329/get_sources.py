import os
import h5py
import pickle
import hashlib
import time
import numpy as np

import lsst.sims.utils as sims_utils
import lsst.sims.utils.htmModule as htm
import lsst.sims.coordUtils as sims_coord
from lsst.sims.catUtils.baseCatalogModels import LocalStarCatalogObj


def process_chunk_of_stars(data, obs_md_data):
    kplr_dexes = np.where(np.char.find(data['varParamStr'], 'kplr')>=0)
    mlt_dexes = np.where(np.char.find(data['varParamStr'], 'MLT')>=0)
    rrly_dexes = np.where(np.char.find(data['varParamStr'], 'applyRRly')>=0)

    return len(kplr_dexes[0]), len(mlt_dexes[0]), len(rrly_dexes[0])


if __name__ == "__main__":

    star_db = LocalStarCatalogObj(database='LSST',
                                  host='localhost',
                                  port=51432,
                                  driver='mssql+pymssql')

    col_names = ['htmID', 'ra', 'decl',
                 'umag', 'gmag', 'rmag',
                 'imag', 'zmag', 'ymag',
                 'ebv', 'varParamStr']

    star_dtype = np.dtype([('htmID', int), ('simobjid', int),
                           ('ra', float), ('dec', float),
                           ('umag', float), ('gmag', float), ('rmag', float),
                           ('imag', float), ('zmag', float), ('ymag', float),
                           ('ebv', float), ('varParamStr', str, 256)])

    obs_fname = "data/obs_params.h5"
    assert os.path.isfile(obs_fname)
    htmid_to_obs_fname = "data/htmid_to_obs_map.pickle"
    assert os.path.isfile(htmid_to_obs_fname)

    with open(htmid_to_obs_fname, 'rb') as in_file:
        htmid_to_obs = pickle.load(in_file)

    md5_hasher = hashlib.md5()
    with open(obs_fname, 'rb') as in_file:
        for chunk in iter(lambda: in_file.read(4096), b""):
            md5_hasher.update(chunk)
    md5_sum = md5_hasher.hexdigest()

    assert md5_sum == htmid_to_obs['md5_sum']

    htmid_of_interest = 9484
    htmid_level_of_interest = htm.levelFromHtmid(htmid_of_interest)
    obs_list = np.array(htmid_to_obs[htmid_of_interest])
    obs_dex_list = obs_list-1
    print(len(obs_list))

    with h5py.File(obs_fname, 'r') as obs_file:
        shld_be = 1 + np.arange(len(obs_file['obsHistID'].value), dtype=int)
        np.testing.assert_array_equal(shld_be, obs_file['obsHistID'].value)
        np.testing.assert_array_equal(obs_file['obsHistID'].value[obs_dex_list],
                                      obs_list)

        obs_subset = {}
        for field_name in obs_file.keys():
            obs_subset[field_name] = obs_file[field_name].value[obs_dex_list]

        trixel_of_interest = htm.trixelFromHtmid(htmid_of_interest)
        trixel_ra, trixel_dec = trixel_of_interest.get_center()
        trixel_radius = 1.05*trixel_of_interest.get_radius()
        trixel_obs = sims_utils.ObservationMetaData(pointingRA=trixel_ra,
                                                    pointingDec=trixel_dec,
                                                    boundType='circle',
                                                    boundLength=trixel_radius)

        star_iter = star_db.query_columns(col_names,
                                          obs_metadata=trixel_obs,
                                          chunk_size=100000)

        ct_tot = 0
        ct_kplr = 0
        ct_mlt = 0
        ct_rrly = 0
        t_start = time.time()
        for chunk in star_iter:
            chunk_htmid = chunk['htmID']>>(2*(21-htmid_level_of_interest))
            valid = np.where(chunk_htmid==htmid_of_interest)
            if len(valid[0]) == 0:
                continue
            chunk = chunk[valid]
            i_k, i_m, i_r = process_chunk_of_stars(chunk, obs_file)
            ct_tot += len(chunk)
            ct_kplr += i_k
            ct_mlt += i_m
            ct_rrly += i_r
        print(ct_tot)
        print(ct_kplr+ct_mlt+ct_rrly)
        print(ct_kplr,ct_mlt,ct_rrly)
        print(time.time()-t_start)
