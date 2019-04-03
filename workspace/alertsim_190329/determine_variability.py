import os
import numpy as np
import h5py

from lsst.sims.catalogs.db import DBObject

import argparse
import time


def parse_mlt(chunk):
    if not hasattr(parse_mlt, 'lc_id_lookup'):
        parse_mlt.lc_id_lookup = {}
        for ii in '0123':
            iii = int(ii)
            parse_mlt.lc_id_lookup[10+iii] = 'early_inactive_%s' % ii
            parse_mlt.lc_id_lookup[20+iii] = 'early_active_%s' % ii
            parse_mlt.lc_id_lookup[30+iii] = 'mid_inactive_%s' % ii
            parse_mlt.lc_id_lookup[40+iii] = 'mid_active_%s' % ii
            parse_mlt.lc_id_lookup[50+iii] = 'late_active_%s' % ii

        parse_mlt.dflux_lookup = {}
        with h5py.File('data/dflux_SNR5_lookup.h5','r') as snr_dflux_file:
            for kk in snr_dflux_file.keys():
                parse_mlt.dflux_lookup[kk] = snr_dflux_file[kk].value

    mlt_dflux_name = 'data/mlt_dflux_lookup.h5'
    assert os.path.isfile(mlt_dflux_name)
    is_var = np.zeros(len(chunk), dtype=int)
    with h5py.File(mlt_dflux_name, 'r') as mlt_dflux_file:
        px_grid = mlt_dflux_file['parallax_grid(mas)'].value
        ebv_01 = np.round(chunk['ebv'], decimals=2)
        unq_ebv_01 = np.unique(ebv_01)
        unq_lc_id = np.unique(chunk['lc_id'])

        # loop over all combinations of ebv and lc_id
        for ebv in unq_ebv_01:
            for lc_id in unq_lc_id:

                # figure out which stars we are evaluating now
                considering = np.where(np.logical_and(np.abs(ebv_01-ebv)<0.0001,
                                                      chunk['lc_id']==lc_id))

                if len(considering[0]) == 0:
                    continue

                # loop over bandpasses, figuring out what the dflux threshold
                # is for each of these stars in each bandpass
                for bp in 'ugrizy':
                    dflux_threshold = np.interp(chunk[bp][considering],
                                                parse_mlt.dflux_lookup['mag_grid'],
                                                parse_mlt.dflux_lookup['%s_dflux' % bp])


                    # interpolate the actual dflux value of these stars at these
                    # parallaxes
                    mlt_name = '%s_ebv_%.2f_%s' % (parse_mlt.lc_id_lookup[lc_id],
                                               ebv, bp)
                    dflux_max = np.interp(chunk['parallax'][considering],
                                          px_grid,
                                          mlt_dflux_file[mlt_name])

                    # any star that might cross the threshold gets marked
                    # as is_var == 1
                    local_is_var = (dflux_max>=dflux_threshold)
                    local_is_var_dexes = considering[0][local_is_var]
                    is_var[local_is_var_dexes] = 1

    return is_var


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--partition', type=str, default=None,
                        help='htmid tag of stars_partition_* table to run')

    args = parser.parse_args()

    try:
        db = DBObject(database='LSST',
                      port=1433,
                      host='epyc.astro.washington.edu',
                      driver='mssql+pymssql')
    except:
        db = DBObject(database='LSST',
                      port=51432,
                      host='localhost',
                      driver='mssql+pymssql')


    table_name = 'stars_partition_%s' % args.partition
    out_name = 'data/isvar_lookup_%s.txt' % args.partition
    if os.path.isfile(out_name):
        os.unlink(out_name)
        #raise RuntimeError("\n%s\nexists\n" % out_name)

    query = "SELECT top 1000000 "
    query += "htmid, simobjid, umag, gmag, rmag, imag, zmag, ymag, "
    query += "parallax, ebv, lc_id, var_type "
    query += "FROM %s" % table_name

    dtype = np.dtype([('htmid', int), ('simobjid', int),
                      ('u', float), ('g', float), ('r', float),
                      ('i', float), ('z', float), ('y', float),
                      ('parallax', float), ('ebv', float),
                      ('lc_id', int), ('var_type', int)])

    t_start = time.time()
    data_iter = db.get_arbitrary_chunk_iterator(query, dtype=dtype,
                                                chunk_size=100000)

    ct_var = 0
    ct_mlt = 0
    with open(out_name, 'w') as out_file:
        for chunk in data_iter:
            is_kplr = np.where(chunk['var_type'] == 1)
            is_mlt = np.where(chunk['var_type'] == 2)

            #print('is_kplr %d' % len(is_kplr[0]))
            #print('is_mlt %d' % len(is_mlt[0]))
            is_var = parse_mlt(chunk[is_mlt])
            ct_var += is_var.sum()
            ct_mlt += len(is_mlt[0])

    print('is_var %d of %d' % (ct_var, ct_mlt))
    print('took %e' % (time.time()-t_start))
