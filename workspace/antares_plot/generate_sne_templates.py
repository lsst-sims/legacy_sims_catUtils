import numpy as np
from lsst.sims.catUtils.supernovae import SNUniverse, SNObject
from lsst.sims.catalogs.db import DBObject
from lsst.sims.photUtils import CosmologyObject

from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed

import time

import multiprocessing as mproc

def get_sne_from_chunk(chunk, cc, t0, x1, bessell_mag, output_dict):
    for i_sn in range(len(chunk)):
        sn = SNObject()
        sn.set(x1=x1[i_sn], c=cc[i_sn], t0=t0[i_sn], z=chunk['z'][i_sn])
        sn.set_MWebv(0.0)
        sn.source.set_peakmag(bessell_mag[i_sn], band='bessellb', magsys='ab')

        time_arr = np.arange(t0[i_sn]-200.0, t0[i_sn]+200.0, 1.0)
        gmag = np.zeros(len(time_arr), dtype=float)
        imag= np.zeros(len(time_arr), dtype=float)
        for i_t, tt in enumerate(time_arr):
            vv = sn.catsimManyBandMags(tt,bp_dict)
            gmag[i_t] = vv[1]
            imag[i_t] = vv[2]
            #gmag[i_t] = sn.catsimBandFlux(tt,bp_dict['g'])
            #imag[i_t] = sn.catsimBandFlux(tt,bp_dict['i'])

        gmag = dummy_sed.magFromFlux(dummy_sed.fluxFromMag(gmag)+dummy_sed.fluxFromMag(chunk['gmag'][i_sn]))
        imag = dummy_sed.magFromFlux(dummy_sed.fluxFromMag(imag)+dummy_sed.fluxFromMag(chunk['imag'][i_sn]))
        gmag -= chunk['gmag'][i_sn]
        imag -= chunk['imag'][i_sn]

        valid = np.where(np.logical_or(np.abs(gmag)>=0.005, np.abs(imag)>=0.005))

        time_arr = time_arr[valid]
        gmag = gmag[valid]
        imag = imag[valid]

        is_bad = np.where(np.logical_or(np.isnan(gmag),
                          np.logical_or(np.isinf(gmag),
                          np.logical_or(np.isnan(imag), np.isinf(imag)))))

        if len(time_arr)>0 and len(is_bad[0])==0:
            id_val = chunk['galid'][i_sn]
            output_dict['t_%d' % id_val] = time_arr
            output_dict['g_%d' % id_val] = gmag
            output_dict['i_%d' % id_val] = imag


if __name__ == "__main__":

    n_proc = 2
    mgr = mproc.Manager()
    output_dict = mgr.dict()

    rng = np.random.RandomState(812432)

    mjd_min = 59580.0-365.25
    mjd_max = mjd_min+2.0*365.25
    dummy_sed = Sed()

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles(bandpassNames=['u','g','i'])

    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    query = 'SELECT galid, redshift, g_ab, i_ab FROM galaxy ORDER BY galid'
    dtype = np.dtype([('galid', int), ('z', float), ('gmag', float), ('imag', float)])
    data_iter = db.get_arbitrary_chunk_iterator(query, dtype=dtype,
                                                chunk_size=1000)

    cosmo = CosmologyObject()
    snu = SNUniverse()
    snu.suppressDimSN = False
    ct_valid = 0
    hundredyear = 1.0/snu.snFrequency
    print('hundredyear %e' % hundredyear)
    print(-hundredyear/2.0+61406.25,hundredyear/2.0+61406.25)
    print(mjd_min,mjd_max)
    t_start = time.time()
    running_ct = 0
    float_ct = 0
    p_list = []
    for chunk in data_iter:

        snu.numobjs = len(chunk)
        snu.badvalues =-999.0
        t0 = rng.uniform(-hundredyear/2.0+61406.25,
                         hundredyear/2.0+61406.25,
                         size=len(chunk))

        valid = np.where(np.logical_and(t0>mjd_min, t0+100.0<mjd_max))
        chunk = chunk[valid]
        print('chunk len %d' % len(chunk))

        x1 = rng.normal(0.0,1.0, size=len(chunk))
        cc = rng.normal(0.0, 0.1, size=len(chunk))
        mu = cosmo.distanceModulus(redshift=chunk['z'])
        bessell_mag = rng.normal(-19.3, 0.3, size=len(chunk)) + mu

        p = mproc.Process(target=get_sne_from_chunk,
                          args=(chunk, cc, t0, x1, bessell_mag, output_dict))
        p.start()
        p_list.append(p)
        if len(p_list) >= n_proc:
            for p in p_list:
                p.join()
            p_list = []
            print('after batch len %d' % (len(output_dict)/3))


