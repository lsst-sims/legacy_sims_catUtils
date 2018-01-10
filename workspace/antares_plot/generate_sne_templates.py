import numpy as np
from lsst.sims.catUtils.supernovae import SNUniverse, SNObject
from lsst.sims.catalogs.db import DBObject
from lsst.sims.photUtils import CosmologyObject

if __name__ == "__main__":

    rng = np.random.RandomState(812432)

    mjd_min = 59580.0-365.25
    mjd_max = mjd_min+2.0*365.25

    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    query = 'SELECT galid, redshift FROM galaxy'
    dtype = np.dtype([('galid', int), ('z', float)])
    data_iter = db.get_arbitrary_chunk_iterator(query, dtype=dtype,
                                                chunk_size=100000)

    cosmo = CosmologyObject()
    snu = SNUniverse()
    snu.suppressDimSN = False
    ct_valid = 0
    hundredyear = 1.0/snu.snFrequency
    print('hundredyear %e' % hundredyear)
    print(-hundredyear/2.0+61406.25,hundredyear/2.0+61406.25)
    print(mjd_min,mjd_max)
    for chunk in data_iter:
        mu = cosmo.distanceModulus(redshift=chunk['z'])
        snu.numobjs = len(chunk)
        snu.badvalues =-999.0
        t0 = rng.uniform(-hundredyear/2.0+61406.25,
                         hundredyear/2.0+61406.25,
                         size=len(chunk))

        valid = np.where(np.logical_and(t0>mjd_min, t0+100.0<mjd_max))
        ct_valid += len(valid[0])
        print('ct_valid %.2e' % ct_valid)
