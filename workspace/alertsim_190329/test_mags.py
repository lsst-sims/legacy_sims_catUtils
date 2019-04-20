from lsst.sims.catUtils.baseCatalogModels.LocalStarModels import LocalStarCatalogObj
from lsst.sims.photUtils import cache_LSST_seds
from lsst.sims.photUtils import Sed, BandpassDict
from lsst.sims.photUtils import getImsimFluxNorm
from lsst.sims.utils import ObservationMetaData
import os
import numpy as np
from lsst.sims.utils import defaultSpecMap

db = LocalStarCatalogObj(database='LSST',
                            host='epyc.astro.washington.edu',
                            port=1433,
                            driver='mssql+pymssql')

cache_LSST_seds(wavelen_min=100.0, wavelen_max=1500.0)

lsst_bp = BandpassDict.loadTotalBandpassesFromFiles()

col_names = ['sedFilename', 'magNorm',
             'ebv',
             'umag', 'gmag', 'rmag',
             'imag', 'zmag', 'ymag']

obs = ObservationMetaData(pointingRA=22.0, pointingDec=-12.0,
                          boundType='circle', boundLength=0.1)

data_iter = db.query_columns(col_names, obs_metadata=obs,
                             chunk_size=100,
                             constraint="var_type=2")


for chunk in data_iter:
    for star in chunk:
        sed = Sed()
        sed.readSED_flambda(os.path.join(os.environ['SIMS_SED_LIBRARY_DIR'],
                                         defaultSpecMap[star['sedFilename']]))

        fnorm = getImsimFluxNorm(sed, star['magNorm'])
        sed.multiplyFluxNorm(fnorm)
        mags = lsst_bp.magListForSed(sed)

        du_nodust = np.abs(mags[0]-star['umag'])
        dg_nodust = np.abs(mags[1]-star['gmag'])
        dr_nodust = np.abs(mags[2]-star['rmag'])
        di_nodust = np.abs(mags[3]-star['imag'])
        dz_nodust = np.abs(mags[4]-star['zmag'])
        dy_nodust = np.abs(mags[5]-star['ymag'])


        a_x, b_x = sed.setupCCM_ab()
        sed.addDust(a_x, b_x, R_v=3.1, ebv=star['ebv'])
        mags = lsst_bp.magListForSed(sed)
        du = np.abs(mags[0]-star['umag'])
        dg = np.abs(mags[1]-star['gmag'])
        dr = np.abs(mags[2]-star['rmag'])
        di = np.abs(mags[3]-star['imag'])
        dz = np.abs(mags[4]-star['zmag'])
        dy = np.abs(mags[5]-star['ymag'])

        print('%.2e %.2e %2e %.2e %.2e %.2e -- %.2e %.2e %2e %.2e %.2e %.2e' %
        (du_nodust, dg_nodust, dr_nodust, di_nodust,
         dz_nodust, dy_nodust,
         du,dg,dr,di,dz,dy))
    break
