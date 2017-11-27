import os
import numpy as np
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST
from lsst.sims.utils import _pupilCoordsFromRaDec, ObservationMetaData
from lsst.sims.utils import trixelFromHtmid
from lsst.sims.catUtils.utils import StellarAlertDBObj
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator

htmid = 2682
trixel = trixelFromHtmid(htmid)

ra, dec = trixel.get_center()
ra = ra%360.0
radius = trixel.get_radius()
cosdec = np.cos(np.radians(dec))

opsim_db = os.path.join('/Users', 'danielsf', 'physics', 'lsst_150412',
                        'Development', 'garage', 'OpSimData',
                        'minion_1016_sqlite.db')

if not os.path.exists(opsim_db):
    opsim_db = os.path.join('/local', 'lsst', 'danielsf', 'OpSimData',
                            'minion_1016_sqlite.db')

obs_gen = ObservationMetaDataGenerator(opsim_db, driver='sqlite')



obs_list = obs_gen.getObservationMetaData(fieldRA=(ra-radius/cosdec, ra+radius/cosdec),
                                          fieldDec=(dec-radius, dec+radius))

db = StellarAlertDBObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                       port=1433, driver='mssql+pymssql')

colnames = ['raJ2000', 'decJ2000',
            'properMotionRa', 'properMotionDec',
            'radialVelocity', 'parallax']

data_iter = db.query_columns_htmid(htmid=htmid, chunk_size=1000)

for i_chunk, chunk in enumerate(data_iter):
    if i_chunk%100 ==0:
        print('    name %d' % i_chunk)
    for obs in obs_list:
        xp, yp = _pupilCoordsFromRaDec(chunk['raJ2000'], chunk['decJ2000'],
                                       pm_ra=chunk['properMotionRa'],
                                       pm_dec=chunk['properMotionDec'],
                                       v_rad=chunk['radialVelocity'],
                                       parallax=chunk['parallax'],
                                       obs_metadata=obs)

        chip_name_list = chipNameFromPupilCoordsLSST(xp, yp)
