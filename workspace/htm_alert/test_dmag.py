import os
import numpy as np
from lsst.sims.utils import _pupilCoordsFromRaDec, ObservationMetaData
from lsst.sims.utils import trixelFromHtmid
from lsst.sims.catUtils.utils import StellarAlertDBObj
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AlertStellarVariabilityCatalog

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

print('n_obs %d' % len(obs_list))

tai_arr = np.zeros(len(obs_list), dtype=float)
for i_obs in range(len(obs_list)):
    tai_arr[i_obs] = obs_list[i_obs].mjd.TAI

db = StellarAlertDBObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                       port=1433, driver='mssql+pymssql')

colnames = ['raJ2000', 'decJ2000','varParamStr', 'parallax', 'ebv']

data_iter = db.query_columns_htmid(htmid=htmid, chunk_size=1000)
phot_cat = AlertStellarVariabilityCatalog(db, obs_list[0],
                                         column_outputs=['lsst_u',
                                                         'lsst_g',
                                                         'lsst_r',
                                                         'lsst_i',
                                                         'lsst_z',
                                                         'lsst_y'])
phot_cat.load_parametrized_light_curves()

for i_chunk, chunk in enumerate(data_iter):
    if i_chunk%100 ==0:
        print('    dmag %d' % i_chunk)

    phot_cat._set_current_chunk(chunk)
    dmag_arr = phot_cat.applyVariability(chunk['varParamStr'],
                                         expmjd=tai_arr)
