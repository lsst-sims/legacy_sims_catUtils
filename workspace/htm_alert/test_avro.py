from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AvroGenerator

import os
from lsst.utils import getPackageDir

opsim_db = os.path.join('/Users', 'danielsf', 'physics', 'lsst_150412',
                        'Development', 'garage', 'OpSimData',
                        'minion_1016_sqlite.db')

if not os.path.exists(opsim_db):
    opsim_db = os.path.join('/local', 'lsst', 'danielsf', 'OpSimData',
                            'minion_1016_sqlite.db')

obs_gen = ObservationMetaDataGenerator(opsim_db, driver='sqlite')

obs_list = obs_gen.getObservationMetaData(night=(0,40))

print('%d obs' % len(obs_list))

avro_gen = AvroGenerator(obs_list)

from lsst.sims.catUtils.baseCatalogModels import StarObj
db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
             port=1433, driver='mssql+pymssql')

import time
t_start = time.time()
avro_gen.alerts_from_db(db)
print('total took %e hours' % ((time.time()-t_start)/3600.0))
