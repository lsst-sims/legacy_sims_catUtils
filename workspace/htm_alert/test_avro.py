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

avro_gen = AvroGenerator(obs_list, n_proc_max=4)

from lsst.sims.catUtils.baseCatalogModels import StarObj

def query_htmid(avro_gen, htmid_list, output_list):

    db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                 port=1433, driver='mssql+pymssql')

    for htmid in htmid_list:
        ct = avro_gen.alert_data_from_htmid(htmid, db)
        output_list.append(ct)

import time
total_obs = len(obs_list)

import multiprocessing as mproc
t_start = time.time()
n_proc = 4

mgr = mproc.Manager()
out_list = mgr.list()

process_list = []
htmid_list=[]
for i_list in range(n_proc):
    htmid_list.append([])
i_list = 0
for htmid in avro_gen.htmid_list[:8]:
    htmid_list[i_list].append(htmid)
    i_list += 1
    if i_list>=n_proc:
        i_list = 0

for i_proc in range(n_proc):
    p = mproc.Process(target=query_htmid,
                      args=(avro_gen, htmid_list[i_proc], out_list))
    p.start()
    process_list.append(p)

for p in process_list:
    p.join()

so_far = 0
for ii in out_list:
    so_far += ii
elapsed = time.time()-t_start
print('%d took %.2e hours' % (so_far, elapsed/3600.0))
print('total will take %.2e hours' % (total_obs*elapsed/(3600.0*so_far)))
