from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AlertDataGenerator
from lsst.sims.catUtils.utils import AlertAgnVariabilityCatalog

import os
from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up

import time
import gc
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--n_proc', type=int, default=1)

args = parser.parse_args()

opsim_db = os.path.join('/Users', 'danielsf', 'physics', 'lsst_150412',
                        'Development', 'garage', 'OpSimData',
                        'minion_1016_sqlite.db')

if not os.path.exists(opsim_db):
    opsim_db = os.path.join('/local', 'lsst', 'danielsf', 'OpSimData',
                            'minion_1016_sqlite.db')

obs_gen = ObservationMetaDataGenerator(opsim_db, driver='sqlite')

from lsst.sims.utils import trixelFromHtmid

trixel = trixelFromHtmid(2661)
ra0, dec0 = trixel.get_center()
ra0 = ra0%360.0

obs_list = obs_gen.getObservationMetaData(night=(30,61))


del obs_gen
sims_clean_up()
gc.collect()

alert_gen = AlertDataGenerator(n_proc_max=1)
alert_gen.subdivide_obs(obs_list)

from lsst.sims.catUtils.utils import AgnAlertDBObj

db = AgnAlertDBObj(database='LSSTCATSIM',
                   host='fatboy.phys.washington.edu',
                   port=1433, driver='mssql+pymssql',
                   cache_connection=False)


htmid_max=-1
n_max=0
for htmid in alert_gen.htmid_list:
    n_obs = len(alert_gen._htmid_dict[htmid])
    if n_obs>n_max:
        n_max = n_obs
        htmid_max = htmid


import time
t_start = time.time()
alert_gen.alert_data_from_htmid(htmid_max, db,
                                chunk_size=1000,
                                write_every=100000,
                                photometry_class=AlertAgnVariabilityCatalog)

elapsed = time.time()-t_start
print("that took %.2e hours" % (elapsed/3600.0))
exit()

def query_htmid(alert_gen, htmid_list, output_list):

    db = AgnAlertDBObj(database='LSSTCATSIM',
                       host='fatboy.phys.washington.edu',
                       port=1433, driver='mssql+pymssql',
                       cache_connection=False)

    t_start = time.time()
    for i_htmid, htmid in enumerate(htmid_list):
        ct = alert_gen.alert_data_from_htmid(htmid, db)
        output_list.append(ct)
        elapsed = time.time()-t_start
        print('%d took %.2e hours; total %.2e' %
        (i_htmid+1,elapsed/3600.0,len(htmid_list)*elapsed/(3600.0*(i_htmid+1))))

total_obs = len(obs_list)

import multiprocessing as mproc
t_start = time.time()

mgr = mproc.Manager()
out_list = mgr.list()

process_list = []
htmid_list=[]
for i_list in range(args.n_proc):
    htmid_list.append([])
i_list = 0

n_total_htmid = len(alert_gen.htmid_list)
n_htmid = 0

for htmid in alert_gen.htmid_list[:args.n_proc*10]:
    n_htmid += 1
    htmid_list[i_list].append(htmid)
    i_list += 1
    if i_list>=args.n_proc:
        i_list = 0

for i_proc in range(args.n_proc):
    p = mproc.Process(target=query_htmid,
                      args=(alert_gen, htmid_list[i_proc], out_list))
    p.start()
    process_list.append(p)

for p in process_list:
    p.join()

so_far = 0
for ii in out_list:
    so_far += ii
elapsed = time.time()-t_start
print('%d took %.2e hours' % (n_htmid, elapsed/3600.0))
print('total will take %.2e hours' % (n_total_htmid*elapsed/(3600.0*n_htmid)))
