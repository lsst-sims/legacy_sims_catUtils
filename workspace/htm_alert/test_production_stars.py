from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AlertDataGenerator

import os
from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up

from lsst.sims.catUtils.utils import StellarAlertDBObj
from lsst.sims.catUtils.utils import AlertStellarVariabilityCatalog


import multiprocessing as mproc

import time
import gc
import argparse


def query_htmid(alert_gen, htmid_list, out_dir, out_prefix,
                log_file_name, lock,
                write_every, chunk_size, dmag_cutoff):

    db = StellarAlertDBObj(database='LSSTCATSIM',
                           host='fatboy.phys.washington.edu',
                           port=1433, driver='mssql+pymssql',
                           cache_connection=False)

    for i_htmid, htmid in enumerate(htmid_list):
        t_start = time.time()
        print('passing in %s' % str(htmid))
        n_rows = alert_gen.alert_data_from_htmid(htmid, db, chunk_size=chunk_size,
                                                 write_every=write_every,
                                                 output_dir=out_dir,
                                                 output_prefix=out_prefix,
                                                 dmag_cutoff=dmag_cutoff,
                                                 photometry_class=AlertStellarVariabilityCatalog,
                                                 lock=lock,
                                                 log_file_name=log_file_name)

        lock.acquire()
        with open(log_file_name, 'a') as out_file:
            elapsed = (time.time()-t_start)/3600.0
            out_file.write('htmid %d nobs %d n_rows %d time %.2e hrs; per_row %.2e -- n_htmid %d of %d\n' %
                           (htmid, alert_gen.n_obs(htmid), n_rows, elapsed,
                            elapsed/n_rows, i_htmid, len(htmid_list)))

        lock.release()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--n_proc', type=int, default=1)
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--out_prefix', type=str, default='stellar')
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--night0', type=int, default=30)
    parser.add_argument('--night1', type=int, default=61)
    parser.add_argument('--write_every', type=int, default=5000000)
    parser.add_argument('--chunk_size', type=int ,default=10000)
    parser.add_argument('--dmag_cutoff', type=float, default=0.001)
    parser.add_argument('--already_done', type=str, default=None)

    args = parser.parse_args()

    if args.out_dir is None:
        raise RuntimeError('must specify out_dir')
    if args.log_file is None:
        raise RuntimeError('must specify log file')
    if os.path.exists(args.log_file):
        raise RuntimeError('%s already exists' % args.log_file)

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    already_done = set()
    if args.already_done is not None:
        with open(args.already_done, 'r') as in_file:
            for line in in_file:
                if 'INDEXING' in line:
                    try:
                        params = line.strip().split()
                        htmid = int(params[1])
                        already_done.add(htmid)
                    except:
                        pass

    opsim_db = os.path.join('/Users', 'danielsf', 'physics', 'lsst_150412',
                            'Development', 'garage', 'OpSimData',
                            'minion_1016_sqlite.db')

    if not os.path.exists(opsim_db):
        opsim_db = os.path.join('/local', 'lsst', 'danielsf', 'OpSimData',
                                'minion_1016_sqlite.db')

    obs_gen = ObservationMetaDataGenerator(opsim_db, driver='sqlite')

    obs_list = obs_gen.getObservationMetaData(night=(args.night0,args.night1))

    del obs_gen
    sims_clean_up()
    gc.collect()

    alert_gen = AlertDataGenerator()
    alert_gen.subdivide_obs(obs_list, htmid_level=6)

    n_tot_obs=0
    for htmid in alert_gen.htmid_list:
        n_tot_obs += alert_gen.n_obs(htmid)

    with open(args.log_file, 'a') as out_file:
        for htmid in alert_gen.htmid_list:
            out_file.write('htmid %d n_obs %d\n' % (htmid, alert_gen.n_obs(htmid)))
        out_file.write('n_htmid %d n_obs(total) %d\n' % (len(alert_gen.htmid_list), n_tot_obs))

    htm_population = {}
    with open('htm_population_lookup.txt', 'r') as in_file:
        for line in in_file:
            p = line.strip().split()
            htm_population[int(p[0])] = int(p[1])

    htmid_list = []
    n_htmid_list = []
    for i_p in range(args.n_proc):
        htmid_list.append([])
        n_htmid_list.append(0)

    for htmid in alert_gen.htmid_list:
        if htmid in already_done:
            continue
        n_obs = alert_gen.n_obs(htmid)*htm_population[htmid]
        n_min = -1
        i_min = -1
        for i_htmid, n_htmid in enumerate(n_htmid_list):
            if i_min<0 or n_htmid<n_min:
                n_min = n_htmid
                i_min = i_htmid
        htmid_list[i_min].append(htmid)
        n_htmid_list[i_min] += n_obs

    t_start = time.time()
    print('htmid_list %s' % str(htmid_list))
    lock = mproc.Lock()
    p_list = []
    for i_p in range(len(htmid_list)):
        p = mproc.Process(target=query_htmid,
                          args = (alert_gen, htmid_list[i_p],
                                  args.out_dir, args.out_prefix,
                                  args.log_file, lock,
                                  args.write_every,
                                  args.chunk_size, args.dmag_cutoff))

        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    with open(args.log_file, 'a') as out_file:
        elapsed = (time.time()-t_start)/3600.0
        out_file.write('all done -- %.2e hours\n' % elapsed)
        out_file.write('night0 = %d; night1 = %d\n' % (args.night0, args.night1))
