from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from EclipsingBinaryModule import EclipsingBinaryAlertDataGenerator

import os
from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up

from EclipsingBinaryModule import OBAFGKObj
from lsst.sims.catUtils.utils import AlertStellarVariabilityCatalog

import multiprocessing as mproc

import time
import gc
import argparse


def query_htmid(alert_gen, htmid_list, out_dir, out_prefix,
                log_file_name, lock,
                write_every, chunk_size, dmag_cutoff):

    try:
        db = OBAFGKObj(database='LSSTCATSIM',
                       host='fatboy.phys.washington.edu',
                       port=1433, driver='mssql+pymssql',
                       cache_connection=False)
    except RuntimeError:
        db = OBAFGKObj(cache_connection=False)

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
                                                 log_file_name=log_file_name,
                                                 constraint='imag<24.0',
                                                 limit=100000)

        lock.acquire()
        with open(log_file_name, 'a') as out_file:
            elapsed = (time.time()-t_start)/3600.0
            out_file.write('htmid %d nobs %d n_rows %d time %.2e hrs; per_row %.2e -- n_htmid %d of %d\n' %
                           (htmid, alert_gen.n_obs(htmid), n_rows, elapsed,
                            elapsed/n_rows, i_htmid, len(htmid_list)))

        lock.release()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--n_proc', type=int, default=1,
                        help='Number of independent process to run')
    parser.add_argument('--out_dir', type=str, default=None,
                        help='Directory in which to write output')
    parser.add_argument('--out_prefix', type=str, default='stellar',
                        help='Prefix for output sqlite files')
    parser.add_argument('--log_file', type=str, default=None,
                        help='Name of file to write progress to')
    parser.add_argument('--night0', type=int, default=30,
                        help='First night of the survey to '
                        'simulate (default=30)')
    parser.add_argument('--night1', type=int, default=61,
                        help='Last night of the survey to '
                        'simulate (default=61)')
    parser.add_argument('--write_every', type=int, default=5000000,
                        help='Write output to sqlite files every write_every rows '
                        '(default = 5e6)')
    parser.add_argument('--chunk_size', type=int ,default=10000,
                        help='Number of rows to bring down from the remote '
                        'database at a time (default=10,000)')
    parser.add_argument('--dmag_cutoff', type=float, default=0.001,
                        help='Minimum delta magnitude to trigger an alert '
                        '(default = 0.001)')
    parser.add_argument('--opsim_db', type=str,
                        default=os.path.join('/local', 'lsst', 'danielsf',
                                             'OpSimData', 'minion_1016_sqlite.db'),
                        help='Path to OpSim database used for survey cadence')

    args = parser.parse_args()

    if args.out_dir is None:
        raise RuntimeError('must specify out_dir')
    if args.log_file is None:
        raise RuntimeError('must specify log file')
    if os.path.exists(args.log_file):
        raise RuntimeError('%s already exists' % args.log_file)

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    # get the list of ObservationMetaData to simulate
    obs_gen = ObservationMetaDataGenerator(args.opsim_db, driver='sqlite')
    g_obs_list = obs_gen.getObservationMetaData(night=(args.night0,args.night1),
                                                telescopeFilter='g')

    i_obs_list = obs_gen.getObservationMetaData(night=(args.night0,args.night1),
                                                telescopeFilter='i')

    obs_list = []
    for obs in g_obs_list:
        obs_list.append(obs)
    for obs in i_obs_list:
        obs_list.append(obs)

    del obs_gen
    sims_clean_up()
    gc.collect()

    # get the list of trixel htmids to simulate
    alert_gen = EclipsingBinaryAlertDataGenerator()
    alert_gen.subdivide_obs(obs_list, htmid_level=6)

    n_tot_obs=0
    for htmid in alert_gen.htmid_list:
        n_tot_obs += alert_gen.n_obs(htmid)

    with open(args.log_file, 'a') as out_file:
        for htmid in alert_gen.htmid_list:
            out_file.write('htmid %d n_obs %d\n' % (htmid, alert_gen.n_obs(htmid)))
        out_file.write('n_htmid %d n_obs(total) %d\n' % (len(alert_gen.htmid_list), n_tot_obs))

    # subdivide the trixels into n_proc lists, trying to make sure
    # each process has an equal number of observations to simulate
    htmid_list = []
    n_htmid_list = []
    for i_p in range(args.n_proc):
        htmid_list.append([])
        n_htmid_list.append(0)

    for htmid in alert_gen.htmid_list:
        n_obs = alert_gen.n_obs(htmid)
        n_min = -1
        i_min = -1
        for i_htmid, n_htmid in enumerate(n_htmid_list):
            if i_min<0 or n_htmid<n_min:
                n_min = n_htmid
                i_min = i_htmid
        htmid_list[i_min].append(htmid)
        n_htmid_list[i_min] += n_obs

    # start the independent processes
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
