from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AlertDataGenerator

import os
from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up

from lsst.sims.catUtils.utils import StellarAlertDBObjMixin
from lsst.sims.catUtils.baseCatalogModels import StarBase
from lsst.sims.catUtils.utils import AlertStellarVariabilityCatalog


import multiprocessing as mproc

import time
import gc
import argparse

class MLTObj(StarBase):
    objid = 'obafgkstars'
    tableid = 'StarMLTForceseek'
    objectTypeId = 14
    doRunTest = True
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 40)]


class MLTAlertDBObj(StellarAlertDBObjMixin, MLTObj):
    pass


def query_htmid(alert_gen, htmid_list, out_dir, out_prefix,
                log_file_name, lock, stdout_lock, ct_lock,
                ct_dict, write_every):

    db = MLTAlertDBObj(database='LSSTCATSIM',
                      host='fatboy.phys.washington.edu',
                      port=1433, driver='mssql+pymssql',
                      cache_connection=False)

    for i_htmid, htmid in enumerate(htmid_list):
        t_start = time.time()
        print('passing in %s' % str(htmid))
        n_rows = alert_gen.alert_data_from_htmid(htmid, db, chunk_size=5000,
                                                 write_every=write_every,
                                                 output_dir=out_dir,
                                                 output_prefix=out_prefix,
                                                 dmag_cutoff=0.001,
                                                 photometry_class=AlertStellarVariabilityCatalog,
                                                 lock=lock,
                                                 stdout_lock=stdout_lock,
                                                 ct_lock=ct_lock,
                                                 ct_dict=ct_dict,
                                                 log_file_name=log_file_name)

        lock.acquire()
        with open(log_file_name, 'a') as out_file:
            elapsed = (time.time()-t_start)/3600.0
            out_file.write('\nFINISHED htmid %d nobs %d n_rows %d time %.2e hrs; per_row %.2e -- n_htmid %d of %d\n' %
                           (htmid, alert_gen.n_obs(htmid), n_rows, elapsed,
                            elapsed/n_rows, i_htmid, len(htmid_list)))

        lock.release()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--n_proc', type=int, default=1)
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--out_prefix', type=str, default='stellar')
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--write_every', type=int, default=5000000)

    args = parser.parse_args()

    if args.out_dir is None:
        raise RuntimeError('must specify out_dir')
    if args.log_file is None:
        raise RuntimeError('must specify log file')
    if os.path.exists(args.log_file):
        raise RuntimeError('%s already exists' % args.log_file)

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    opsim_db = os.path.join('/Users', 'danielsf', 'physics', 'lsst_150412',
                            'Development', 'garage', 'OpSimData',
                            'minion_1016_sqlite.db')

    if not os.path.exists(opsim_db):
        opsim_db = os.path.join('/local', 'lsst', 'danielsf', 'OpSimData',
                                'minion_1016_sqlite.db')

    obs_gen = ObservationMetaDataGenerator(opsim_db, driver='sqlite')

    obs_list = obs_gen.getObservationMetaData(night=(30,61))

    del obs_gen
    sims_clean_up()
    gc.collect()

    alert_gen = AlertDataGenerator()
    alert_gen.subdivide_obs(obs_list)

    htmid_list = []
    n_htmid_list = []
    for i_p in range(args.n_proc):
        htmid_list.append([])
        n_htmid_list.append(0)

    n_tot_obs=0
    for htmid in alert_gen.htmid_list:
        n_tot_obs += alert_gen.n_obs(htmid)

    with open(args.log_file, 'a') as out_file:
        for htmid in alert_gen.htmid_list:
            out_file.write('htmid %d n_obs %d\n' % (htmid, alert_gen.n_obs(htmid)))
        out_file.write('n_htmid %d n_obs(total) %d\n' % (len(alert_gen.htmid_list), n_tot_obs))

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

    print('htmid_list %s' % str(htmid_list))
    lock = mproc.Lock()
    stdout_lock = mproc.Lock()
    mgr = mproc.Manager()
    ct_dict = mgr.dict()
    ct_lock = mproc.Lock()
    ct_dict['number_writing'] = 0
    ct_dict['allowed_to_write'] = 5
    p_list = []
    for i_p in range(len(htmid_list)):
        p = mproc.Process(target=query_htmid,
                          args = (alert_gen, htmid_list[i_p],
                                  args.out_dir, args.out_prefix,
                                  args.log_file, lock, stdout_lock,
                                  ct_lock, ct_dict, args.write_every))

        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    with open(args.log_file, 'a') as out_file:
        out_file.write('all done')
