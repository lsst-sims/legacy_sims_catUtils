import numpy as np
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AlertDataGenerator

import os
from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up

from lsst.sims.catUtils.utils import AvroAlertGenerator
import time

if __name__ == "__main__":

    data_dir = 'sql_data'
    out_dir = 'avro_out_dir'
    sql_prefix_list = ['stellar', 'agn']
    physics_dir = os.path.join('/Users', 'danielsf', 'physics')
    lsst_dir = os.path.join(physics_dir, 'lsst_171025')
    schema_dir = os.path.join(lsst_dir, 'Development', 'sample-avro-alert',
                              'schema')
    diasource_schema = os.path.join(schema_dir, 'diasource.avsc')

    opsim_db = os.path.join('/Users', 'danielsf', 'physics', 'lsst_150412',
                            'Development', 'garage', 'OpSimData',
                            'minion_1016_sqlite.db')

    if not os.path.exists(opsim_db):
        opsim_db = os.path.join('/local', 'lsst', 'danielsf', 'OpSimData',
                                'minion_1016_sqlite.db')

    obs_gen = ObservationMetaDataGenerator(opsim_db, driver='sqlite')

    obs_list = obs_gen.getObservationMetaData(night=15)

    del obs_gen
    sims_clean_up()

    print('%d obs' % len(obs_list))

    alert_gen = AlertDataGenerator()
    alert_gen.subdivide_obs(obs_list, htmid_level=6)

    obshistid_list = []
    tai_list = []
    for obs in obs_list:
       obshistid_list.append(obs.OpsimMetaData['obsHistID'])
       tai_list.append(obs.mjd.TAI)
    obshistid_list = np.array(obshistid_list)
    tai_list = np.array(tai_list)
    sorted_dex = np.argsort(obshistid_list)
    sorted_dex_2 = np.argsort(tai_list)
    np.testing.assert_array_equal(sorted_dex, sorted_dex_2)
    tai_list = tai_list[sorted_dex]
    obshistid_list = obshistid_list[sorted_dex]

    obshistid_to_htmid = {}
    for htmid in alert_gen.htmid_list:
        for obs in alert_gen.obs_from_htmid(htmid):
            obshistid = obs.OpsimMetaData['obsHistID']
            if obshistid not in obshistid_to_htmid:
                obshistid_to_htmid[obshistid] = []
            obshistid_to_htmid[obshistid].append(htmid)

    avro_gen = AvroAlertGenerator()
    avro_gen.load_schema(schema_dir)

    out_file_root = os.path.join('avro_out_dir', 'test_avro')
    avro_gen.write_alerts(obshistid_list[0], data_dir,
                          sql_prefix_list,
                          obshistid_to_htmid[obshistid_list[0]],
                          out_file_root)

