# This script will provide classes to process the hdf5 files produced
# by the AlertDataGenerator and write them as json objects

try:
    import avro.schema
    from avro.io import DatumWriter
    from avro.datafile import DataFileWriter
except ImportError:
    pass

from lsst.sims.catalogs.db import DBObject
import os
import numpy as np
import json
import warnings
import time

__all__ = ["AvroAlertGenerator"]


################
# The methods
# combine_schemasand load_single_avsc
# are copied from Maria Patterson's
# validateAvroNestedSchema.py script in
# https://github.com:lsst-dim/sample-avro-alert
# I am copying them here since sample-avro-alert
# is not a proper eups package designed to be
# distributed and installed with the LSST stack
################

def combine_schemas(schema_files):
    """Combine multiple nested schemas into a single schema.
    """
    known_schemas = avro.schema.Names()

    for s in schema_files:
        schema = load_single_avsc(s, known_schemas)
    return schema


def load_single_avsc(file_path, names):
    """Load a single avsc file.
    """
    with open(file_path) as file_text:
        json_data = json.load(file_text)
    schema = avro.schema.SchemaFromJSONData(json_data, names)
    return schema


class AvroAlertGenerator(object):

    def __init__(self):
        self._diasource_schema = None
        self._diasource_ct = {}
        self._rng = np.random.RandomState(7123)
        self._n_bit_shift = 10

    def load_schema(self, schema_dir):
        file_names = [os.path.join(schema_dir, 'diasource.avsc'),
                      os.path.join(schema_dir, 'diaobject.avsc'),
                      os.path.join(schema_dir, 'ssobject.avsc'),
                      os.path.join(schema_dir, 'cutout.avsc'),
                      os.path.join(schema_dir, 'alert.avsc')]

        self._alert_schema = combine_schemas(file_names)

    def _create_sources(self, obshistid, diasource_data):

        bp_name_dict = {0: 'u', 1: 'g', 2: 'r', 3: 'i', 4: 'z', 5: 'y'}

        avro_diasource_list = []

        tot_flux = diasource_data['dflux'] + diasource_data['quiescent_flux']
        full_noise = tot_flux/diasource_data['tot_snr']
        quiescent_noise = diasource_data['quiescent_flux']/diasource_data['quiescent_snr']
        diff_noise = np.sqrt(full_noise**2 + quiescent_noise**2)
        diff_snr = np.abs(diasource_data['dflux']/diff_noise)

        for i_source in range(len(diasource_data)):
            diasource = diasource_data[i_source]
            if diasource['uniqueId'] not in self._diasource_ct:
                self._diasource_ct[diasource['uniqueId']] = 1

            avro_diasource = {}
            avro_diasource['diaSourceId'] = np.long((diasource['uniqueId'] << self._n_bit_shift) +
                                                 self._diasource_ct[diasource['uniqueId']])
            self._diasource_ct[diasource['uniqueId']] += 1
            avro_diasource['ccdVisitId'] = np.long((diasource['chipNum'] *10**7) + obshistid)
            avro_diasource['diaObjectId'] = np.long(diasource['uniqueId'])

            avro_diasource['midPointTai'] = diasource['TAI']
            avro_diasource['filterName'] =bp_name_dict[diasource['band']]
            avro_diasource['ra'] = diasource['ra']
            avro_diasource['decl'] = diasource['dec']
            avro_diasource['flags'] = self._rng.randint(10,1000)


            avro_diasource['x'] = diasource['xPix']
            avro_diasource['y'] = diasource['yPix']
            avro_diasource['snr'] = diff_snr[i_source]
            avro_diasource['psFlux'] = diasource['dflux']

            ra_dec_cov = {}
            ra_dec_cov['raSigma'] = self._rng.random_sample()*0.001
            ra_dec_cov['declSigma'] = self._rng.random_sample()*0.001
            ra_dec_cov['ra_decl_Cov'] = self._rng.random_sample()*0.001

            avro_diasource['ra_decl_Cov'] = ra_dec_cov

            x_y_cov = {}
            x_y_cov['xSigma'] = self._rng.random_sample()*0.001*3600.0/0.2
            x_y_cov['ySigma'] = self._rng.random_sample()*0.001*3600.0/0.2
            x_y_cov['x_y_Cov'] = self._rng.random_sample()*0.001

            avro_diasource['x_y_Cov'] = x_y_cov

            avro_diasource['totFlux'] = diasource['quiescent_flux'] + diasource['dflux']
            avro_diasource['totFluxErr'] = full_noise[i_source]  # should this be quiescent SNR?
            avro_diasource['diffFlux'] = diasource['dflux']
            avro_diasource['diffFluxErr'] = diff_noise[i_source]

            avro_diasource_list.append(avro_diasource)

        return avro_diasource_list

    def _create_objects(self, diaobject_data):
        diaobject_dict = {}
        for i_object in range(len(diaobject_data)):
            diaobject = diaobject_data[i_object]

            avro_diaobject = {}
            avro_diaobject['flags'] = np.long(self._rng.randint(10,1000))
            avro_diaobject['diaObjectId'] = np.long(diaobject['uniqueId'])
            avro_diaobject['ra'] = diaobject['ra']
            avro_diaobject['decl'] = diaobject['dec']

            ra_dec_cov = {}
            ra_dec_cov['raSigma'] = self._rng.random_sample()*0.001
            ra_dec_cov['declSigma'] = self._rng.random_sample()*0.001
            ra_dec_cov['ra_decl_Cov'] = self._rng.random_sample()*0.001

            avro_diaobject['ra_decl_Cov'] = ra_dec_cov
            avro_diaobject['radecTai'] = diaobject['TAI']

            # placeholders; we need to add this data to the sqlite files
            # generated by AlertDataGenerator
            avro_diaobject['pmRa'] = diaobject['pmRA']
            avro_diaobject['pmDecl'] = diaobject['pmDec']
            avro_diaobject['parallax'] = diaobject['parallax']
            pm_parallax_cov = {}
            for field in ('pmRaSigma', 'pmDeclSigma', 'parallaxSigma', 'pmRa_pmDecl_Cov',
                          'pmRa_parallax_Cov', 'pmDecl_parallax_Cov'):
                 pm_parallax_cov[field] = 0.0
            avro_diaobject['pm_parallax_Cov'] = pm_parallax_cov

            avro_diaobject['pmParallaxLnL'] = self._rng.random_sample()
            avro_diaobject['pmParallaxChi2'] = self._rng.random_sample()
            avro_diaobject['pmParallaxNdata'] = 0

            diaobject_dict[diaobject['uniqueId']] = avro_diaobject
        return diaobject_dict


    def write_alerts(self, obshistid, data_dir, prefix_list, htmid_list, out_dir, out_prefix):

        dmag_cutoff = 0.005

        with DataFileWriter(open(os.path.join(out_dir, "%s_%d.avro" % (out_prefix, obshistid)), "wb"),
                            DatumWriter(), self._alert_schema) as data_writer:

            diasource_query = 'SELECT alert.uniqueId, alert.xPix, alert.yPix, '
            diasource_query += 'alert.chipNum, alert.dflux, alert.snr, alert.ra, alert.dec, '
            diasource_query += 'meta.band, meta.TAI, quiescent.flux, quiescent.snr '
            diasource_query += 'FROM alert_data as alert '
            diasource_query += 'INNER JOIN metadata AS meta ON alert.obshistId=meta.obshistId '
            diasource_query += 'INNER JOIN quiescent_flux AS quiescent ON quiescent.uniqueId=alert.uniqueID '
            diasource_query += 'AND quiescent.band=meta.band '
            diasource_query += 'WHERE alert.obshistId=%d ' % obshistid
            diasource_query += 'ORDER BY alert.uniqueId'

            diasource_dtype = np.dtype([('uniqueId', int), ('xPix', float), ('yPix', float),
                                        ('chipNum', int), ('dflux', float), ('tot_snr', float),
                                        ('ra', float), ('dec', float), ('band', int), ('TAI', float),
                                        ('quiescent_flux', float), ('quiescent_snr', float)])

            diaobject_query = 'SELECT uniqueId, ra, dec, TAI, pmRA, pmDec, parallax '
            diaobject_query += 'FROM baseline_astrometry'

            diaobject_dtype = np.dtype([('uniqueId', int), ('ra', float), ('dec', float),
                                        ('TAI', float), ('pmRA', float), ('pmDec', float),
                                        ('parallax', float)])

            t_start = time.time()
            alert_ct = 0
            for htmid in htmid_list:
                for prefix in prefix_list:
                    db_name = os.path.join(data_dir, '%s_%d_sqlite.db' % (prefix, htmid))
                    print('processing %s' % db_name)
                    if not os.path.exists(db_name):
                        warnings.warn('%s does not exist' % db_name)
                        continue

                    db_obj = DBObject(db_name, driver='sqlite')

                    diaobject_data = db_obj.execute_arbitrary(diaobject_query,
                                                              dtype=diaobject_dtype)

                    diaobject_dict = self._create_objects(diaobject_data)

                    diasource_data = db_obj.execute_arbitrary(diasource_query,
                                                              dtype=diasource_dtype)

                    dmag = 2.5*np.log10(1.0+diasource_data['dflux']/diasource_data['quiescent_flux'])
                    valid_alerts = np.where(np.abs(dmag)>=dmag_cutoff)
                    print('num alerts %d num valid %d'% (len(diasource_data), len(valid_alerts[0])))
                    diasource_data = diasource_data[valid_alerts]
                    avro_diasource_list = self._create_sources(obshistid, diasource_data)

                    for i_source in range(len(avro_diasource_list)):
                        alert_ct += 1
                        unq = diasource_data[i_source]['uniqueId']
                        diaobject = diaobject_dict[unq]
                        diasource = avro_diasource_list[i_source]

                        avro_alert = {}
                        avro_alert['alertId'] = np.long((obshistid<<20) + alert_ct)
                        avro_alert['l1dbId'] = np.long(unq)
                        avro_alert['diaSource'] = diasource
                        avro_alert['diaObject'] = diaobject

                        data_writer.append(avro_alert)
                    print('ct_source %d time %.2e hrs' % (alert_ct, (time.time()-t_start)/3600.0))

        print('wrote %d' % alert_ct)
