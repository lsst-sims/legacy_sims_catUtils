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
    """
    This class reads in the sqlite files created by the AlertDataGenerator
    and converts them into avro files separated by obsHistID (the unique
    integer identifying each pointing in an OpSim run).
    """

    def __init__(self):
        self._diasource_schema = None
        self._diasource_ct = {}
        self._rng = np.random.RandomState(7123)
        self._n_bit_shift = 10

    def load_schema(self, schema_dir):
        """
        Load the schema for the avro files.  Currently, these are in

        https://github.com/lsst-dm/sample-avro-alert/tree/master/schema
        """
        file_names = [os.path.join(schema_dir, 'diasource.avsc'),
                      os.path.join(schema_dir, 'diaobject.avsc'),
                      os.path.join(schema_dir, 'ssobject.avsc'),
                      os.path.join(schema_dir, 'cutout.avsc'),
                      os.path.join(schema_dir, 'alert.avsc')]

        self._alert_schema = combine_schemas(file_names)

    def _create_sources(self, obshistid, diasource_data):
        """
        Create a list of diaSources that adhere to the corresponding
        avro schema.

        Parameters
        ----------
        obshistid is an integer corresponding to the OpSim pointing
        being simulated

        diasource_data is numpy recarray containing all of the data
        for the diaSources being formatted

        Returns
        -------
        A list of dicts, each of which is ready to be written as
        an avro-formatted diaSource.
        """

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
        """
        Create a dict of diaObjects formatted according to the
        appropriate avro schema

        Parameters
        ----------
        diaobject_data is a numpy recarray containing all of the
        data needed for the diaObject

        Returns
        -------
        A dict keyed on uniqueId (the CatSim unique identifier for each
        astrophysical source).  Each value is a properly formatted
        diaObject corresponding to its key.
        """
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


    def write_alerts(self, obshistid, data_dir, prefix_list,
                     htmid_list, out_dir, out_prefix,
                     dmag_cutoff, lock=None, log_file_name=None):
        """
        Write the alerts for an obsHistId to a properly formatted avro file.

        Parameters
        ----------
        obshistid is the integer uniquely identifying the OpSim pointing
        being simulated

        data_dir is the directory containing the sqlite files created by
        the AlertDataGenerator

        prefix_list is a list of prefixes for those sqlite files.

        htmid_list is the list of htmids identifying the trixels that overlap
        this obshistid's field of view. For each htmid in htmid_list and each
        prefix in prefix_list, this method will process the files
            data_dir/prefix_htmid_sqlite.db
        searching for alerts that correspond to this obshistid

        out_dir is the directory to which the avro files should be written

        out_prefix is the prefix of the avro file names

        dmag_cutoff is the minimum delta magnitude needed to trigger an alert

        lock is an optional multiprocessing.Lock() for use when running many
        instances of this method. It prevents multiple processes from writing to
        the logfile or stdout at once.

        log_file_name is the name of an optional text file to which progress is
        written.
        """


        out_name = os.path.join(out_dir,'%s_%d.avro' % (out_prefix, obshistid))
        if os.path.exists(out_name):
            os.unlink(out_name)

        with DataFileWriter(open(out_name, "wb"),
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

        if lock is not None:
            lock.acquire()

        elapsed = (time.time()-t_start)/3600.0

        msg = 'finished obshistid %d; %d alerts in %.2e hrs' % (obshistid, alert_ct, elapsed)

        print(msg)

        if log_file_name is not None:
            with open(log_file_name, 'a') as out_file:
                out_file.write(msg)
                out_file.write('\n')

        if lock is not None:
            lock.release()
