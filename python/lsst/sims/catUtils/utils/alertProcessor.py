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
import warnings
import time

__all__ = ["AvroAlertGenerator"]


class AvroAlertGenerator(object):

    def __init__(self):
        self._diasource_schema = None
        self._diasource_ct = {}
        self._rng = np.random.RandomState(7123)

    def load_diasource_schema(self, file_name):
        with open(file_name, "rb") as input_file:
            self._diasource_schema = avro.schema.Parse(input_file.read())

    def create_diasources(self, obshistid, data_dir, prefix_list, htmid_list, n_obshistid,
                          out_file_root):

        dmag_cutoff = 0.005
        bp_name_dict = {0: 'u', 1: 'g', 2: 'r', 3: 'i', 4: 'z', 5: 'y'}

        n_bit_shift = int(np.ceil(np.log(n_obshistid)/np.log(2.0)))
        print('n bit shift %d' % (n_bit_shift))
        ct_source = 0

        with DataFileWriter(open("%s_%d.avro" % (out_file_root, obshistid), "wb"),
                            DatumWriter(), self._diasource_schema) as data_writer:

            query = 'SELECT alert.uniqueId, alert.xPix, alert.yPix, '
            query += 'alert.chipNum, alert.dflux, alert.snr, alert.ra, alert.dec, '
            query += 'meta.band, meta.TAI, quiescent.flux, quiescent.snr '
            query += 'FROM alert_data as alert '
            query += 'INNER JOIN metadata AS meta ON alert.obshistId=meta.obshistId '
            query += 'INNER JOIN quiescent_flux AS quiescent ON quiescent.uniqueId=alert.uniqueID '
            query += 'AND quiescent.band=meta.band '
            query += 'WHERE alert.obshistId=%d ' % obshistid
            query += 'ORDER BY alert.uniqueId'

            alert_dtype = np.dtype([('uniqueId', int), ('xPix', float), ('yPix', float),
                                    ('chipNum', int), ('dflux', float), ('tot_snr', float),
                                    ('ra', float), ('dec', float), ('band', int), ('TAI', float),
                                    ('quiescent_flux', float), ('quiescent_snr', float)])

            t_start = time.time()
            for htmid in htmid_list:
                for prefix in prefix_list:
                    db_name = os.path.join(data_dir, '%s_%d_sqlite.db' % (prefix, htmid))
                    print('processing %s' % db_name)
                    if not os.path.exists(db_name):
                        warnings.warn('%s does not exist' % db_name)
                        continue
                    db_obj = DBObject(db_name, driver='sqlite')
                    diasource_data = db_obj.execute_arbitrary(query, dtype=alert_dtype)

                    dmag = 2.5*np.log10(1.0+diasource_data['dflux']/diasource_data['quiescent_flux'])
                    valid_alerts = np.where(np.abs(dmag)>=dmag_cutoff)
                    print('num alerts %d num valid %d'% (len(diasource_data), len(valid_alerts[0])))
                    diasource_data = diasource_data[valid_alerts]

                    tot_flux = diasource_data['dflux'] + diasource_data['quiescent_flux']
                    full_noise = tot_flux/diasource_data['tot_snr']
                    quiescent_noise = diasource_data['quiescent_flux']/diasource_data['quiescent_snr']
                    diff_noise = np.sqrt(full_noise**2 + quiescent_noise**2)
                    diff_snr = np.abs(diasource_data['dflux']/diff_noise)

                    for i_source in range(len(diasource_data)):
                        diasource = diasource_data[i_source]
                        if diasource['uniqueId'] not in self._diasource_ct:
                            self._diasource_ct[diasource['uniqueId']] = 1

                        avro_source = {}
                        avro_source['diaSourceId'] = np.long((diasource['uniqueId'] << n_bit_shift) +
                                                             self._diasource_ct[diasource['uniqueId']])
                        self._diasource_ct[diasource['uniqueId']] += 1
                        avro_source['ccdVisitId'] = np.long((diasource['chipNum'] *10**7) + obshistid)
                        avro_source['diaObjectId'] = np.long(diasource['uniqueId'])

                        avro_source['midPointTai'] = diasource['TAI']
                        avro_source['filterName'] =bp_name_dict[diasource['band']]
                        avro_source['ra'] = diasource['ra']
                        avro_source['decl'] = diasource['dec']
                        avro_source['flags'] = self._rng.randint(10,1000)


                        avro_source['x'] = diasource['xPix']
                        avro_source['y'] = diasource['yPix']
                        avro_source['snr'] = diff_snr[i_source]
                        avro_source['psFlux'] = diasource['dflux']

                        ra_dec_cov = {}
                        ra_dec_cov['raSigma'] = self._rng.random_sample()*0.001
                        ra_dec_cov['declSigma'] = self._rng.random_sample()*0.001
                        ra_dec_cov['ra_decl_Cov'] = self._rng.random_sample()*0.001

                        avro_source['ra_decl_Cov'] = ra_dec_cov

                        x_y_cov = {}
                        x_y_cov['xSigma'] = self._rng.random_sample()*0.001*3600.0/0.2
                        x_y_cov['ySigma'] = self._rng.random_sample()*0.001*3600.0/0.2
                        x_y_cov['x_y_Cov'] = self._rng.random_sample()*0.001

                        avro_source['x_y_Cov'] = x_y_cov

                        avro_source['totFlux'] = diasource['quiescent_flux'] + diasource['dflux']
                        avro_source['totFluxErr'] = full_noise[i_source]  # should this be quiescent SNR?
                        avro_source['diffFlux'] = diasource['dflux']
                        avro_source['diffFluxErr'] = diff_noise[i_source]
                        avro_source['fpBkgd'] = self._rng.random_sample()
                        avro_source['fpBkgdErr'] = self._rng.random_sample()

                        data_writer.append(avro_source)
                        ct_source += 1
                    print('ct_source %d time %.2e hrs' % (ct_source, (time.time()-t_start)/3600.0))

        print('wrote %d' % ct_source)
