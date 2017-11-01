# This script will provide classes to process the hdf5 files produced
# by the AlertDataGenerator and write them as json objects

try:
    import avro.schema
    from avro.io import DatumWriter
    from avro.datafile import DataFileWriter
except ImportError:
    pass

import h5py
import os
import numpy as np

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator


__all__ = ["AlertProcessor"]


class AlertProcessor(object):

    def __init__(self):
        self._open_file_handle = None
        self._open_file_name = None
        self._diasource_schema = None
        self._rng = np.random.RandomState(7123)

    def __del__(self):
        self._close_file()

    def load_diasource_schema(self, file_name):
        with open(file_name, "rb") as input_file:
            self._diasource_schema = avro.schema.parse(input_file.read())

    def _close_file(self):
        if self._open_file_handle is None:
            if self._open_file_name is not None:
                raise RuntimeError("open file is None; file name is %s" %
                                   self._open_file_name)

        if self._open_file_name is not None:
            if self._open_file_handle is None:
                raise RuntimeError("open file is %s; file name is None" %
                                   self._open_file_handle)

        if self._open_file_handle is not None:
            self._open_file_handle.close()
            self._open_file_handle = None
            self._open_file_name = None


    def _open_file(self, file_name):

        if self._open_file_handle is None:
            if self._open_file_name is not None:
                raise RuntimeError("open file is None; file name is %s" %
                                   self._open_file_name)

        if self._open_file_name is not None:
            if self._open_file_handle is None:
                raise RuntimeError("open file is %s; file name is None" %
                                   self._open_file_handle)

        if file_name == self._open_file_name:
            return

        self._close_file()
        self._open_file_handle = h5py.File(file_name, 'r')
        self._open_file_name = file_name

    def _process_obs(self, obshistid, tai, bandpass, out_file_root):
        hdf_file_name = self._obshistid_to_file_name[obshistid]
        self._open_file(hdf_file_name)

        bp_name_dict = {0: 'u', 1: 'g', 2: 'r', 3: 'i', 4: 'z', 5: 'y'}

        with DataFileWriter(open("%s_%d.avro" % (out_file_root, obshistid), "wb"),
                            DatumWriter(), self._diasource_schema) as data_writer:

            for batch_ct in self._obshistid_ct_list[obshistid]:
                field_root = '%d_%d' % (obshistid, batch_ct)
                id_arr = self._open_file_handle['%s_uniqueId' % field_root].value
                ra_arr = self._open_file_handle['%s_raICRS' % field_root].value
                dec_arr = self._open_file_handle['%s_decICRS' % field_root].value
                flux_arr = self._open_file_handle['%s_flux' % field_root].value
                dflux_arr = self._open_file_handle['%s_dflux' % field_root].value
                snr_arr = self._open_file_handle['%s_SNR' % field_root].value
                chipnum_arr = self._open_file_handle['%s_chipNum' % field_root].value
                xpix_arr = self._open_file_handle['%s_xPix' % field_root].value
                ypix_arr = self._open_file_handle['%s_yPix' % field_root].value

                for unqId, ra, dec, flux, dflux, snr, chipnum, xpix, ypix in \
                zip(id_arr, ra_arr, dec_arr, flux_arr, dflux_arr, snr_arr, chipnum_arr,
                    xpix_arr, ypix_arr):

                    source = {}
                    source['diaSourceId'] = unqId << 22 + obshistid
                    source['ccdVisitId'] = chipnum*10**7 + obshistid
                    source['diaObjectId'] = unqId
                    source['midPointTai'] = tai
                    source['filtername'] = bp_name_dict[bandpass]
                    source['ra'] = ra
                    source['decl'] = dec

                    ra_dec_cov = {}
                    ra_dec_cov['raSigma'] = self._rng.random_sample()*0.001
                    ra_dec_cov['declSigma'] = self._rng.random_sample()*0.001
                    ra_dec_cov['ra_decl_Cov'] = self._rng.random_sample()*0.001

                    source['ra_decl_Cov'] = ra_dec_cov

                    source['x'] = xpix
                    source['y'] = ypix

                    x_y_cov = {}
                    x_y_cov['xSigma'] = self._rng.random_sample()*0.001*3600.0/0.2
                    x_y_cov['ySigma'] = self._rng.random_sample()*0.001*3600.0/0.2
                    x_y_cov['x_y_Cov'] = self._rng.random_sample()*0.001

                    source['x_y_Cov'] = x_y_cov

                    apFlux = {}
                    apFlux['apMeanSb01'] = dflux * (1.0+self._rng.random_sample()-0.5)
                    apFlux['apMeanSb02'] = dflux * (1.0+self._rng.random_sample()-0.5)
                    apFlux['apMeanSb03'] = dflux * (1.0+self._rng.random_sample()-0.5)
                    apFlux['apMeanSb04'] = dflux * (1.0+self._rng.random_sample()-0.5)
                    apFlux['apMeanSb05'] = dflux * (1.0+self._rng.random_sample()-0.5)
                    apFlux['apMeanSb06'] = dflux * (1.0+self._rng.random_sample()-0.5)
                    apFlux['apMeanSb07'] = dflux * (1.0+self._rng.random_sample()-0.5)
                    apFlux['apMeanSb08'] = dflux * (1.0+self._rng.random_sample()-0.5)
                    apFlux['apMeanSb09'] = dflux * (1.0+self._rng.random_sample()-0.5)
                    apFlux['apMeanSb10'] = dflux * (1.0+self._rng.random_sample()-0.5)

                    source['apFlux'] = apFlux

                    apFluxErr = {}
                    apFluxErr['apMeanSb01Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr
                    apFluxErr['apMeanSb02Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr
                    apFluxErr['apMeanSb03Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr
                    apFluxErr['apMeanSb04Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr
                    apFluxErr['apMeanSb05Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr
                    apFluxErr['apMeanSb06Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr
                    apFluxErr['apMeanSb07Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr
                    apFluxErr['apMeanSb08Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr
                    apFluxErr['apMeanSb09Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr
                    apFluxErr['apMeanSb10Err'] = dflux*(1.0+self._rng.random_sample()-0.5)/snr

                    source['apFluxErr'] = apFluxErr

                    source['snr'] = snr
                    source['psFlux'] = dflux*(1.0 + self._rng.random_sample()*0.2)
                    source['psRa'] = ra + (self._rng.random_sample()-0.5)*0.001
                    source['psDecl'] = dec + (self._rng.random_sample()-0.5)*0.001

                    ps_cov = {}
                    ps_cov['psFluxSigma'] = dflux/snr
                    ps_cov['psRaSigma'] = self._rng.random_sample()*0.001
                    ps_cov['psDeclSigma'] = self._rng.random_sample()*0.001
                    ps_cov['psFlux_psRa_Cov'] = self._rng.random_sample()*0.001
                    ps_cov['psFlux_psDecl_Cov'] = self._rng.random_sample()*0.001
                    ps_cov['psRa_psDecl_Cov'] = self._rng.random_sample()*0.001

                    source['ps_Cov'] = ps_cov

                    chi2 = self._rng.random_sample()*10.0
                    lnl = -0.5*chi2
                    source['psLnL'] = lnl
                    source['psChi2'] = chi2
                    source['psNdata'] = self._rng.randint(10, 100)

                    source['trailFlux'] = dflux + self._rng.random_sample()*0.1
                    source['trailRa'] = ra + self._rng.random_sample()*0.1
                    source['trailDecl'] = dec + self._rng.random_sample()*0.1
                    source['trailLength']  = self._rng.random_sample*0.2
                    source['trailAngle'] = self._rng.random_sample()*360.0

                    trail_cov = {}
                    trail_cov['trailFluxSigma'] = self._rng.random_sample()*0.01
                    trail_cov['trailRaSigma'] = self._rng.random_sample()*0.01
                    trail_cov['trailDeclSigma'] = self._rng.random_sample()*0.01
                    trail_cov['trailLengthSigma'] = self._rng.random_sample()*0.01
                    trail_cov['trailAngleSigma'] = self._rng.random_sample()*0.1
                    trail_cov['trailFlux_trailRa_Cov'] = self._rng.random_sample()*0.01
                    trail_cov['trailFlux_trailDecl_Cov'] = self._rng.random_sample()*0.01
                    trail_cov['trailFlux_trailLength_Cov'] = self._rng.random_sample()*0.01
                    trail_cov['trailFlux_trailAngle_Cov'] = self._rng.random_sample()*0.01
                    trail_cov['trailRa_trailDecl_Cov'] = self._rng.random_sample()*0.01
                    trail_cov['trailRa_trailLength_Cov'] = self._rng.random_sample()*0.01
                    trail_cov['trailRa_trailAngle_Cov'] = self._rng.random_sample()*0.01
                    trail_cov['trailDecl_trailLength_Cov'] = self._rng.random_sample()*0.01
                    trail_cov['trailDecl_trailAngle_Cov'] = self._rng.random_sample()*0.01
                    trail_cov['trailLength_trailAngle_Cov'] = self._rng.random_sample()*0.01

                    source['trail_Cov'] = trail_cov

                    chi2 = self._rng.random_sample()*100.0
                    lnl = -0.5*chi2
                    source['trailLnL'] = lnl
                    source['trailChi2'] = chi2
                    source['trailNdata'] = self._rng.randint(10, 100)

                    source['dipMeanFlux'] = dflux * (1.0+self._rng.random_sample()*0.2)
                    source['dipFluxDiff'] = dflux * (1.0+self._rng.random_sample()*0.2)
                    source['dipRa'] = ra + self._rng.random_sample()*0.001
                    source['dipDecl'] = dec + self._rng.random_sample()*0.001
                    source['dipLength'] = self._rng.random_sample()*0.1
                    source['dipAngle'] = self._rng.random_sample()*360.0

                    dip_cov = {}
                    dip_cov['dipMeanFluxSigma'] = dflux*(1.0+self._rng.random_sample()*0.2)/snr
                    dip_cov['dipFluxDiffSigma'] = dflux*(1.0+self._rng.random_sample()*0.2)/snr
                    dip_cov['dipRaSigma'] = self._rng.random_sample()*0.001
                    dip_cov['dipDeclSigma'] = self._rng.random_sample()*0.001
                    dip_cov['dipLengthSigma'] = self._rng.random_sample()*0.1
                    dip_cov['dipAngleSigma'] = self._rng.random_sample()*10.0
                    dip_cov['dipMeanFlux_dipFluxDiff_Cov'] = self._rng.random_sample()
                    dip_cov['dipMeanFlux_dipRa_Cov'] = self._rng.random_sample()
                    dip_cov['dipMeanFlux_dipDecl_Cov'] = self._rng.random_sample()
                    dip_cov['dipMeanFlux_dipLength_Cov'] = self._rng.random_sample()
                    dip_cov['dipMeanFlux_dipAngle_Cov'] = self._rng.random_sample()
                    dip_cov['dipFluxDiff_dipRa_Cov'] = self._rng.random_sample()
                    dip_cov['dipFluxDiff_dipDecl_Cov'] = self._rng.random_sample()
                    dip_cov['dipFluxDiff_dipLength_Cov'] = self._rng.random_sample()
                    dip_cov['dipFluxDiff_dipAngle_Cov'] = self._rng.random_sample()
                    dip_cov['dipRa_dipDecl_Cov'] = self._rng.random_sample()
                    dip_cov['dipRa_dipLength_Cov'] = self._rng.random_sample()
                    dip_cov['dipRa_dipAngle_Cov'] = self._rng.random_sample()
                    dip_cov['dipDecl_dipLength_Cov'] = self._rng.random_sample()
                    dip_cov['dipDecl_dipAngle_Cov'] = self._rng.random_sample()
                    dip_cov['dipLength_dipAngle_Cov'] = self._rng.random_sample()

                    source['dip_Cov'] = dip_cov

                    chi2 = self._rng.random_sample()*100.0
                    lnl = -0.5*chi2
                    source['dipLnL'] = lnl
                    source['dipChi2'] = chi2
                    source['dpNdata'] = self._rng.randint(10, 100)
                    source['totFlux'] = flux
                    source['totFluxErr'] = flux*self._rng.random_sample()  # should this be quiescent SNR?
                    source['diffFlux'] = dflux
                    source['diffFluxErr'] = dflux/snr
                    source['fpBkgd'] = self._rng.random_sample()
                    source['fpBkgdErr'] = self._rng.random_sample()
                    source['ixx'] = self._rng.random_sample()
                    source['iyy'] = self._rng.random_sample()
                    source['ixy'] = self._rng.random_sample()

                    i_cov = {}
                    i_cov['ixxSigma'] = self._rng.random_sample()
                    i_cov['iyySigma'] = self._rng.random_sample()
                    i_cov['ixySigam'] = self._rng.random_sample()
                    i_cov['ixx_iyy_Cov'] = self._rng.random_sample()
                    i_cov['ixx_ixy_Cov'] = self._rng.random_sample()
                    i_cov['iyy_ixy_Cov'] = self._rng.random_sample()

                    source['i_cov'] = i_cov

                    source['ixxPSF'] = self._rng.random_sample()
                    source['iyyPSF'] = self._rng.random_sample()
                    source['ixyPSF'] = self._rng.random_sample()
                    source['extendedness'] = self._rng.random_sample()
                    source['spuriousness'] = self._rng.random_sample()
                    source['flags'] = self._rng.randint(0,1000)

                    data_writer.append(source)




    def process(self, opsimdb, hdf5_dir, opsim_driver='sqlite'):

        if self._diasource_schema is None:
            raise RuntimeError("Need to specify diasource_schema")

        obs_gen = ObservationMetaDataGenerator(opsimdb, driver=opsim_driver)
        hdf5_list = []
        for file_name in os.listdir(hdf5_dir):
            if file_name.endswith('hdf5'):
                hdf5_list.append(os.path.join(hdf5_dir, file_name))

        # First, scan the hdf5 files and figure out what observations
        # are available.  Sort them into chronological order.
        obshistid_list = []
        tai_list = []
        bandpass_list = []
        self._obshistid_to_file_name = {}
        self_obhistid_ct_list = {}
        for file_name in hdf5_list:
            with h5py.File(file_name, 'r') as input_file:
                for obshistid, tai, bp in zip(input_file['obshistID'].value,
                                              input_file['TAI'].value,
                                              input_file['bandpass'].value):

                    obshistid_list.append(obhistid)
                    tai_list.append(tai)
                    bandpass_list.append(bp)
                    self._obshistid_to_file_name[obshistid] = file_name
                    self._obshistid_ct_list[obshistid] = []

                # figure out how many chunks each obsHistID is stored in
                for arr_key in input_file.keys():
                    params = arr_key.split('_')
                    try:
                        obshistid = int(params[0])
                    except ValueError:
                        continue
                    self._obshistid_ct_list[obshistid].append(int(params[1]))

        obshistid_list = np.array(obshistid_list)
        tai_list = np.array(tai_list)
        bandpass_list = np.array(bandpass_list)

        sorted_dex = np.argsort(tai_list)
        tai_list = tai_list[sorted_dex]
        obshistid_list = obshistid_list[sorted_dex]
        bandpass_list = bandpass_list[sorted_dex]
