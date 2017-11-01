# This script will provide classes to process the hdf5 files produced
# by the AlertDataGenerator and write them as json objects

import h5py
import os
import numpy as np

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator

__all__ = ["AlertProcessor"]


class AlertProcessor(object):

    def __init__(self):
        self._open_file_handle = None
        self._open_file_name = None

    def __del__(self)
        self._close_file()

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

    def _process_obs(self, obshistid, tai, bandpass, out_file_handle):
        hdf_file_name = self._obshistid_to_file_name[obshistid]
        self._open_file(hdf_file_name)



    def process(self, opsimdb, hdf5_dir, opsim_driver='sqlite'):

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
