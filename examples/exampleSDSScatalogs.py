import os
import numpy

from lsst.utils import getPackageDir
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.mixins import PhotometryStars, PhotometryGalaxies
from lsst.sims.catUtils.mixins import EBVmixin
from lsst.sims.catalogs.measures.instance import compound, InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.utils import ObservationMetaData

from lsst.sims.catUtils.baseCatalogModels import *

class sdssGalaxies(InstanceCatalog,EBVmixin,PhotometryGalaxies):

    catalog_type = 'sdssGalaxies'
    column_outputs = ['galid','sdss_u','sdss_g','sdss_r','sdss_i','sdss_z',
                      'sdss_bulge_u','sdss_bulge_g','sdss_bulge_r','sdss_bulge_i','sdss_bulge_z',
                      'sdss_disk_u','sdss_disk_g','sdss_disk_r','sdss_disk_i','sdss_disk_z',
                      'sdss_agn_u','sdss_agn_g','sdss_agn_r','sdss_agn_i','sdss_agn_z']

    @compound('sdss_bulge_u', 'sdss_bulge_g', 'sdss_bulge_r', 'sdss_bulge_i', 'sdss_bulge_z')
    def get_sdss_bulge_mags(self):
        """
        An example getter for SDSS bulge magnitudes
        """

        # load a BandpassDict of SDSS bandpasses, if not done already
        if not hasattr(self, 'sdssBandpassDict'):
            bandpassNames = ['u','g','r','i','z']
            bandpassDir = os.path.join(getPackageDir('throughputs'),'sdss')
            bandpassRoot = 'sdss_'

            self.sdssBandpassDict = BandpassDict.loadTotalBandpassesFromFiles(bandpassNames,
                                                                 bandpassRoot = bandpassRoot,
                                                                 bandpassDir = bandpassDir)

        # actually calculate the magnitudes
        return self._magnitudeGetter('bulge', self.sdssBandpassDict,
                                     self.get_sdss_bulge_mags._colnames)


    @compound('sdss_disk_u', 'sdss_disk_g', 'sdss_disk_r', 'sdss_disk_i', 'sdss_disk_z')
    def get_sdss_disk_mags(self):
        """
        An example getter for SDSS disk magnitudes
        """

        # load a BandpassDict of SDSS bandpasses, if not done already
        if not hasattr(self, 'sdssBandpassDict'):
            bandpassNames = ['u','g','r','i','z']
            bandpassDir = os.path.join(getPackageDir('throughputs'),'sdss')
            bandpassRoot = 'sdss_'

            self.sdssBandpassDict = BandpassDict.loadTotalBandpassesFromFiles(bandpassNames,
                                                                 bandpassRoot = bandpassRoot,
                                                                 bandpassDir = bandpassDir)

        # actually calculate the magnitudes
        return self._magnitudeGetter('disk', self.sdssBandpassDict,
                                     self.get_sdss_disk_mags._colnames)


    @compound('sdss_agn_u', 'sdss_agn_g', 'sdss_agn_r', 'sdss_agn_i', 'sdss_agn_z')
    def get_sdss_agn_mags(self):
        """
        An example getter for SDSS AGN magnitudes
        """

        # load a BandpassDict of SDSS bandpasses, if not done already
        if not hasattr(self, 'sdssBandpassDict'):
            bandpassNames = ['u','g','r','i','z']
            bandpassDir = os.path.join(getPackageDir('throughputs'),'sdss')
            bandpassRoot = 'sdss_'

            self.sdssBandpassDict = BandpassDict.loadTotalBandpassesFromFiles(bandpassNames,
                                                                 bandpassRoot = bandpassRoot,
                                                                 bandpassDir = bandpassDir)

        # actually calculate the magnitudes
        return self._magnitudeGetter('agn', self.sdssBandpassDict,
                                     self.get_sdss_agn_mags._colnames)



    @compound('sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z', 'sdss_y')
    def get_sdss_total_mags(self):
        """
        An example getter for total galaxy magnitudes in SDSS bands
        """
        idList = self.column_by_name('uniqueId')
        numObj = len(idList)
        output = []

        # Go through each of the calculated columns.  Find the bulge, disk, and agn
        # magnitudes corresponding to each bandpass.  Sum them using the
        # sum_magnitudes method
        for columnName in self.get_sdss_total_mags._colnames:
            if columnName not in self._actually_calculated_columns:
                sub_list = [numpy.NaN]*numObj
            else:
                bandpass = columnName[-1]
                bulge = self.column_by_name('sdss_bulge_%s' % bandpass)
                disk = self.column_by_name('sdss_disk_%s' % bandpass)
                agn = self.column_by_name('sdss_agn_%s' % bandpass)
                sub_list = self.sum_magnitudes(bulge=bulge, disk=disk, agn=agn)

            output.append(sub_list)
        return numpy.array(output)



class sdssStars(InstanceCatalog,PhotometryStars):

    catalog_type = 'sdssStars'
    column_outputs = ['id','sdss_u','sdss_g','sdss_r','sdss_i','sdss_z']

    @compound('sdss_u','sdss_g','sdss_r','sdss_i','sdss_z')
    def get_sdss_magnitudes(self):
        """
        An example getter for stellar magnitudes in SDSS bands
        """

        # Load a BandpassDict of SDSS bandpasses, if not done already
        if not hasattr(self, 'sdssBandpassDict'):
            bandpassNames = ['u','g','r','i','z']
            bandpassDir = os.path.join(getPackageDir('throughputs'),'sdss')
            bandpassRoot = 'sdss_'

            self.sdssBandpassDict = BandpassDict.loadTotalBandpassesFromFiles(bandpassNames,
                                                                 bandpassRoot = bandpassRoot,
                                                                 bandpassDir = bandpassDir)

        # Actually calculate the magnitudes
        return self._magnitudeGetter(self.sdssBandpassDict, self.get_sdss_magnitudes._colnames)


if __name__ == "__main__":

    obs_metadata_pointed = ObservationMetaData(mjd=2013.23, boundType='circle',
                                               pointingRA=200.0, pointingDec=-30.0, boundLength=1.0)

    obs_metadata_pointed.metadata = {}
    obs_metadata_pointed.metadata['Opsim_filter'] = 'i'
    dbObj = CatalogDBObject.from_objid('rrlystars')
    sdssStars = sdssStars(dbObj, obs_metadata = obs_metadata_pointed)
    sdssStars.write_catalog("example_sdss_stars.txt")

    obs_metadata_pointed = ObservationMetaData(mjd=50000.0, boundType='circle',
                             pointingRA=0.0, pointingDec=0.0, boundLength=0.01)

    obs_metadata_pointed.metadata = {}
    obs_metadata_pointed.metadata['Opsim_filter'] = 'i'
    dbObj = CatalogDBObject.from_objid('galaxyBase')
    sdssGalaxies = sdssGalaxies(dbObj, obs_metadata = obs_metadata_pointed)
    sdssGalaxies.write_catalog("example_sdss_galaxies.txt")
