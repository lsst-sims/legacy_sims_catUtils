import os

from lsst.sims.photUtils.Photometry import PhotometryStars, PhotometryGalaxies
from lsst.sims.photUtils.EBV import EBVmixin
from lsst.sims.catalogs.measures.instance import compound, InstanceCatalog
from lsst.sims.catalogs.generation.db import DBObject, ObservationMetaData

from lsst.sims.catUtils.baseCatalogModels import *

class sdssGalaxies(InstanceCatalog,EBVmixin,PhotometryGalaxies):
    
    catalog_type = 'sdssGalaxies'
    column_outputs = ['galid','sdss_u','sdss_g','sdss_r','sdss_i','sdss_z',
                      'sdss_uBulge','sdss_gBulge','sdss_rBulge','sdss_iBulge','sdss_zBulge',
                      'sdss_uDisk','sdss_gDisk','sdss_rDisk','sdss_iDisk','sdss_zDisk',
                      'sdss_uAgn','sdss_gAgn','sdss_rAgn','sdss_iAgn','sdss_zAgn']
    
    @compound('sdss_u', 'sdss_g', 'sdss_r', 
              'sdss_i', 'sdss_z',
              'sdss_uBulge', 'sdss_gBulge', 'sdss_rBulge', 'sdss_iBulge', 'sdss_zBulge',
              'sdss_uDisk', 'sdss_gDisk', 'sdss_rDisk', 'sdss_iDisk', 'sdss_zDisk',
              'sdss_uAgn', 'sdss_gAgn', 'sdss_rAgn', 'sdss_iAgn', 'sdss_zAgn')
    def get_all_sdss_mags(self):
        """
        example getter for sdss galaxy magnitudes
        
        bandPassRoot is the root of the names of the files in which
        the bandpasses are stored
        
        """
        idNames = self.column_by_name('galid')
        bandPassNames = ['u','g','r','i','z']
        bandPassDir = os.path.join(os.getenv('THROUGHPUTS_DIR'),'sdss')
        bandPassRoot = 'sdss_'
        
        """
        Here is where we need some code to load a list of bandPass objects
        into self.bandPassList and then call self.setupPhiArray_dict()
        so that the bandPasses are available to the mixin.  Ideally, we
        would only do this once for the whole catalog
        """
        if self.bandPassList is None or self.phiArray is None:
            self.loadBandPassesFromFiles(bandPassNames, 
                     bandPassRoot = bandPassRoot,
                     bandPassDir = bandPassDir)
            self.setupPhiArray_dict()
        
        return self.meta_magnitudes_getter(idNames)



class sdssStars(InstanceCatalog,PhotometryStars):
    
    catalog_type = 'sdssStars'
    column_outputs = ['id','sdss_u','sdss_g','sdss_r','sdss_i','sdss_z']
    
    @compound('sdss_u','sdss_g','sdss_r','sdss_i','sdss_z')
    def get_sdss_magnitudes(self):
        """
        example getter for sdss stellar magnitudes
        
        bandPassRoot is the root of the names of the files in which
        the bandpasses are stored
        """
        idNames = self.column_by_name('id')
        bandPassNames = ['u','g','r','i','z']
        bandPassDir = os.path.join(os.getenv('THROUGHPUTS_DIR'),'sdss')
        bandPassRoot = 'sdss_'
        
        """
        Here is where we need some code to load a list of bandPass objects
        into self.bandPassList and then call self.setupPhiArray_dict()
        so that the bandPasses are available to the mixin.  Ideally, we
        would only do this once for the whole catalog
        """
        if self.bandPassList is None or self.phiArray is None:
            self.loadBandPassesFromFiles(bandPassNames, 
                     bandPassRoot = bandPassRoot,
                     bandPassDir = bandPassDir)
            self.setupPhiArray_dict()
        
        return self.meta_magnitudes_getter(idNames)


obs_metadata_pointed = ObservationMetaData(mjd=2013.23, circ_bounds=dict(ra=200., dec=-30, radius=1.))
obs_metadata_pointed.metadata = {}
obs_metadata_pointed.metadata['Opsim_filter'] = 'i'
dbObj = DBObject.from_objid('rrlystars')
sdssStars = sdssStars(dbObj, obs_metadata = obs_metadata_pointed)
sdssStars.write_catalog("example_sdss_stars.txt")

obs_metadata_pointed = ObservationMetaData(mjd=50000.0, circ_bounds=dict(ra=0., dec=0., radius = 0.01))
obs_metadata_pointed.metadata = {}
obs_metadata_pointed.metadata['Opsim_filter'] = 'i'
dbObj = DBObject.from_objid('galaxyBase')
sdssGalaxies = sdssGalaxies(dbObj, obs_metadata = obs_metadata_pointed)
sdssGalaxies.write_catalog("example_sdss_galaxies.txt")
