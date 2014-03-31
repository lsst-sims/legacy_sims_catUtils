"""Instance Catalog"""
import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.catalogs.measures.example_utils import ValidationUtils
from lsst.sims.catalogs.measures.astrometry.Astrometry import Astrometry
from lsst.sims.catalogs.measures.photometry.Photometry import PhotometryGalaxies, PhotometryStars

class RefCatalogGalaxyBase(InstanceCatalog, Astrometry, PhotometryGalaxies):
    comment_char = ''
    catalog_type = 'ref_catalog_galaxy'
    column_outputs = ['uniqueId', 'objId', 'galid', 'meanRaJ2000', 'meanDecJ2000', 'redshift',
                      'majorAxisDisk', 'minorAxisDisk', 'positionAngleDisk',
                      'majorAxisBulge', 'minorAxisBulge', 'positionAngleBulge', 
                      'DiskLSSTu', 'DiskLSSTg', 'DiskLSSTr', 'DiskLSSTi', 'DiskLSSTz', 'DiskLSSTy',
                      'BulgeLSSTu', 'BulgeLSSTg', 'BulgeLSSTr', 'BulgeLSSTi', 'BulgeLSSTz', 'BulgeLSSTy',
                      'u_ab', 'g_ab', 'r_ab', 'i_ab', 'z_ab', 'y_ab'] 
    default_formats = {'S':'%s', 'f':'%.8f', 'i':'%i'}
    transformations = {'meanRaJ2000':numpy.degrees, 'meanDecJ2000':numpy.degrees, 'majorAxisDisk':numpy.degrees, 'minorAxisDisk':numpy.degrees,
                       'positionAngleDisk':numpy.degrees, 'majorAxisBulge':numpy.degrees, 'minorAxisBulge':numpy.degrees, 
                       'positionAngleBulge':numpy.degrees}

    def get_objectId(self): 
        return self.column_by_name(self.refIdCol)

    def get_meanRaJ2000(self):
        return self.column_by_name('raJ2000')

    def get_meanDecJ2000(self):
        return self.column_by_name('decJ2000')

class GalaxyPhotometry(RefCatalogGalaxyBase):
    catalog_type = 'galaxy_photometry_cat'
    column_outputs = ['uniqueId', 'objId', 'galid', 'meanRaJ2000', 'meanDecJ2000', 'redshift',
                      'majorAxisDisk', 'minorAxisDisk', 'positionAngleDisk',
                      'majorAxisBulge', 'minorAxisBulge', 'positionAngleBulge', 
		      'internalAvBulge', 'internalAvDisk',
                      'uRecalc', 'gRecalc', 'rRecalc', 'iRecalc', 'zRecalc', 'yRecalc',
                      'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
                      'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
                      'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn',
                      'u_ab', 'g_ab', 'r_ab', 'i_ab', 'z_ab', 'y_ab'] 
