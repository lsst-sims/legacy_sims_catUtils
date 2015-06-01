"""Instance Catalog"""
import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.catUtils import ValidationUtils
from lsst.sims.coordUtils.Astrometry import AstrometryGalaxies
from lsst.sims.photUtils.Photometry import PhotometryGalaxies, PhotometryStars

__all__ = ["RefCatalogGalaxyBase", "GalaxyPhotometry", "RefCatalogStarBase"]

class RefCatalogGalaxyBase(InstanceCatalog, AstrometryGalaxies, PhotometryGalaxies):
    comment_char = ''
    catalog_type = 'ref_catalog_galaxy'
    column_outputs = ['uniqueId', 'objId', 'galid', 'meanRaJ2000', 'meanDecJ2000', 'redshift',
                      'majorAxisDisk', 'minorAxisDisk', 'positionAngleDisk',
                      'majorAxisBulge', 'minorAxisBulge', 'positionAngleBulge',
                      'DiskLSSTu', 'DiskLSSTg', 'DiskLSSTr', 'DiskLSSTi', 'DiskLSSTz', 'DiskLSSTy',
                      'BulgeLSSTu', 'BulgeLSSTg', 'BulgeLSSTr', 'BulgeLSSTi', 'BulgeLSSTz', 'BulgeLSSTy',
                      'u_ab', 'g_ab', 'r_ab', 'i_ab', 'z_ab', 'y_ab']
    default_formats = {'S':'%s', 'f':'%.8f', 'i':'%i'}
    transformations = {'meanRaJ2000':numpy.degrees, 'meanDecJ2000':numpy.degrees,
                       'majorAxisDisk':arcsecFromRadians, 'minorAxisDisk':arcsecFromRadians,
                       'positionAngleDisk':numpy.degrees,
                       'majorAxisBulge':arcsecFromRadians, 'minorAxisBulge':arcsecFromRadians,
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
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
                      'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
                      'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn',
                      'u_ab', 'g_ab', 'r_ab', 'i_ab', 'z_ab', 'y_ab']

    transformations = {'majorAxisDisk':arcsecFromRadians, 'minorAxisDisk':arcsecFromRadians,
                       'majorAxisBulge':arcsecFromRadians, 'minorAxisBulge':arcsecFromRadians}

class RefCatalogStarBase(InstanceCatalog, AstrometryGalaxies, PhotometryGalaxies):
    comment_char = ''
    catalog_type = 'ref_catalog_star'
    column_outputs = ['uniqueId', 'objId', 'id', 'meanRaJ2000', 'meanDecJ2000',
                      'umag', 'gmag', 'rmag', 'imag', 'zmag', 'ymag']
    default_formats = {'S':'%s', 'f':'%.8f', 'i':'%i'}
    transformations = {'meanRaJ2000':numpy.degrees, 'meanDecJ2000':numpy.degrees}

    def get_objectId(self):
        return self.column_by_name(self.refIdCol)

    def get_meanRaJ2000(self):
        return self.column_by_name('raJ2000')

    def get_meanDecJ2000(self):
        return self.column_by_name('decJ2000')
