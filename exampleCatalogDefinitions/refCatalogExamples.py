"""Instance Catalog"""
import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.catalogs.measures.example_utils import ValidationUtils

class PhotometryMixin(object):
    

    def get_ug_color(self):
        u = self.column_by_name('umag')
        g = self.column_by_name('gmag')
        return u - g

    def get_gr_color(self):
        g = self.column_by_name('gmag')
        r = self.column_by_name('rmag')
        return g - r

    @compound('uRecalc', 'gRecalc', 'rRecalc', 'iRecalc', 'zRecalc', 'yRecalc',
              'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
              'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
              'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn')
    def get_allMags(self):
        retMags, retbMags, retdMags, retaMags = ValidationUtils.calcTotalLSSTMags(self, bandpasses=['u','g','r','i','z','y'])
        return (numpy.array(retMags['u'], dtype=float), 
                numpy.array(retMags['g'], dtype=float),
                numpy.array(retMags['r'], dtype=float),
                numpy.array(retMags['i'], dtype=float),
                numpy.array(retMags['z'], dtype=float),
                numpy.array(retMags['y'], dtype=float),
                numpy.array(retbMags['u'], dtype=float), 
                numpy.array(retbMags['g'], dtype=float),
                numpy.array(retbMags['r'], dtype=float),
                numpy.array(retbMags['i'], dtype=float),
                numpy.array(retbMags['z'], dtype=float),
                numpy.array(retbMags['y'], dtype=float),
                numpy.array(retdMags['u'], dtype=float), 
                numpy.array(retdMags['g'], dtype=float),
                numpy.array(retdMags['r'], dtype=float),
                numpy.array(retdMags['i'], dtype=float),
                numpy.array(retdMags['z'], dtype=float),
                numpy.array(retdMags['y'], dtype=float),
                numpy.array(retaMags['u'], dtype=float), 
                numpy.array(retaMags['g'], dtype=float),
                numpy.array(retaMags['r'], dtype=float),
                numpy.array(retaMags['i'], dtype=float),
                numpy.array(retaMags['z'], dtype=float),
                numpy.array(retaMags['y'], dtype=float))



class AstrometryMixin(object):
    @compound('ra_corr', 'dec_corr')
    def get_points(self):
        return (self.column_by_name('raJ2000') + 0.001,
                self.column_by_name('decJ2000') - 0.001)


class RefCatalogGalaxyBase(InstanceCatalog, AstrometryMixin, PhotometryMixin):
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
