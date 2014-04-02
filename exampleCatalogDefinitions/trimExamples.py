"""Instance Catalog"""
import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound, cached
from lsst.sims.coordUtils.Astrometry import Astrometry
from lsst.sims.photUtils.Photometry import PhotometryStars
from lsst.sims.photUtils.EBV import EBVmixin

class TrimCatalogPoint(InstanceCatalog, Astrometry, PhotometryStars,EBVmixin):
    catalog_type = 'trim_catalog_POINT'
    column_outputs = ['prefix', 'uniqueId','raTrim','decTrim','magNorm','sedFilepath',
                      'redshift','shear1','shear2','kappa','raOffset','decOffset',
                      'spatialmodel','galacticExtinctionModel','galacticAv','galacticRv',
                      'internalExtinctionModel']
    default_columns = [('redshift', 0., float),('shear1', 0., float), ('shear2', 0., float), 
                       ('kappa', 0., float), ('raOffset', 0., float), ('decOffset', 0., float), 
                       ('galacticExtinctionModel', 'CCM', (str,3)),
                       ('internalExtinctionModel', 'none', (str,4))]
    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}
    delimiter = " "
    filtMap = dict([(c, i) for i,c in enumerate('ugrizy')])
    transformations = {'raTrim':numpy.degrees, 'decTrim':numpy.degrees}
    headerTransformations = {'Unrefracted_RA':numpy.degrees, 'Unrefracted_Dec':numpy.degrees, 
                       'Opsim_moonra':numpy.degrees, 'Opsim_moondec':numpy.degrees, 
                       'Opsim_rotskypos':numpy.degrees, 'Opsim_rottelpos':numpy.degrees, 
                       'Opsim_sunalt':numpy.degrees, 'Opsim_moonalt':numpy.degrees, 
                       'Opsim_dist2moon':numpy.degrees, 'Opsim_altitude':numpy.degrees, 
                       'Opsim_azimuth':numpy.degrees, 'Opsim_filter':filtMap.get}

    def get_prefix(self):
        chunkiter = xrange(len(self.column_by_name(self.refIdCol)))
        return numpy.array(['object' for i in chunkiter], dtype=(str, 6))
    def get_sedFilepath(self):
        return numpy.array([self.specFileMap[k] 
                         for k in self.column_by_name('sedFilename')])
    def get_raTrim(self):
        return self.column_by_name('ra_corr')
    def get_decTrim(self):
        return self.column_by_name('dec_corr')
    def get_spatialmodel(self):
        chunkiter = xrange(len(self._current_chunk))
        return numpy.array([self.db_obj.getSpatialModel() for i in
               chunkiter], dtype=(str, 7))
                        
    def write_header(self, file_handle):
        md = self.obs_metadata.metadata
        if md is None:
            raise RuntimeError("Can't write a trim without a full metadata dictionary")
        for k in md:
            typ = md[k][1].kind
            templ = self.default_formats.get(typ, None)
            if templ is None:
                warnings.warn("Using raw formatting for header key %s "+\
                              "with type %s" % (k, typ))
                templ = "%s"
            templ = "%s "+templ
            if k in self.headerTransformations.keys():
                outval = self.headerTransformations[k](md[k][0])
            else:
                outval = md[k][0]
            file_handle.write(templ%(k, outval)+"\n") 


class TrimCatalogZPoint(TrimCatalogPoint, Astrometry, PhotometryStars,EBVmixin):
    catalog_type = 'trim_catalog_ZPOINT'
    column_outputs = ['prefix', 'uniqueId','raTrim','decTrim','magNorm','sedFilepath',
                      'redshift','shear1','shear2','kappa','raOffset','decOffset',
                      'spatialmodel','galacticExtinctionModel','galacticAv','galacticRv',
                      'internalExtinctionModel']
    default_columns = [('shear1', 0., float), ('shear2', 0., float), ('kappa', 0., float),
                       ('raOffset', 0., float), ('decOffset', 0., float), ('spatialmodel', 'ZPOINT', (str, 6)),
                       ('galacticExtinctionModel', 'CCM', (str,3)),
                       ('galacticAv', 0.1, float),
                       ('internalExtinctionModel', 'none', (str,4))]
    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}
    delimiter = " "
    transformations = {'raTrim':numpy.degrees, 'decTrim':numpy.degrees}


class TrimCatalogSersic2D(TrimCatalogZPoint, Astrometry, PhotometryStars,EBVmixin):
    catalog_type = 'trim_catalog_SERSIC2D'
    column_outputs = ['prefix', 'uniqueId','objId','raTrim','decTrim','magNorm','sedFilepath',
                      'redshift','shear1','shear2','kappa','raOffset','decOffset',
                      'spatialmodel','majorAxis','minorAxis','positionAngle','sindex',
                      'galacticExtinctionModel','galacticAv','galacticRv',
                      'internalExtinctionModel','internalAv','internalRv']
    default_columns = [('shear1', 0., float), ('shear2', 0., float), ('kappa', 0., float),
                       ('raOffset', 0., float), ('decOffset', 0., float), ('galacticAv', 0.1, float),
                       ('galacticExtinctionModel', 'CCM', (str,3)),
                       ('internalExtinctionModel', 'CCM', (str,3)), ('internalAv', 0., float),
                       ('internalRv', 3.1, float) ]
    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}
    delimiter = " "
    transformations = {'raTrim':numpy.degrees, 'decTrim':numpy.degrees, 'positionAngle':numpy.degrees, 
    'majorAxis':numpy.degrees, 'minorAxis':numpy.degrees} 

