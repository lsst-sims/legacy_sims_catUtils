"""Instance Catalog"""
import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.utils import radiansToArcsec
from lsst.sims.coordUtils.Astrometry import AstrometryStars, AstrometryGalaxies
from lsst.sims.photUtils.Photometry import PhotometryStars, PhotometryGalaxies
from lsst.sims.photUtils.EBV import EBVmixin

__all__ = ["PhosimInputBase", "PhoSimCatalogPoint", "PhoSimCatalogZPoint",
           "PhoSimCatalogSersic2D"]

class PhosimInputBase(InstanceCatalog):
    filtMap = dict([(c, i) for i,c in enumerate('ugrizy')])

    cannot_be_null = ['sedFilepath']

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
        return numpy.array([self.specFileMap[k] if self.specFileMap.has_key(k) else None
                         for k in self.column_by_name('sedFilename')])
    def get_spatialmodel(self):
        chunkiter = xrange(len(self._current_chunk))
        return numpy.array([self.spatialModel for i in
               chunkiter], dtype=(str, 8))

    def write_header(self, file_handle):
        md = self.obs_metadata.phoSimMetadata
        if md is None:
            raise RuntimeError("Can't write a phoSim catalog without a full phoSimMetadata dictionary")
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

class PhoSimCatalogPoint(PhosimInputBase, AstrometryStars, PhotometryStars, EBVmixin):
    catalog_type = 'phoSim_catalog_POINT'
    column_outputs = ['prefix', 'uniqueId','raPhoSim','decPhoSim','magNorm','sedFilepath',
                      'redshift','shear1','shear2','kappa','raOffset','decOffset',
                      'spatialmodel','galacticExtinctionModel','galacticAv','galacticRv',
                      'internalExtinctionModel']
    default_columns = [('redshift', 0., float),('shear1', 0., float), ('shear2', 0., float),
                       ('kappa', 0., float), ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticExtinctionModel', 'CCM', (str,3)),
                       ('internalExtinctionModel', 'none', (str,4))]
    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}
    delimiter = " "
    spatialModel = "point"
    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees}


class PhoSimCatalogZPoint(PhosimInputBase, AstrometryGalaxies, PhotometryGalaxies, EBVmixin):
    catalog_type = 'phoSim_catalog_ZPOINT'
    column_outputs = ['prefix', 'uniqueId','raPhoSim','decPhoSim','magNorm','sedFilepath',
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
    spatialModel = "point"
    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees}


class PhoSimCatalogSersic2D(PhoSimCatalogZPoint):
    catalog_type = 'phoSim_catalog_SERSIC2D'
    column_outputs = ['prefix', 'uniqueId','raPhoSim','decPhoSim','magNorm','sedFilepath',
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
    spatialModel = "sersic2d"
    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees, 'positionAngle':numpy.degrees,
    'majorAxis':radiansToArcsec, 'minorAxis':radiansToArcsec}

