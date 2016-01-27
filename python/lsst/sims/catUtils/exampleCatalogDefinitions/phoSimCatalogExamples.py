"""Instance Catalog"""
import numpy
from lsst.sims.utils import SpecMap, defaultSpecMap
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.catUtils.mixins import AstrometryStars, AstrometryGalaxies, \
                                      EBVmixin

__all__ = ["PhosimInputBase", "PhoSimCatalogPoint", "PhoSimCatalogZPoint",
           "PhoSimCatalogSersic2D", "PhoSimSpecMap"]


PhoSimSpecMap = SpecMap(fileDict=defaultSpecMap.fileDict,
                        dirDict={'(^lte)':'starSED/phoSimMLT'})


class PhosimInputBase(InstanceCatalog):
    filtMap = dict([(c, i) for i,c in enumerate('ugrizy')])

    specFileMap = PhoSimSpecMap

    cannot_be_null = ['sedFilepath']

    headerTransformations = {'pointingRA':numpy.degrees, 'pointingDec':numpy.degrees,
                       'Opsim_moonra':numpy.degrees, 'Opsim_moondec':numpy.degrees,
                       'Opsim_rotskypos':numpy.degrees, 'Opsim_rottelpos':numpy.degrees,
                       'Opsim_sunalt':numpy.degrees, 'Opsim_moonalt':numpy.degrees,
                       'Opsim_dist2moon':numpy.degrees, 'Opsim_altitude':numpy.degrees,
                       'Opsim_azimuth':numpy.degrees, 'Opsim_filter':filtMap.get,
                       'Unrefracted_Altitude':numpy.degrees, 'Unrefracted_Azimuth':numpy.degrees}


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

    def get_phoSimMagNorm(self):
        """
        This getter returns the magnitude normalization expected by PhoSim (the magnitude at
        500 nm).

        To account for variability, the getter adds the delta_lsst_x column from the Variability
        mixin where 'x' is the bandpass defined by self.observation_metadata.bandpass (assuming
        that self.observation_metadata.bandpass is not list-like; if it is list-like, then no
        variability is added to the magNorm value).

        Redshift is currently ignored.  That may or may  not be appropriate.  This requires
        further investigation into the behavior of PhoSim.
        """

        magNorm = self.column_by_name('magNorm')
        varName = None
        if self.obs_metadata is not None:
            if self.obs_metadata.bandpass is not None:
                if not hasattr(self.obs_metadata.bandpass, '__iter__'):
                    varName = 'delta_lsst_' + self.obs_metadata.bandpass

        if varName is not None and varName in self._all_available_columns:
            magNorm += self.column_by_name(varName)

        return magNorm


    def write_header(self, file_handle):
        header_name_map = {'pointingRA': 'Unrefracted_RA', 'pointingDec': 'Unrefracted_Dec'}
        md = self.obs_metadata.phoSimMetaData
        if md is None:
            raise RuntimeError("Can't write a phoSim catalog without a full phoSimMetaData dictionary")
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

            if k in header_name_map:
                k_out = header_name_map[k]
            else:
                k_out = k

            file_handle.write(templ%(k_out, outval)+"\n")

class PhoSimCatalogPoint(PhosimInputBase, AstrometryStars, EBVmixin):
    catalog_type = 'phoSim_catalog_POINT'
    column_outputs = ['prefix', 'uniqueId','raPhoSim','decPhoSim','phoSimMagNorm','sedFilepath',
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


class PhoSimCatalogZPoint(PhosimInputBase, AstrometryGalaxies, EBVmixin):
    catalog_type = 'phoSim_catalog_ZPOINT'
    column_outputs = ['prefix', 'uniqueId','raPhoSim','decPhoSim','phoSimMagNorm','sedFilepath',
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
    column_outputs = ['prefix', 'uniqueId','raPhoSim','decPhoSim','phoSimMagNorm','sedFilepath',
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
    'majorAxis':arcsecFromRadians, 'minorAxis':arcsecFromRadians}

