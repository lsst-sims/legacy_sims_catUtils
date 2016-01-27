"""Instance Catalog"""
import numpy
from lsst.sims.utils import SpecMap, defaultSpecMap
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.utils import arcsecFromRadians, _observedFromICRS
from lsst.sims.catUtils.mixins import AstrometryStars, AstrometryGalaxies, \
                                      AstrometrySSM, EBVmixin

from lsst.sims.utils import sphericalFromCartesian, cartesianFromSpherical
from lsst.sims.utils import rotationMatrixFromVectors

__all__ = ["PhosimInputBase","PhoSimAstrometryBase",
           "PhoSimCatalogPoint", "PhoSimCatalogZPoint",
           "PhoSimCatalogSersic2D", "PhoSimCatalogSSM", "PhoSimSpecMap"]


PhoSimSpecMap = SpecMap(fileDict=defaultSpecMap.fileDict,
                        dirDict={'(^lte)':'starSED/phoSimMLT'})


class PhoSimAstrometryBase(object):

    def _dePrecess(self, ra_in, dec_in, obs_metadata):
        """
        Calculate the displacement between the boresite and the boresite
        corrected for precession, nutation, and aberration (not refraction).

        Convert boresite and corrected boresite to Cartesian coordinates.

        Calculate the rotation matrix to go between those Cartesian vectors.

        Convert [ra_in, dec_in] into Cartesian coordinates.

        Apply the rotation vector to those Cartesian coordinates.

        Convert back to ra, dec-like coordinates

        @param [in] ra_in is a numpy array of RA in radians

        @param [in] dec_in is a numpy array of Dec in radians

        @param [in] obs_metadata is an ObservationMetaData

        @param [out] ra_out is a numpy array of de-precessed RA in radians

        @param [out] dec_out is a numpy array of de-precessed Dec in radians
        """

        if len(ra_in)==0:
            return numpy.array([[],[]])

        xyz_bore = cartesianFromSpherical(numpy.array([obs_metadata._pointingRA]),
                                          numpy.array([obs_metadata._pointingDec]))

        precessedRA, precessedDec = _observedFromICRS(numpy.array([obs_metadata._pointingRA]),
                                                      numpy.array([obs_metadata._pointingDec]),
                                                      obs_metadata=obs_metadata, epoch=2000.0,
                                                      includeRefraction=False)

        xyz_precessed = cartesianFromSpherical(precessedRA, precessedDec)

        norm = numpy.sqrt(numpy.power(xyz_bore[0],2).sum())
        xyz_bore = xyz_bore/norm

        norm = numpy.sqrt(numpy.power(xyz_precessed[0],2).sum())
        xyz_precessed = xyz_precessed/norm

        rotMat = rotationMatrixFromVectors(xyz_precessed[0], xyz_bore[0])

        xyz_list = cartesianFromSpherical(ra_in, dec_in)

        xyz_de_precessed = numpy.array([numpy.dot(rotMat, xx) for xx in xyz_list])
        ra_deprecessed, dec_deprecessed = sphericalFromCartesian(xyz_de_precessed)
        return numpy.array([ra_deprecessed, dec_deprecessed])


class PhosimInputBase(InstanceCatalog, PhoSimAstrometryBase):

    filtMap = dict([(c, i) for i,c in enumerate('ugrizy')])

    specFileMap = PhoSimSpecMap

    cannot_be_null = ['sedFilepath']

    delimiter = " "

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
            magNorm_out = magNorm + self.column_by_name(varName)
            return magNorm_out
        else:
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

    spatialModel = "point"

    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees}


    @compound('raPhoSim','decPhoSim')
    def get_phoSimCoordinates(self):
        """Getter for unrefracted observed RA, Dec as expected by PhoSim"""
        raObs, decObs = self.observedStellarCoordinates(includeRefraction = False)
        return self._dePrecess(raObs, decObs, self.obs_metadata)


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

    spatialModel = "point"

    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees}


    @compound('raPhoSim','decPhoSim')
    def get_phoSimCoordinates(self):
        """Getter for unrefracted observed RA, Dec as expected by PhoSim"""
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')
        raObs, decObs = _observedFromICRS(ra, dec, includeRefraction = False, obs_metadata=self.obs_metadata,
                                          epoch=self.db_obj.epoch)

        return self._dePrecess(raObs, decObs, self.obs_metadata)


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

    spatialModel = "sersic2d"

    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees, 'positionAngle':numpy.degrees,
                       'majorAxis':arcsecFromRadians, 'minorAxis':arcsecFromRadians}


    @compound('raPhoSim','decPhoSim')
    def get_phoSimCoordinates(self):
        """Getter for unrefracted observed RA, Dec as expected by PhoSim"""
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')
        raObs, decObs = _observedFromICRS(ra, dec, includeRefraction = False, obs_metadata=self.obs_metadata,
                                          epoch=self.db_obj.epoch)

        return self._dePrecess(raObs, decObs, self.obs_metadata)


class PhoSimCatalogSSM(PhosimInputBase, AstrometrySSM):

    catalog_type = 'phoSim_catalog_SSM'

    column_outputs = ['prefix', 'uniqueId','raPhoSim','decPhoSim','phoSimMagNorm','sedFilepath',
                      'redshift','shear1','shear2','kappa','raOffset','decOffset',
                      'spatialmodel','galacticExtinctionModel', 'internalExtinctionModel']

    default_columns = [('redshift', 0., float),('shear1', 0., float), ('shear2', 0., float),
                       ('kappa', 0., float), ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticExtinctionModel', 'none', (str,3)),
                       ('internalExtinctionModel', 'none', (str,4))]

    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}

    spatialModel = "point"

    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees}

    @compound('raPhoSim', 'decPhoSim')
    def get_phoSimCoordinates(self):
        raObs, decObs = self.observedSSMCoordinates(includeRefraction = False)
        return self._dePrecess(raObs, decObs, self.obs_metadata)
