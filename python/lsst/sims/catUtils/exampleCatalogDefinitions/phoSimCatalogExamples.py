"""Instance Catalog"""
import numpy
from lsst.sims.utils import SpecMap, defaultSpecMap
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.utils import arcsecFromRadians, _observedFromICRS, altAzPaFromRaDec
from lsst.sims.catUtils.mixins import AstrometryStars, AstrometryGalaxies, \
                                      AstrometrySSM, EBVmixin, PhoSimAstrometryStars, \
                                      PhoSimAstrometryGalaxies, PhoSimAstrometrySSM

__all__ = ["write_phoSim_header", "PhosimInputBase",
           "PhoSimCatalogPoint", "PhoSimCatalogZPoint",
           "PhoSimCatalogSersic2D", "PhoSimCatalogSSM", "PhoSimSpecMap"]


PhoSimSpecMap = SpecMap(fileDict=defaultSpecMap.fileDict,
                        dirDict={'(^lte)':'starSED/phoSimMLT'})


def write_phoSim_header(obs, file_handle):
    """
    Write the data contained in an ObservationMetaData as a header in a
    PhoSim InstanceCatalog.

    obs is the ObservationMetaData

    file_handle points to the catalog being written.
    """
    try:
        file_handle.write('Unrefracted_RA %.9g\n' % obs.pointingRA)
        file_handle.write('Unrefracted_Dec %.9g\n' % obs.pointingDec)
        file_handle.write('Opsim_expmjd %.9g\n' % obs.mjd.TAI)
        alt, az, pa = altAzPaFromRaDec(obs.pointingRA, obs.pointingDec, obs)
        file_handle.write('Opsim_altitude %.9g\n' % alt)
        file_handle.write('Opsim_azimuth %.9g\n' % az)
        airmass = 1.0/numpy.cos(0.5*numpy.pi-numpy.radians(alt))
        file_handle.write('airmass %.9g\n' % airmass)
        file_handle.write('Opsim_filter %d\n' %
                      {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}[obs.bandpass])

        file_handle.write('Opsim_rotskypos %.9g\n' % obs.rotSkyPos)
    except:
        print "The ObservationMetaData you tried to write a PhoSim header from"
        print "lacks one of the required parameters"
        print "(pointingRA, pointingDec, mjd, bandpass, rotSkyPos)"
        raise


    already_written = ('Unrefracted_RA', 'Unrefracted_Dec',
                       'Opsim_expmjd', 'Opsim_altitude',
                       'Opsim_azimuth', 'airmass', 'Opsim_filter',
                       'Opsim_rotskypos')

    for kk in obs._phoSimMetadata:
        if kk in already_written:
            raise RuntimeError("This ObservationMetaData has conflicting values for "
                               "%s" % kk)

        file_handle.write('%s %.9g\n' % (kk, obs._phoSimMetadata[kk]))



class PhosimInputBase(InstanceCatalog):

    filtMap = dict([(c, i) for i,c in enumerate('ugrizy')])

    specFileMap = PhoSimSpecMap

    cannot_be_null = ['sedFilepath']

    delimiter = " "


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
        write_phoSim_header(self.obs_metadata, file_handle)


class PhoSimCatalogPoint(PhosimInputBase, PhoSimAstrometryStars, EBVmixin):

    catalog_type = 'phoSim_catalog_POINT'

    column_outputs = ['prefix', 'uniqueId','raPhoSim','decPhoSim','phoSimMagNorm','sedFilepath',
                      'redshift','shear1','shear2','kappa','raOffset','decOffset',
                      'spatialmodel','galacticExtinctionModel','galacticAv','galacticRv',
                      'internalExtinctionModel']

    default_columns = [('redshift', 0., float),('shear1', 0., float), ('shear2', 0., float),
                       ('kappa', 0., float), ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticExtinctionModel', 'CCM', (str,3)), ('galacticRv', 3.1, float),
                       ('internalExtinctionModel', 'none', (str,4))]

    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}

    spatialModel = "point"

    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees}


class PhoSimCatalogZPoint(PhosimInputBase, PhoSimAstrometryGalaxies, EBVmixin):

    catalog_type = 'phoSim_catalog_ZPOINT'

    column_outputs = ['prefix', 'uniqueId','raPhoSim','decPhoSim','phoSimMagNorm','sedFilepath',
                      'redshift','shear1','shear2','kappa','raOffset','decOffset',
                      'spatialmodel','galacticExtinctionModel','galacticAv','galacticRv',
                      'internalExtinctionModel']

    default_columns = [('shear1', 0., float), ('shear2', 0., float), ('kappa', 0., float),
                       ('raOffset', 0., float), ('decOffset', 0., float), ('spatialmodel', 'ZPOINT', (str, 6)),
                       ('galacticExtinctionModel', 'CCM', (str,3)),
                       ('galacticAv', 0.1, float), ('galacticRv', 3.1, float),
                       ('internalExtinctionModel', 'none', (str,4))]

    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}

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
                       ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticAv', 0.1, float), ('galacticRv', 3.1, float),
                       ('galacticExtinctionModel', 'CCM', (str,3)),
                       ('internalExtinctionModel', 'CCM', (str,3)), ('internalAv', 0., float),
                       ('internalRv', 3.1, float) ]

    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}

    spatialModel = "sersic2d"

    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees, 'positionAngle':numpy.degrees,
                       'majorAxis':arcsecFromRadians, 'minorAxis':arcsecFromRadians}


    @compound('raPhoSim','decPhoSim')
    def get_phoSimCoordinates(self):
        """Getter for RA, Dec coordinates expected by PhoSim.

        These are observed RA, Dec coordinates with the effects of nutation, aberration,
        and precession subtracted out by the PhosimInputBase._dePrecess() method.
        This preserves the relative effects of nutation, aberration, and precession while
        re-aligning the catalog with the boresite RA, Dec so that astrometric solutions
        make sense."""

        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')
        raObs, decObs = _observedFromICRS(ra, dec, includeRefraction = False, obs_metadata=self.obs_metadata,
                                          epoch=self.db_obj.epoch)

        return self._dePrecess(raObs, decObs, self.obs_metadata)


class PhoSimCatalogSSM(PhosimInputBase, PhoSimAstrometrySSM):

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
