"""Instance Catalog"""
from __future__ import print_function
from builtins import range
import numpy as np
from lsst.sims.utils import SpecMap, defaultSpecMap
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import compound, cached
from lsst.sims.utils import arcsecFromRadians, _observedFromICRS, altAzPaFromRaDec
from lsst.sims.catUtils.mixins import (EBVmixin, PhoSimAstrometryStars,
                                       PhoSimAstrometryGalaxies,
                                       PhoSimAstrometrySSM,
                                       FrozenSNCat)

__all__ = ["write_phoSim_header", "PhosimInputBase",
           "PhoSimCatalogPoint", "PhoSimCatalogZPoint", "PhoSimCatalogSN",
           "PhoSimCatalogSersic2D", "PhoSimCatalogSSM", "PhoSimSpecMap",
           "DefaultPhoSimHeaderMap", 'DefaultPhoSimInstanceCatalogCols']


PhoSimSpecMap = SpecMap(fileDict=defaultSpecMap.fileDict,
                        dirDict={'(^lte)': 'starSED/phoSimMLT'})


# This is a dict of transformations mapping from data in OpSim to header
# values expected by PhoSim.  The dict is keyed on the names of PhoSim
# parameters.  The values of the dict are either straight values, in which
# case those values are written to the PhoSim InstanceCatalog header, or tuples
# containing the name of the OpSim column which should be use to calculate
# the PhoSim parameter and any transformation needed to go from OpSim
# units to PhoSim units.
#
# PhoSim header parameters are documented here
# https://bitbucket.org/phosim/phosim_release/wiki/Instance%20Catalog
#
# OpSim columns are documented here
# https://www.lsst.org/scientists/simulations/opsim/summary-table-column-descriptions-v335

DefaultPhoSimHeaderMap = {'rottelpos': ('rotTelPos', np.degrees),
                          'obshistid': ('obsHistID', None),
                          'moonra': ('moonRA', np.degrees),
                          'moondec': ('moonDec', np.degrees),
                          'moonphase': ('moonPhase', None),
                          'moonalt': ('moonAlt', np.degrees),
                          'dist2moon': ('dist2Moon', np.degrees),
                          'sunalt': ('sunAlt', np.degrees),
                          'seeing': ('rawSeeing', None),
                          'vistime': ('visitExpTime', lambda x: x + 3.0),
                          'nsnap': 2,
                          'seed': ('obsHistID', None)}

# This variable contains all of the columns in a phosim Instance Catalog
# Can be used to reconstruct the information encoded in a phosim instance
# Catalog
DefaultPhoSimInstanceCatalogCols = ('object', 'uniqueID', 'RA', 'DEC'\
                                    , 'MAG_NORM', 'SED_NAME', 'REDSHIFT'\
                                    , 'GAMMA1', 'GAMMA2', 'MU', 'DELTA_RA'\
                                    , 'DELTA_DEC', 'SOURCE_TYPE'\
                                    , 'source_pars', 'DUST_REST_NAME'\
                                    , 'dust_pars_1a', 'dust_pars_1b'\
                                    , 'DUST_LAB_NAME', 'dust_pars_2a'\
                                    , 'dust_pars_2b') 

def evaluate_phosim_header(param, phosim_header_map, obs):
    """
    Find the value of a PhoSim header parameter given a PhoSimHeaderMap,
    and an ObservationMetaData

    @param [in] param is the name of the header parameter desired.

    @param [in] phosim_header_map is a phoSimHeaderMap dict

    @param [in] obs is an ObservationMetaData

    @param [out] returns the value of the PhoSim header parameter
    """

    if param not in phosim_header_map:
        return None

    if isinstance(phosim_header_map[param], tuple):
        transform = phosim_header_map[param]
        if obs.OpsimMetaData is None:
            return None

        if transform[0] not in obs.OpsimMetaData:
            return None

        if transform[1] is None:
            return obs.OpsimMetaData[transform[0]]

        return transform[1](obs.OpsimMetaData[transform[0]])

    return phosim_header_map[param]


def write_phoSim_header(obs, file_handle, phosim_header_map):
    """
    Write the data contained in an ObservationMetaData as a header in a
    PhoSim InstanceCatalog.

    obs is the ObservationMetaData

    file_handle points to the catalog being written.
    """

    if phosim_header_map is None:
        raw_opsim_contents = ""

        if obs.OpsimMetaData is not None:
            sorted_opsim_metadata = list(obs.OpsimMetaData.keys())
            sorted_opsim_metadata.sort()
            for tag in sorted_opsim_metadata:
                raw_opsim_contents += "%s\n" % tag

        raise RuntimeError("You have tried to write a PhoSim InstanceCatalog without specifying "
                           "a phoSimHeaderMap in your call to write_catalog().\n"
                           "\n"
                           "A phoSimHeaderMap is a dict that maps between columns stored in the "
                           "OpSim database from which an ObservationMetaData was created and "
                           "header parameters expected by PhoSim.  The dict is keyed on the names of PhoSim "
                           "parameters.  The values of the dict are either straight values, in which "
                           "case those values are written to the PhoSim InstanceCatalog header, or tuples "
                           "containing the name of the OpSim column which should be use to calculate "
                           "the PhoSim parameter and any transformation needed to go from OpSim "
                           "units to PhoSim units (the transform is 'None' if no transform is needed).\n"
                           "\n"
                           "To add a phoSimHeaderMap, simple assign it to the 'phoSimHeaderMap' "
                           "member variable of your catalog with (for instance):\n"
                           "\n"
                           "myCat.phoSimHeaderMap = my_phosim_header_map\n"
                           "\n"
                           "Before calling write_catalog()\n"
                           "\n"
                           "The header parameters expected by PhoSim can be found at\n"
                           "\n"
                           "https://bitbucket.org/phosim/phosim_release/wiki/Instance%20Catalog\n"
                           "\n"
                           "The contents of the OpSim database's Summary table can be found at\n"
                           "\n"
                           "https://www.lsst.org/scientists/simulations/opsim/summary-table-column-descriptions-v335\n"
                           "\n"
                           "If you do not wish to use any of these columns, you can just pass in an empty "
                           "dict.  If you want to use a pre-made mapping, use the DefaultPhoSimHeaderMap "
                           "imported from lsst.sims.catUtils.exampleCatalogDefinitions\n"
                           "\n"
                           "Note: do not specify ra, dec, alt, az, mjd, filter, or rotSkyPos in your "
                           "phoSimHeaderMap.  These are handled directly by the ObservationMetaData "
                           "to ensure self-consistency.\n"
                           "\n"
                           "For reference, the OpSim columns you can choose to map (i.e. those contained "
                           "in your ObservationMetaData) are:\n\n" +
                           raw_opsim_contents +
                           "\n(Even if your ObservationMetaData contains no extra OpSim Columns "
                           "you may wish to consider adding default PhoSim parameters through "
                           "the phoSimHeaderMap)")

    try:

        # PhoSim wants the MJD at the middle of the visit (i.e. between the two exposures
        # in our two-snap model).  OpSim gives the MJD at the start of the visit.
        # below we calculate the change in MJD necessary to transform the OpSim value
        # into the PhoSim value
        vistime = evaluate_phosim_header('vistime', phosim_header_map, obs)
        if vistime is not None:
            delta_t = 0.5*vistime
        else:
            delta_t = 16.5  # half of the default 33 seconds

        file_handle.write('rightascension %.7f\n' % obs.pointingRA)
        file_handle.write('declination %.7f\n' % obs.pointingDec)
        file_handle.write('mjd %.7f\n' % (obs.mjd.TAI + delta_t/86400.0))
        alt, az, pa = altAzPaFromRaDec(obs.pointingRA, obs.pointingDec, obs, includeRefraction=False)
        file_handle.write('altitude %.7f\n' % alt)
        file_handle.write('azimuth %.7f\n' % az)
        file_handle.write('filter %d\n' %
                          {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4, 'y': 5}[obs.bandpass])

        file_handle.write('rotskypos %.7f\n' % obs.rotSkyPos)
    except:
        print("\n\n")
        print("The ObservationMetaData you tried to write a PhoSim header from")
        print("lacks one of the required parameters")
        print("(pointingRA, pointingDec, mjd, bandpass, rotSkyPos)\n")
        raise

    # sort the header map keys so that PhoSim headers generated with the same
    # map always have parameters in the same order.
    sorted_header_keys = list(phosim_header_map.keys())
    sorted_header_keys.sort()

    for name in sorted_header_keys:
        val = evaluate_phosim_header(name, phosim_header_map, obs)
        if val is not None:
            if isinstance(val, float) or isinstance(val, np.float):
                file_handle.write('%s %.7f\n' % (name, val))
            elif isinstance(val, int) or isinstance(val, np.int):
                file_handle.write('%s %d\n' % (name, val))
            elif isinstance(val, int):
                file_handle.write('%s %ld\n' % (name, val))
            else:
                file_handle.write('%s %s\n' % (name, str(val)))


class PhosimInputBase(InstanceCatalog):

    filtMap = dict([(c, i) for i, c in enumerate('ugrizy')])

    specFileMap = PhoSimSpecMap

    cannot_be_null = ['sedFilepath']

    delimiter = " "

    @cached
    def get_prefix(self):
        chunkiter = range(len(self.column_by_name(self.refIdCol)))
        return np.array(['object' for i in chunkiter], dtype=(str, 6))

    @cached
    def get_sedFilepath(self):
        return np.array([self.specFileMap[k] if k in self.specFileMap else None
                         for k in self.column_by_name('sedFilename')])

    @cached
    def get_spatialmodel(self):
        chunkiter = range(len(self._current_chunk))
        return np.array([self.spatialModel for i in chunkiter], dtype=(str, 8))

    @cached
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

        try:
            string_class = basestring
        except:
            string_class = str

        magNorm = self.column_by_name('magNorm')
        varName = None
        if self.obs_metadata is not None:
            if self.obs_metadata.bandpass is not None:
                if isinstance(self.obs_metadata.bandpass, string_class):
                    varName = 'delta_lsst_' + self.obs_metadata.bandpass

        if varName is not None and varName in self._all_available_columns:
            magNorm_out = magNorm + self.column_by_name(varName)
            return magNorm_out
        else:
            return magNorm

    def write_header(self, file_handle):
        if not hasattr(self, 'phoSimHeaderMap'):
            self.phoSimHeaderMap = None

        write_phoSim_header(self.obs_metadata, file_handle, self.phoSimHeaderMap)


class PhoSimCatalogPoint(PhosimInputBase, PhoSimAstrometryStars):

    catalog_type = 'phoSim_catalog_POINT'

    column_outputs = ['prefix', 'uniqueId', 'raPhoSim', 'decPhoSim', 'phoSimMagNorm', 'sedFilepath',
                      'redshift', 'gamma1', 'gamma2', 'kappa', 'raOffset', 'decOffset',
                      'spatialmodel', 'internalExtinctionModel',
                      'galacticExtinctionModel', 'galacticAv', 'galacticRv']

    default_columns = [('redshift', 0., float), ('gamma1', 0., float), ('gamma2', 0., float),
                       ('kappa', 0., float), ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticExtinctionModel', 'CCM', (str, 3)), ('galacticRv', 3.1, float),
                       ('internalExtinctionModel', 'none', (str, 4))]

    default_formats = {'S': '%s', 'f': '%.9g', 'i': '%i'}

    spatialModel = "point"

    transformations = {'raPhoSim': np.degrees, 'decPhoSim': np.degrees}


class PhoSimCatalogZPoint(PhosimInputBase, PhoSimAstrometryGalaxies, EBVmixin):

    catalog_type = 'phoSim_catalog_ZPOINT'

    column_outputs = ['prefix', 'uniqueId', 'raPhoSim', 'decPhoSim', 'phoSimMagNorm', 'sedFilepath',
                      'redshift', 'gamma1', 'gamma2', 'kappa', 'raOffset', 'decOffset',
                      'spatialmodel', 'internalExtinctionModel',
                      'galacticExtinctionModel', 'galacticAv', 'galacticRv']

    default_columns = [('gamma1', 0., float), ('gamma2', 0., float), ('kappa', 0., float),
                       ('raOffset', 0., float), ('decOffset', 0., float),
                       ('spatialmodel', 'ZPOINT', (str, 6)),
                       ('galacticExtinctionModel', 'CCM', (str, 3)),
                       ('galacticAv', 0.1, float), ('galacticRv', 3.1, float),
                       ('internalExtinctionModel', 'none', (str, 4))]

    default_formats = {'S': '%s', 'f': '%.9g', 'i': '%i'}

    spatialModel = "point"

    transformations = {'raPhoSim': np.degrees, 'decPhoSim': np.degrees}


class PhoSimCatalogSN(PhoSimCatalogZPoint, FrozenSNCat, EBVmixin):
    """
    Mixin for PhoSim Instance Catalogs in PhoSim
    """
    catalog_type = 'phoSim_SNcatalog'
    writeSEDFile = True

    def get_sedFilepath(self):
        return self.column_by_name('TsedFilepath')

    def get_phoSimMagNorm(self):
        return self.column_by_name('magNorm')


class PhoSimCatalogSersic2D(PhoSimCatalogZPoint):

    catalog_type = 'phoSim_catalog_SERSIC2D'

    column_outputs = ['prefix', 'uniqueId', 'raPhoSim', 'decPhoSim', 'phoSimMagNorm', 'sedFilepath',
                      'redshift', 'gamma1', 'gamma2', 'kappa', 'raOffset', 'decOffset',
                      'spatialmodel', 'majorAxis', 'minorAxis', 'positionAngle', 'sindex',
                      'internalExtinctionModel', 'internalAv', 'internalRv',
                      'galacticExtinctionModel', 'galacticAv', 'galacticRv',]

    default_columns = [('gamma1', 0., float), ('gamma2', 0., float), ('kappa', 0., float),
                       ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticAv', 0.1, float), ('galacticRv', 3.1, float),
                       ('galacticExtinctionModel', 'CCM', (str, 3)),
                       ('internalExtinctionModel', 'CCM', (str, 3)), ('internalAv', 0., float),
                       ('internalRv', 3.1, float)]

    default_formats = {'S': '%s', 'f': '%.9g', 'i': '%i'}

    spatialModel = "sersic2d"

    transformations = {'raPhoSim': np.degrees, 'decPhoSim': np.degrees, 'positionAngle': np.degrees,
                       'majorAxis': arcsecFromRadians, 'minorAxis': arcsecFromRadians}

    @compound('raPhoSim', 'decPhoSim')
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

    column_outputs = ['prefix', 'uniqueId', 'raPhoSim', 'decPhoSim', 'phoSimMagNorm', 'sedFilepath',
                      'redshift', 'gamma1', 'gamma2', 'kappa', 'raOffset', 'decOffset',
                      'spatialmodel', 'internalExtinctionModel', 'galacticExtinctionModel']

    default_columns = [('redshift', 0., float), ('gamma1', 0., float), ('gamma2', 0., float),
                       ('kappa', 0., float), ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticExtinctionModel', 'none', (str, 4)),
                       ('internalExtinctionModel', 'none', (str, 4))]

    default_formats = {'S': '%s', 'f': '%.9g', 'i': '%i'}

    spatialModel = "point"

    transformations = {'raPhoSim': np.degrees, 'decPhoSim': np.degrees}
