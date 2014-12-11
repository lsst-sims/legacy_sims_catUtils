import os
import eups
from lsst.sims.photUtils import Bandpass, Sed
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catUtils.baseCatalogModels import StarObj, GalaxyAgnObj, \
                                                 GalaxyDiskObj, GalaxyBulgeObj

__all__ = ["calcADUwrapper", "SearchReversion", "testGalaxyBulge",
           "testGalaxyDisk", "testGalaxyAgn", "testStars"]


def calcADUwrapper(sedName=None, magNorm=None, redshift=None, internalAv=None, internalRv=None,
                   galacticAv=None, galacticRv=None, bandpass=None):

    imsimband = Bandpass()
    imsimband.imsimBandpass()
    sedDir = eups.productDir('sims_sed_library')
    sedFile = os.path.join(sedDir, sedName)
    sed = Sed()
    sed.readSED_flambda(sedFile)
    fNorm = sed.calcFluxNorm(magNorm, imsimband)
    sed.multiplyFluxNorm(fNorm)
    if internalAv is not None and internalRv is not None:
        if internalAv != 0.0 and internalRv != 0.0:
            a_int, b_int = sed.setupCCMab()
            sed.addCCMDust(a_int, b_int, A_v=internalAv, R_v=internalRv)
    
    if redshift is not None and redshift!=0.0:
        sed.redshiftSED(redshift, dimming=False)
    
    a_int, b_int = sed.setupCCMab()
    sed.addCCMDust(a_int, b_int, A_v=galacticAv, R_v=galacticRv)
    
    adu = sed.calcADU(bandpass)
    
    return adu


class SearchReversion(CatalogDBObject):
    """
    This is a mixin which is used below to force the galaxy CatalogDBObjects created for
    this unittest to use the methods defined in the CatalogDBObject class.  This is because
    we are using classes which inherit from GalaxyTileObj but do not actually want to use
    the tiled query routines.

    We also use this mixin for our stellar database object.  This is because StarObj
    implements a query search based on htmid, which the test database for this unit
    test will not have.
    """

    def _get_column_query(self, *args, **kwargs):
        return CatalogDBObject._get_column_query(self,*args, **kwargs)

    def _final_pass(self, *args, **kwargs):
        return CatalogDBObject._final_pass(self,*args, **kwargs)

    def query_columns(self, *args, **kwargs):
        return CatalogDBObject.query_columns(self, *args, **kwargs)

class testGalaxyBulge(SearchReversion, GalaxyBulgeObj):
    """
    A class for storing galaxy bulges
    """
    objid = 'testBulgeDBObj'
    objectTypeId = 88

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject daughter classes assume that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = GalaxyBulgeObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testGalaxyDisk(SearchReversion, GalaxyDiskObj):
    objid = 'testDiskDBObj'
    objectTypeId = 89

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject daughter classes assume that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = GalaxyDiskObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testGalaxyAgn(SearchReversion, GalaxyAgnObj):
    objid = 'testAgnDBObj'
    objectTypeId = 90

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject daughter classes assume that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = GalaxyAgnObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testStars(SearchReversion, StarObj):
    objid = 'testStarDBObj'
    objectTypeId = 91

    #The code below removes the definitions of galacticAv and magNorm
    #from this database object.  The definitions of those columns which
    #are implemented in StarObj rely on mathematical functions which
    #are not defined in sqlite.

    columns = StarObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'galacticAv':
            _to_remove.append(entry)
        elif entry[0] == 'magNorm':
            _to_remove.append(entry)

    for target in _to_remove:
        columns.remove(target)

