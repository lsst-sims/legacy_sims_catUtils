import os
import copy
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.baseCatalogModels import StarObj, GalaxyTileObj, GalaxyAgnObj, \
                                                 GalaxyDiskObj, GalaxyBulgeObj

__all__ = ["SearchReversion", "testGalaxyTileDBObj",
           "testGalaxyBulgeDBObj", "testGalaxyDiskDBObj", "testGalaxyAgnDBObj",
           "testStarsDBObj"]

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

class testGalaxyTileDBObj(SearchReversion, GalaxyTileObj):
    objid = 'testGalaxyDBObj'
    objectTypeId = 87

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject class assumes that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = copy.deepcopy(GalaxyTileObj.columns)
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testGalaxyBulgeDBObj(SearchReversion, GalaxyBulgeObj):
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

    columns = copy.deepcopy(GalaxyBulgeObj.columns)
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testGalaxyDiskDBObj(SearchReversion, GalaxyDiskObj):
    objid = 'testDiskDBObj'
    objectTypeId = 89

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject daughter classes assume that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = copy.deepcopy(GalaxyDiskObj.columns)
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testGalaxyAgnDBObj(SearchReversion, GalaxyAgnObj):
    objid = 'testAgnDBObj'
    objectTypeId = 90

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject daughter classes assume that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = copy.deepcopy(GalaxyAgnObj.columns)
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testStarsDBObj(SearchReversion, StarObj):
    objid = 'testStarDBObj'
    objectTypeId = 91

    #The code below removes the definitions of galacticAv and magNorm
    #from this database object.  The definitions of those columns which
    #are implemented in StarObj rely on mathematical functions which
    #are not defined in sqlite.

    columns = copy.deepcopy(StarObj.columns)
    _to_remove = []
    for entry in columns:
        if entry[0] == 'galacticAv':
            _to_remove.append(entry)
        elif entry[0] == 'magNorm':
            _to_remove.append(entry)

    for target in _to_remove:
        columns.remove(target)

