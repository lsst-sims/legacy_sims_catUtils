import warnings
import numpy
from .BaseCatalogModels import BaseCatalogObj
from lsst.sims.catalogs.db import ChunkIterator, CatalogDBObject
from lsst.sims.utils import ObservationMetaData

__all__ = ["SolarSystemObj", "CometObj", "NEOObj", "MBAObj", "MiscSolarSystemObj"]

class SolarSystemObj(BaseCatalogObj):
    """
    This is the base CatalogDBObject from which all other Solar System CatalogDBObjects
    derive.  It will query all of the Solar System object tables on fatboy.

    Note: Solar System objects only exist for 50093.14 < mjd < 51923.0.  If you query outside
    of that date range, you will get an empty catalog.
    """
    objid = 'ssm'
    tableid = ''

    ssmtable = 'fSSMAllBase'
    # The variable ssmtable specifies which Table Valued Function to call on
    # fatboy.  Valid values are:
    # fSSMCometBase -- to query comets
    # fSSMNEOBase -- to query Near Earth Objects
    # fSSMMBABase -- to query Main Belt Asteroids
    # fSSMnonMBABase -- to query objects that are none of the above
    # fSSMAllBase -- to query the union of all of the above

    objectTypeId = 40
    doRunTest = True
    testObservationMetaData = ObservationMetaData(boundType = 'circle',
                                                  pointingRA = 0.0, pointingDec = 0.0,
                                                  boundLength = 0.1, mjd=51200., bandpassName='r', m5=22.0)

    #Turn off default column mapping since we are querying a dynamic resource
    generateDefaultColumnMap = False

    #: Numpy can't cast a NoneType to an integer.  This works with floats
    #: as None is cast to nan, but for integers this raises and exception.
    #: Typically it's not an issue as ints are usually ids of some sort,
    #: but in the case of the base galaxy catalog, it's possible for the
    #: varsimobjid to be None if the object does not contain an AGN.
    #: I'm over riding the _postprocess_results method to take care of this.
    #: I could also have refactored my database table so that no integer values
    #: contain NULL values.
    dbDefaultValues = {'varsimobjid':-1, 'myid':-1}

    #: The following maps column names to database schema.  The tuples
    #: must be at least length 2.  If column name is the same as the name
    #: in the DB the mapping element may be None.  The rest of the tuple
    #: should be formatted like a numpy.dtype.  If ommitted, the dtype
    #: is assumed to be float.
    columns = [('objid', 'ssmid', int),
            ('raJ2000', 'ra*PI()/180.'),
            ('decJ2000', 'decl*PI()/180.'),
            ('sedFilename', 'sed_filename', str, 40),
            ('dist', 'dist', numpy.float),
            ('velRa', 'dradt*PI()/180.'),
            ('velDec', 'ddecldt*PI()/180.'),
            ('magNorm', 'magNorm', numpy.float),
            ('umag', 'umag', numpy.float),
            ('gmag', 'gmag', numpy.float),
            ('rmag', 'rmag', numpy.float),
            ('imag', 'imag', numpy.float),
            ('zmag', 'zmag', numpy.float),
            ('ymag', 'ymag', numpy.float)]

    def _get_column_query(self, colnames=None):
        raise NotImplementedError("We are calling a stored procedure so "
                                  "no need to loop over columns")

    def _get_table(self):
        #This needs to be overridden since __init__() calls the default impl.
        return None

    def getIdColKey(self):
        return 'objid'

    def query_columns(self, colnames=None, chunk_size=None, obs_metadata=None, constraint=None):
        """Execute a query

        **Parameters**

            * colnames : list or None
              a list of valid column names, corresponding to entries in the
              `columns` class attribute.  If not specified, all columns are
              queried.
            * chunksize : int (optional)
              if specified, then return an iterator object to query the database,
              each time returning the next `chunksize` elements.  If not
              specified, all matching results will be returned.
            * obs_metadata : object (optional)
              object containing information on the observation including the region of the sky
              to query and time of the observation.
            * constraint : string (optional)
              if specified, the predicate is added to the query verbatim using AND

        **Returns**

            * result : structured array or iterator
              If chunksize is not specified, then result is a structured array of all
              items which match the specified query with columns named by the column
              names in the columns class attribute.  If chunksize is specified,
              then result is an iterator over structured arrays of the given size.

        """
        if colnames is None:
            colnames = [k for k in self.columnMap.keys()]

        mappedcolnames = ["%s as %s"%(self.columnMap[x], x) for x in colnames]
        mappedcolnames = ",".join(mappedcolnames)

        if obs_metadata is not None and obs_metadata.bounds is not None:
            if obs_metadata.bounds.boundType == 'circle':
                regionStr = 'REGION CIRCLE J2000 %f %f %f'%(obs_metadata.bounds.RAdeg,
                                                            obs_metadata.bounds.DECdeg,
                                                            60.*obs_metadata.bounds.radiusdeg)
            elif obs_metadata.bounds.boundType == 'box':
                regionStr = 'REGION RECT J2000 %f %f %f %f'%(obs_metadata.bounds.RAminDeg,
                                                             obs_metadata.bounds.DECminDeg,
                                                             obs_metadata.bounds.RAmaxDeg,
                                                             obs_metadata.bounds.DECmaxDeg)
            else:
                raise RuntimeError("SolarSystemObj does not know about boundType %s "
                                   % obs_metadata.bounds.boundType)
        else:
            regionStr = 'REGION CIRCLE J2000 180. 0. 10800.'
            warnings.warn("Searching over entire sky "
                          "since no bounds specified. "
                          "This could be a very bad idea "
                          "if the database is large")

        query = "select %s from [LSSTCATSIM].[dbo].%s("%(mappedcolnames, self.ssmtable)+\
                "%f, '%s')"%(obs_metadata.mjd.TAI, regionStr)

        if constraint is not None:
            query += "where %s"%(constraint)

        return ChunkIterator(self, query, chunk_size)


class CometObj(SolarSystemObj):
    """
    This CatalogDBObject class queries the table of comets on fatboy

    Note: Solar System objects only exist for 5093.14 < mjd < 51923.0.  If you query outside
    of that date range, you will get an empty catalog.
    """
    ssmtable = 'fSSMCometBase'

    # This object needs its own testObservationMetaData.  The
    # testObservationMetaData defined in SolarSystemObj is too small.
    # It returns zero objects, which causes an error in testAllObjects.py
    testObservationMetaData = ObservationMetaData(boundType = 'circle',
                                                  pointingRA = 0.0, pointingDec = 0.0,
                                                  boundLength = 0.5, mjd=51200., bandpassName='r', m5=22.0)


class NEOObj(SolarSystemObj):
    """
    This CatalogDBObject class queries the table of Near Earth Objects on fatboy

    Note: Solar System objects only exist for 5093.14 < mjd < 51923.0.  If you query outside
    of that date range, you will get an empty catalog.
    """
    ssmtable = 'fSSMNEOBase'

    # This object needs its own testObservationMetaData.  The
    # testObservationMetaData defined in SolarSystemObj is too small.
    # It returns zero objects, which causes an error in testAllObjects.py
    testObservationMetaData = ObservationMetaData(boundType = 'circle',
                                                  pointingRA = 0.0, pointingDec = 0.0,
                                                  boundLength = 0.5, mjd=51200., bandpassName='r', m5=22.0)



class MBAObj(SolarSystemObj):
    """
    This CatalogDBObject class queries the table of Main Belt Asteroids on fatboy

    Note: Solar System objects only exist for 5093.14 < mjd < 51923.0.  If you query outside
    of that date range, you will get an empty catalog.
    """
    ssmtable = 'fSSMMBABase'


class MiscSolarSystemObj(SolarSystemObj):
    """
    This CatalogDBObject class queriess the table of Solar Systems objects on fatboy
    that are not comets, near Earth objects, or main belt asteroids

    Note: Solar System objects only exist for 5093.14 < mjd < 51923.0.  If you query outside
    of that date range, you will get an empty catalog.
    """
    ssmtable = 'fSSMnonMBABase'
