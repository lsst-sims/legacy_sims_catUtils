import warnings
import numpy
from .BaseCatalogModels import BaseCatalogObj
from lsst.sims.catalogs.db import ChunkIterator, CompoundCatalogDBObject
from lsst.sims.utils import ObservationMetaData

__all__ = ["GalaxyObj", "GalaxyTileObj", "GalaxyBulgeObj",
           "GalaxyDiskObj", "GalaxyAgnObj", "ImageAgnObj", "LensGalaxyObj",
           "GalaxyTileCompoundObj"]


class GalaxyObj(BaseCatalogObj):
    """
    Note: building a catalog out of this object will directly call the
    'galaxy' table.  This table only contains objects for

    -2.5 deg < RA < 2.5 deg, -2.5 deg < Dec < 2.5 deg

    In order to cover the whole sky, call one of the objects that
    inherits from GalaxyTileObj
    """
    objid = 'galaxyBase'

    #: This is the base table for the galaxies
    # tableid = 'final_clone_db'
    tableid = 'galaxy'
    idColKey = 'id'
    raColName = '((CAST(ra AS NUMERIC(9,6))%360.)+360.)%360.'
    decColName = 'dec'
    objectTypeId = 24

    doRunTest = True
    testObservationMetaData = ObservationMetaData(boundType='circle', pointingRA=0.0, pointingDec=0.0,
                                                  boundLength=0.01, mjd=52000., bandpassName='r', m5 = 22.0)

    #: Numpy can't cast a NoneType to an integer.  This works with floats
    #: as None is cast to nan, but for integers this raises and exception.
    #: Typically it's not an issue as ints are usually ids of some sort,
    #: but in the case of the base galaxy catalog, it's possible for the
    #: varsimobjid to be None if the object does not contain an AGN.
    #: I'm over riding the _postprocess_results method to take care of this.
    #: I could also have refactored my database table so that no integer values
    #: contain NULL values.
    dbDefaultValues = {'varsimobjid': -1, 'myid': -1}

    #: The following maps column names to database schema.  The tuples
    #: must be at least length 2.  If column name is the same as the name
    #: in the DB the mapping element may be None.  The rest of the tuple
    #: should be formatted like a numpy.dtype.  If ommitted, the dtype
    #: is assumed to be float.
    columns = [('galid', None, str, 30),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'dec*PI()/180.'),
               ('raJ2000Bulge', 'bra*PI()/180.'),
               ('decJ2000Bulge', 'bdec*PI()/180.'),
               ('raJ2000Disk', 'dra*PI()/180.'),
               ('decJ2000Disk', 'ddec*PI()/180.'),
               ('raJ2000Agn', 'agnra*PI()/180.'),
               ('decJ2000Agn', 'agndec*PI()/180.'),
               ('magNormBulge', 'magnorm_bulge'),
               ('magNormDisk', 'magnorm_disk'),
               ('magNormAgn', 'magnorm_agn'),
               ('sedFilenameBulge', 'sedname_bulge', str, 40),
               ('sedFilenameDisk', 'sedname_disk', str, 40),
               ('sedFilenameAgn', 'sedname_agn', str, 40),
               ('majorAxisBulge', 'a_b*PI()/648000.'),
               ('minorAxisBulge', 'b_b*PI()/648000.'),
               ('positionAngleBulge', 'pa_bulge*PI()/180.'),
               ('sindexBulge', 'bulge_n', int),
               ('majorAxisDisk', 'a_d*PI()/648000.'),
               ('minorAxisDisk', 'b_d*PI()/648000.'),
               ('positionAngleDisk', 'pa_disk*PI()/180.'),
               ('sindexDisk', 'disk_n', int),
               ('internalExtinctionModelBulge', 'ext_model_b', str, 3),
               ('internalAvBulge', 'av_b'),
               ('internalRvBulge', 'rv_b'),
               ('internalExtinctionModelDisk', 'ext_model_d', str, 3),
               ('internalAvDisk', 'av_d'),
               ('internalRvDisk', 'rv_d'),
               ('lsst_u', 'u_ab'),
               ('lsst_g', 'g_ab'),
               ('lsst_r', 'r_ab'),
               ('lsst_i', 'i_ab'),
               ('lsst_z', 'z_ab'),
               ('lsst_y', 'y_ab')]

    def _final_pass(self, results):
        """This is to map the values from 0 - 2*PI() as ra goes negative currently"""
        for ra in ('raJ2000', 'raJ2000Bulge', 'raJ2000Disk', 'raJ2000Agn'):
            if ra in results.dtype.names:
                results[ra] = results[ra]%(numpy.pi*2.)
        return results


class GalaxyTileObj(BaseCatalogObj):
    """
    This is the parent class for galaxy BaseCatalogObjs that sample the whole
    sky (rather than just a very small patch as in GalaxyObj)
    """

    objid = 'galaxyTiled'
    #: This is the base table for the galaxies
    tableid = 'galaxy'
    # componentSubset must be one of "ALL" (default), "AGN", "BULGE", "DISK"
    componentSubset = 'ALL'
    raColName = 'ra'
    decColName = 'dec'
    objectTypeId = 25

    doRunTest = True
    testObservationMetaData = ObservationMetaData(boundType='circle', pointingRA=173.0, pointingDec=-60.0,
                                                  boundLength=0.01, mjd=52000., bandpassName='r', m5=22.0)

    #: Numpy can't cast a NoneType to an integer.  This works with floats
    #: as None is cast to nan, but for integers this raises and exception.
    #: Typically it's not an issue as ints are usually ids of some sort,
    #: but in the case of the base galaxy catalog, it's possible for the
    #: varsimobjid to be None if the object does not contain an AGN.
    #: I'm over riding the _postprocess_results method to take care of this.
    #: I could also have refactored my database table so that no integer values
    #: contain NULL values.
    dbDefaultValues = {'varsimobjid': -1, 'myid': -1}

    #: The following maps column names to database schema.  The tuples
    #: must be at least length 2.  If column name is the same as the name
    #: in the DB the mapping element may be None.  The rest of the tuple
    #: should be formatted like a numpy.dtype.  If ommitted, the dtype
    #: is assumed to be float.
    columns = [('galtileid', None, numpy.int64),
               ('galid', None, str, 30),
               ('raJ2000', 'ra'),
               ('decJ2000', 'dec'),
               ('raJ2000Bulge', 'bra*PI()/180.'),
               ('decJ2000Bulge', 'bdec*PI()/180.'),
               ('raJ2000Disk', 'dra*PI()/180.'),
               ('decJ2000Disk', 'ddec*PI()/180.'),
               ('raJ2000Agn', 'agnra*PI()/180.'),
               ('decJ2000Agn', 'agndec*PI()/180.'),
               ('magNormBulge', 'magnorm_bulge'),
               ('magNormDisk', 'magnorm_disk'),
               ('magNormAgn', 'magnorm_agn'),
               ('sedFilenameBulge', 'sedname_bulge', str, 40),
               ('sedFilenameDisk', 'sedname_disk', str, 40),
               ('sedFilenameAgn', 'sedname_agn', str, 40),
               ('majorAxisBulge', 'a_b*PI()/648000.'),
               ('minorAxisBulge', 'b_b*PI()/648000.'),
               ('positionAngleBulge', 'pa_bulge*PI()/180.'),
               ('sindexBulge', 'bulge_n', int),
               ('majorAxisDisk', 'a_d*PI()/648000.'),
               ('minorAxisDisk', 'b_d*PI()/648000.'),
               ('positionAngleDisk', 'pa_disk*PI()/180.'),
               ('sindexDisk', 'disk_n', int),
               ('internalExtinctionModelBulge', 'ext_model_b', str, 3),
               ('internalAvBulge', 'av_b'),
               ('internalRvBulge', 'rv_b'),
               ('internalExtinctionModelDisk', 'ext_model_d', str, 3),
               ('internalAvDisk', 'av_d'),
               ('internalRvDisk', 'rv_d'),
               ('lsst_u', 'u_ab'),
               ('lsst_g', 'g_ab'),
               ('lsst_r', 'r_ab'),
               ('lsst_i', 'i_ab'),
               ('lsst_z', 'z_ab'),
               ('lsst_y', 'y_ab')]

    def _get_column_query(self, colnames=None):
        raise NotImplementedError("We are calling a stored procedure so "
                                  "no need to loop over columns")

    def _final_pass(self, results):
        """Modify the results of raJ2000 and decJ2000 to be in radians
        **Parameters**

            * results : Structured array of results from query

        **Returns**

            * results : Modified structured array

        """

        results['raJ2000'] = numpy.radians(results['raJ2000'])
        results['decJ2000'] = numpy.radians(results['decJ2000'])
        return results

    def getIdColKey(self):
        return 'galtileid'

    def query_columns(self, colnames=None, chunk_size=None, obs_metadata=None, constraint=None,
                      limit=None):
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
            * limit:
              This kwarg is not actually used.  It exists to preserve the same interface
              as other definitions of query_columns elsewhere in CatSim.  If not None,
              a warning will be emitted, pointing out to the user that 'limit' is not used.

        **Returns**

            * result : structured array or iterator
              If chunksize is not specified, then result is a structured array of all
              items which match the specified query with columns named by the column
              names in the columns class attribute.  If chunksize is specified,
              then result is an iterator over structured arrays of the given size.

        """
        if colnames is None:
            colnames = [k for k in self.columnMap.keys()]

        # We know that galtileid comes back with the query, but we don't want
        # to add it to the query since it's generated on the fly.
        #
        # 25 August 2015
        # The code below has been modified to remove all column names
        # that contain 'galtileid.'  This is to accommodate the
        # CompoundInstanceCatalog and CompoundDBObject classes, which
        # mangle column names such that they include the objid of the
        # specific CatalogDBObject that is asking for them.
        for name in colnames:
            if 'galtileid' in name:
                colnames.remove(name)

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
                raise RuntimeError("GalaxyTileObj does not know about boundType %s "
                                   % obs_metadata.bounds.boundType)
        else:
            regionStr = 'REGION CIRCLE J2000 180. 0. 10800.'
            warnings.warn("Searching over entire sky "
                          "since no bounds specified. "
                          "This could be a very bad idea "
                          "if the database is large")

        query = """EXECUTE [LSSTCATSIM].[dbo].[GalaxySearch2015]
                   @ApertureStr = '%s', @ColumnNames = '%s',
                   @ComponentSubset = '%s' """ % (regionStr, mappedcolnames, self.componentSubset)

        if constraint is not None:
            query += ", @WhereClause = '%s'"%(constraint)

        if limit is not None:
            warnings.warn("You specified a row number limit in your query of a GalaxyTileObj "
                          "daughter class.  Because of the way GalaxyTileObj is searched, row "
                          "number limits are not possible.  If you really want to limit the number "
                          "of rows returned by you query, consider using GalaxyObj (note that you "
                          "will have to you limit your search to -2.5<RA<2.5 -2.25<Dec<2.25 -- both in "
                          "degrees -- as this is the only region where galaxies exist in GalaxyObj).")

        return ChunkIterator(self, query, chunk_size)


class GalaxyTileCompoundObj(CompoundCatalogDBObject, GalaxyTileObj):
    """
    This is a daughter class of CompoundCatalogDBObject that specifically
    inherits from GalaxyTileObj.  This is necessary because GalaxyTileObj
    implements a custom query_columns that executes a stored procedure
    on fatboy (as opposed to the generic query_columns implemented in
    CatalogDBObject, which converts self.columns into a SQL
    query in the most naive way possible).
    """

    doRunTest = False  # because this is not like other CatalogDBObjects
                       # and should not be tested in the same way

    _table_restriction = ['galaxy']

    def _final_pass(self, results):
        for name in results.dtype.fields:
            if 'raJ2000' in name or 'decJ2000' in name:
                results[name] = numpy.radians(results[name])

        return results


class GalaxyBulgeObj(GalaxyTileObj):
    objid = 'galaxyBulge'
    componentSubset = 'BULGE'
    raColName = 'ra'
    decColName = 'dec'
    objectTypeId = 26
    doRunTest = True
    testObservationMetaData = ObservationMetaData(boundType='circle', pointingRA=10.0, pointingDec=-45.0,
                                                  boundLength=0.01, mjd=53000., bandpassName='r', m5=22.0)
    #: The following maps column names to database schema.  The tuples
    #: must be at least length 2.  If column name is the same as the name
    #: in the DB the mapping element may be None.  The rest of the tuple
    #: should be formatted like a numpy.dtype.  If ommitted, the dtype
    #: is assumed to be float.
    columns = [('galtileid', None, numpy.int64),
               ('galid', None, str, 30),
               ('componentra', 'bra*PI()/180.'),
               ('componentdec', 'bdec*PI()/180.'),
               #: This is actually a problem with the stored procedure.
               #: We need to be able to map columns other than
               #: just ra/dec to raJ2000/decJ2000.  This gets
               #: important when we start perturbing the three galaxy components
               ('raJ2000', 'ra'),
               ('decJ2000', 'dec'),
               ('magNorm', 'magnorm_bulge'),
               ('magNormBulge', 'magnorm_bulge'),
               ('sedFilename', 'sedname_bulge', str, 40),
               ('sedFilenameBulge', 'sedname_bulge', str, 40),
               ('majorAxis', 'a_b*PI()/648000.'),
               ('minorAxis', 'b_b*PI()/648000.'),
               ('positionAngle', 'pa_bulge*PI()/180.'),
               ('sindex', 'bulge_n', int),
               ('halfLightRadius', 'BulgeHalfLightRadius*PI()/648000.'),
               ('internalExtinctionModel', 'ext_model_b', str, 3),
               ('internalAv', 'av_b'),
               ('internalRv', 'rv_b'),
               ('internalAvBulge', 'av_b'),
               ('internalRvBulge', 'rv_b'),
               ('lsst_u', 'u_ab'),
               ('lsst_g', 'g_ab'),
               ('lsst_r', 'r_ab'),
               ('lsst_i', 'i_ab'),
               ('lsst_z', 'z_ab'),
               ('lsst_y', 'y_ab')]


class GalaxyDiskObj(GalaxyTileObj):
    objid = 'galaxyDisk'
    componentSubset = 'DISK'
    raColName = 'ra'
    decColName = 'dec'
    objectTypeId = 27
    doRunTest = True
    testObservationMetaData = ObservationMetaData(boundType='circle', pointingRA=66.0, pointingDec=-80.0,
                                                  boundLength=0.01, mjd=53730., bandpassName='r', m5=22.0)
    #: The following maps column names to database schema.  The tuples
    #: must be at least length 2.  If column name is the same as the name
    #: in the DB the mapping element may be None.  The rest of the tuple
    #: should be formatted like a numpy.dtype.  If ommitted, the dtype
    #: is assumed to be float.
    columns = [('galtileid', None, numpy.int64),
               ('galid', None, str, 30),
               ('componentra', 'dra*PI()/180.'),
               ('componentdec', 'ddec*PI()/180.'),
               #: This is actually a problem with the stored procedure.
               #: We need to be able to map columns other than
               #: just ra/dec to raJ2000/decJ2000.  This gets
               #: important when we start perturbing the three galaxy components
               ('raJ2000', 'ra'),
               ('decJ2000', 'dec'),
               ('magNorm', 'magnorm_disk'),
               ('magNormDisk', 'magnorm_disk'),
               ('sedFilename', 'sedname_disk', str, 40),
               ('sedFilenameDisk', 'sedname_disk', str, 40),
               ('majorAxis', 'a_d*PI()/648000.'),
               ('minorAxis', 'b_d*PI()/648000.'),
               ('positionAngle', 'pa_disk*PI()/180.'),
               ('sindex', 'disk_n', int),
               ('halfLightRadius', 'DiskHalfLightRadius*PI()/648000.'),
               ('internalExtinctionModel', 'ext_model_d', str, 3),
               ('internalAv', 'av_d'),
               ('internalRv', 'rv_d'),
               ('internalAvDisk', 'av_d'),
               ('internalRvDisk', 'rv_d'),
               ('lsst_u', 'u_ab'),
               ('lsst_g', 'g_ab'),
               ('lsst_r', 'r_ab'),
               ('lsst_i', 'i_ab'),
               ('lsst_z', 'z_ab'),
               ('lsst_y', 'y_ab')]


class GalaxyAgnObj(GalaxyTileObj):
    objid = 'galaxyAgn'
    componentSubset = 'AGN'
    raColName = 'ra'
    decColName = 'dec'
    objectTypeId = 28
    doRunTest = True
    testObservationMetaData = ObservationMetaData(boundType='circle', pointingRA=234.0, pointingDec=-15.0,
                                                  boundLength=0.01, mjd=51000., bandpassName='r', m5=22.0)
    #: The following maps column names to database schema.  The tuples
    #: must be at least length 2.  If column name is the same as the name
    #: in the DB the mapping element may be None.  The rest of the tuple
    #: should be formatted like a numpy.dtype.  If ommitted, the dtype
    #: is assumed to be float.
    columns = [('galtileid', None, numpy.int64),
               ('galid', None, str, 30),
               ('componentra', 'agnra*PI()/180.'),
               ('componentdec', 'agndec*PI()/180.'),
               #: This is actually a problem with the stored procedure.
               #: We need to be able to map columns other than
               #: just ra/dec to raJ2000/decJ2000.  This gets
               #: important when we start perturbing the three galaxy components
               ('raJ2000', 'ra'),
               ('decJ2000', 'dec'),
               ('magNorm', 'magnorm_agn'),
               ('magNormAgn', 'magnorm_agn'),
               ('sedFilename', 'sedname_agn', str, 40),
               ('sedFilenameAgn', 'sedname_agn', str, 40),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('lsst_u', 'u_ab'),
               ('lsst_g', 'g_ab'),
               ('lsst_r', 'r_ab'),
               ('lsst_i', 'i_ab'),
               ('lsst_z', 'z_ab'),
               ('lsst_y', 'y_ab')]


class ImageAgnObj(BaseCatalogObj):
    objid = 'imageagn'
    tableid = 'image'
    idColKey = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    objectTypeId = 29
    doRunTest = True
    #: all sky since this is a small set.
    testObservationMetaData = ObservationMetaData(mjd=53000., bandpassName='r', m5=22.0)

    dbDefaultValues = {'varsimobjid': -1, 'myid': -1}
    #: The following maps column names to database schema.  The tuples
    #: must be at least length 2.  If column name is the same as the name
    #: in the DB the mapping element may be None.  The rest of the tuple
    #: should be formatted like a numpy.dtype.  If ommitted, the dtype
    #: is assumed to be float.
    columns = [('galid', 'id', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'dec*PI()/180.'),
               ('magNorm', 'magnorm_agn'),
               ('sedFilename', 'sedname_agn', str, 40),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('lsst_u', 'u_ab'),
               ('lsst_g', 'g_ab'),
               ('lsst_r', 'r_ab'),
               ('lsst_i', 'i_ab'),
               ('lsst_z', 'z_ab'),
               ('lsst_y', 'y_ab')]


class LensGalaxyObj(BaseCatalogObj):
    objid = 'lensgalaxy'
    tableid = 'lens'
    idColKey = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    objectTypeId = 30
    doRunTest = True
    #: all sky since this is a small set.
    testObservationMetaData = ObservationMetaData(mjd=53000., bandpassName='r', m5=22.0)

    dbDefaultValues = {'varsimobjid': -1, 'myid': -1, 'variabilityParameters': None}
    #: The following maps column names to database schema.  The tuples
    #: must be at least length 2.  If column name is the same as the name
    #: in the DB the mapping element may be None.  The rest of the tuple
    #: should be formatted like a numpy.dtype.  If ommitted, the dtype
    #: is assumed to be float.
    columns = [('galid', 'id', int),
               ('raJ2000', 'ra_bulge*PI()/180.'),
               ('decJ2000', 'dec_bulge*PI()/180.'),
               ('magNorm', 'magnorm_bulge'),
               ('sedFilename', 'sedname_bulge', str, 40),
               ('lsst_u', 'u_ab'),
               ('lsst_g', 'g_ab'),
               ('lsst_r', 'r_ab'),
               ('lsst_i', 'i_ab'),
               ('lsst_z', 'z_ab'),
               ('lsst_y', 'y_ab')]
