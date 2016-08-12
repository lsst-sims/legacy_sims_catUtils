import warnings
import os
from lsst.sims.catalogs.db import CatalogDBObject,  ChunkIterator
from sqlalchemy.sql import select, func, column
import lsst.pex.config as pexConfig
from lsst.utils import getPackageDir

__all__ = ["BaseCatalogObj", "BaseCatalogConfig"]

class BaseCatalogConfig(pexConfig.Config):
    host = pexConfig.Field(
        dtype = str,
        doc = "Name of the host",
        default = "fatboy-private.phys.washington.edu",
    )
    port = pexConfig.Field(
        dtype = str,
        doc = "Port number of database",
        default = "1433",
    )
    database = pexConfig.Field(
        dtype = str,
        doc = "Name of database. For 'sqlite', the filename is the database name",
        default = "LSSTCATSIM",
    )
    driver = pexConfig.Field(
        dtype = str,
        doc = "Name of the database backend. Takes format of dialect+driver ",
        default = "mssql+pymssql",
    )

class BaseCatalogObj(CatalogDBObject):
    """Base class for Catalogs that query the default
    UW CATSIM database
    """

    config = BaseCatalogConfig()

    #load $SIMS_CATUTILS_DIR/config/db.py
    config.load(os.path.join(getPackageDir("sims_catUtils"), "config", "db.py"))

    host = config.host
    port = config.port
    database = config.database
    driver = config.driver

    def query_columns(self, colnames=None, chunk_size=None,
                      obs_metadata=None, constraint=None,
                      limit=None):
        """Execute a query from the primary catsim database

        Execute a query, taking advantage of the spherical geometry library and
        htmid indexes on all catalog tables in the UW catsim database

        **Parameters**

            * colnames : list or None
              a list of valid column names, corresponding to entries in the
              `columns` class attribute.  If not specified, all columns are
              queried.
            * chunk_size : int (optional)
              if specified, then return an iterator object to query the database,
              each time returning the next `chunk_size` elements.  If not
              specified, all matching results will be returned.
            * obs_metadata : object (optional)
              an observation metadata object which has a "filter" method, which
              will add a filter string to the query.
            * constraint : str (optional)
              a string which is interpreted as SQL and used as a predicate on the query
            * limit : int (optional)
              limits the number of rows returned by the query

        **Returns**

            * result : list or iterator
              If chunk_size is not specified, then result is a list of all
              items which match the specified query.  If chunk_size is specified,
              then result is an iterator over lists of the given size.
        """
        query = self._get_column_query(colnames)

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
                raise RuntimeError("CatalogObject does not know about boundType %s "
                                   % obs_metadata.bounds.boundType)
        else:
            regionStr = 'REGION CIRCLE J2000 180. 0. 10800.'
            warnings.warn("Searching over entire sky "
                          "since no bounds specified. "
                          "This could be a very bad idea "
                          "if the database is large")

        if obs_metadata is not None and regionStr is not None:
            #add spatial constraints to query.

            #Hint sql engine to seek on htmid
            if not self.tableid.endswith('forceseek'):
                query = query.with_hint(self.table, ' WITH(FORCESEEK)', 'mssql')

            #aliased subquery for htmid ranges covering the search region
            htmid_range_alias = select([column('htmidstart'), column('htmidend')]).\
            select_from(func.fHtmCoverRegion(regionStr)).alias()

            #SQL is not case sensitive but python is:
            if 'htmID' in self.columnMap:
                htmidName = 'htmID'
            elif 'htmid' in self.columnMap:
                htmidName = 'htmid'
            else:
                htmidName = 'htmId'

            #Range join on htmid ranges
            query = query.join(htmid_range_alias,
                       self.table.c[htmidName].between(htmid_range_alias.c.htmidstart,
                                                     htmid_range_alias.c.htmidend)
                       )
            query = query.filter(func.sph.fRegionContainsXYZ(func.sph.fSimplifyString(regionStr),
                     self.table.c.cx, self.table.c.cy, self.table.c.cz) == 1)


        if constraint is not None:
            query = query.filter(constraint)

        if limit is not None:
            query = query.limit(limit)

        return ChunkIterator(self, query, chunk_size)

