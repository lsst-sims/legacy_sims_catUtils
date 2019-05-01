from sqlalchemy import text
from lsst.utils import getPackageDir
import lsst.pex.config as pexConfig
from lsst.sims.utils import htmModule as htm
from lsst.sims.catalogs.db import CatalogDBObject, ChunkIterator

__all__ = ["LocalStarCatalogObj"]


class LocalStarCatalogConfig(pexConfig.Config):
    host = pexConfig.Field(
        dtype = str,
        doc = "Name of the host",
        default = "localhost",
    )
    port = pexConfig.Field(
        dtype = str,
        doc = "Port number of database",
        default = "1433",
    )
    database = pexConfig.Field(
        dtype = str,
        doc = "Name of database. For 'sqlite', the filename is the database name",
        default = "LSST",
    )
    driver = pexConfig.Field(
        dtype = str,
        doc = "Name of the database backend. Takes format of dialect+driver ",
        default = "mssql+pymssql",
    )


class _HiddenStarCatalogObj(CatalogDBObject):

    objid = 'hiddenstar'
    tableid = None
    idColKey = 'id'
    raColName = 'ra'
    decColName = 'decl'
    objectTypeId = 4

    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 40)]


class LocalStarChunkIterator(ChunkIterator):

    _partition_lim = ((0,11000000000000,8700000000000),
                      (11000000000000, 11600000000000, 11000000000000),
                      (11600000000000, 11800000000000, 11600000000000),
                      (11800000000000, 12000000000000, 11800000000000),
                      (12000000000000, 12200000000000, 12000000000000),
                      (12200000000000, 12500000000000, 12200000000000),
                      (12500000000000, 14000000000000, 12500000000000),
                      (14000000000000, 100000000000000000000, 14000000000000))

    def __init__(self, dbobj, colnames, obs_metadata, chunk_size, constraint,
                 database, host, port, driver):
        """
        Parameters
        ----------
        dbobj -- a CatalogDBObject connected to the 'galaxies' table on fatboy

        colnames -- a list of the columns to query

        chunk_size -- size of chunks to return

        constraint -- a string specifying a SQL 'WHERE' clause

        database -- the name of the database to connect to

        host -- the name of the host to connect to

        port -- the port to connect to

        driver -- the sqlalchemy driver to use
        """
        self.database = database
        self.host = host
        self.port = port
        self.driver = driver

        self.arbitrarySQL = False
        self._chunk_size = chunk_size
        self._obs_metadata = obs_metadata
        self._constraint = constraint
        self._colnames = colnames

        half_space = htm.halfSpaceFromRaDec(obs_metadata.pointingRA,
                                            obs_metadata.pointingDec,
                                            obs_metadata.boundLength)

        self._trixel_search_level = 12
        self._trixel_bounds = half_space.findAllTrixels(self._trixel_search_level)

        self._tables_to_query = set()
        self._tables_to_query.add('starsRRLy')
        where_clause = '('
        global_min_21 = None
        global_max_21 = None
        for bound in self._trixel_bounds:
            min_21 = bound[0] << 2*(21-self._trixel_search_level)
            max_21 = (bound[1]+1) << 2*(21-self._trixel_search_level)
            if global_min_21 is None or min_21<global_min_21:
                global_min_21 = min_21
            if global_max_21 is None or max_21>global_max_21:
                global_max_21 = max_21

            if where_clause != '(':
                where_clause += ' OR '

            if min_21 == max_21:
                where_clause += 'htmid==%d' % min_21
            else:
                where_clause += '(htmid>=%d AND htmid<=%d)' % (min_21, max_21)

            for part in self._partition_lim:
                part_name = 'stars_partition_%d' % part[2]
                if min_21>=part[0] and min_21<part[1]:
                    self._tables_to_query.add(part_name)
                elif max_21>=part[0] and max_21<part[1]:
                    self._tables_to_query.add(part_name)
                elif min_21<=part[0] and max_21>part[1]:
                    self._tables_to_query.add(part_name)

        where_clause += ')'
        self._htmid_where_clause = '(htmid>=%d AND htmid<=%d AND ' % (global_min_21, global_max_21)
        self._htmid_where_clause += where_clause
        self._htmid_where_clause += ')'

        self._active_query = None

    def _load_next_star_db(self, colnames):
        table_name = self._tables_to_query.pop()
        db = _HiddenStarCatalogObj(table=table_name,
                                   database=self.database,
                                   host=self.host,
                                   port=self.port,
                                   driver=self.driver)
        self.dbobj = db
        column_query = db._get_column_query(colnames)
        column_query = column_query.filter(text(self._htmid_where_clause))
        if self._constraint is not None:
            column_query = column_query.filter(text(self._constraint))
        exec_query = db.connection.session.execute(column_query)
        return exec_query

    def __next__(self):

        if self._active_query is None or self._active_query.closed:
            if len(self._tables_to_query) == 0:
                raise StopIteration
            self._active_query = self._load_next_star_db(self._colnames)

        if self._chunk_size is None:
            chunk = self._active_query.fetchall()
        elif self._chunk_size is not None:
            chunk = self._active_query.fetchmany(self._chunk_size)
            if len(chunk) == 0:
                self._active_query.close()
                return self.__next__()

        return self._postprocess_results(chunk)


class LocalStarCatalogObj(CatalogDBObject):

    config = LocalStarCatalogConfig()

    database = config.database
    host = config.host
    port = config.port
    driver = config.driver

    idColKey = 'id'
    objid = 'epycStarBase'
    tableid = 'stars_partition_8700000000000'  # just a placeholder

    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 40)]

    def get_table(self):
        """
        We don't actually want this to do anything, since this CatalogDBObject
        does its own search over all of the available star tables
        """
        pass

    def query_columns(self, colnames=None, chunk_size=None,
                      obs_metadata=None, constraint=None,
                      limit=None):

        if obs_metadata.boundType != 'circle':
            raise RuntimeError("Cannot use boundType %s in this catalog; only 'circle'"
                               % str(obs_metadata.boundType))

        return LocalStarChunkIterator(self, colnames, obs_metadata, chunk_size,
                                      constraint,
                                      self.database,
                                      self.host,
                                      self.port,
                                      self.driver)
