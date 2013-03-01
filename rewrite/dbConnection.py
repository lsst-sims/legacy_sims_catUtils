import warnings
import math

from sqlalchemy.orm import scoped_session, sessionmaker, mapper
from sqlalchemy.sql import expression
from sqlalchemy import (create_engine, ThreadLocalMetaData, MetaData,
                        Table, Column, BigInteger)


DEFAULT_ADDRESS = "mssql+pymssql://LSST-2:L$$TUser@fatboy.npl.washington.edu:1433/LSST"


class ObservationMetaData(object):
    """Observation Metadata
    
    This class contains any metadata for a query which is associated with
    a particular telescope pointing, including bounds in RA and DEC, and
    the time of the observation.

    Parameters
    ----------
    circ_bounds : dict (optional)
        a dictionary with the keys 'ra', 'dec', and 'radius' measured, in
        degrees
    box_bounds : dict (optional)
        a dictionary with the keys 'ra_min', 'ra_max', 'dec_min', 'dec_max',
        measured in degrees
    MJD : list (optional)
        a list of MJD values to query

    Examples
    --------
    >>> data = ObservationMetaData.from_obshistid(88544919)
    >>> print data.MJD
    """
    @classmethod
    def from_obshistid(cls, obshistid):
        raise NotImplementedError()
    
    def __init__(self, circ_bounds=None, box_bounds=None, mjd=None):
        self.circ_bounds = circ_bounds
        self.box_bounds = box_bounds
        self.mjd = mjd

    def filter(self, query, RAname='ra', DECname='decl', MJDname='mjd'):
        """Filter the query by the associated metadata"""
        on_clause = self.to_SQL(RAname, DECname, MJDname)
        if on_clause:
            query = query.filter(on_clause)
        return query

    def to_SQL(self, RAname='ra', DECname='decl', MJDname='mjd'):
        constraint = ""
        if self.box_bounds is not None:
            bb = self.box_bounds
            constraint += self.box_bound_constraint(bb['ra_min'],
                                                    bb['ra_max'],
                                                    bb['dec_min'],
                                                    bb['dec_max'],
                                                    RAname, DECname)
        if self.circ_bounds is not None:
            cb = self.circ_bounds
            constraint += self.circle_bound_constraint(cb['ra'], cb['dec'],
                                                       cb['radius'],
                                                       RAname, DECname)
        if self.mjd is not None:
            constraint += self.mjd_constraint(self.mjd, MJDname)
            
        return constraint

    @staticmethod
    def mdj_constraint(MJD, MJDname):
        raise NotImplementedError("haven't implemented MJD bound yet")

    @staticmethod
    def box_bound_constraint(RAmin, RAmax, DECmin, DECmax,
                             RAname='ra', DECname='decl'):
        (RAmin, RAmax, DECmin, DECmax) = map(math.radians,
                                             (RAmin, RAmax, DECmin, DECmax))

        if RAmin < 0 and RAmax > 2. * math.pi:
            bound = "%s between %f and %f" % (DECname, DECmin, DECmax)

        elif RAmin < 0 and RAmax <= 2 * math.pi:
            # XXX is this right?  It seems strange.
            bound = ("%s not between %f and %f and %s between %f and %f"
                     % (RAname, RAmin % (2 * math.pi), RAmax,
                        DECname, DECmin, DECmax))

        elif RAmin >= 0 and RAmax > 2. * math.pi:
            bound = ("%s not between %f and %f and %s between %f and %f" 
                     % (RAname, RAmin, RAmax % (2 * math.pi),
                        DECname, DECmin, DECmax))

        else:
            bound = ("%s between %f and %f and %s between %f and %f"
                     % (RAname, RAmin, RAmax, DECname, DECmin, DECmax))

        return bound

    @staticmethod
    def circle_bound_constraint(RA, DEC, radius,
                                RAname='ra', DECname='decl'):
        RAmax = RA + radius / math.cos(math.radians(DEC))
        RAmin = RA - radius / math.cos(math.radians(DEC))
        DECmax = DEC + radius
        DECmin = DEC - radius
        return ObservationMetaData.box_bound_constraint(RAmin, RAmax,
                                                        DECmin, DECmax,
                                                        RAname, DECname)    


#------------------------------------------------------------
# Iterator for database chunks
class ChunkIterator(object):
    """Iterator for query chunks"""
    def __init__(self, exec_query, chunksize):
        self.dbobj = dbobj
        self.exec_query = dbobj.session.execute(query)
        self.chunksize = chunksize

    def __iter__(self):
        return self

    def next(self):
        chunk = self.exec_query.fetchmany(self.chunksize)
        if len(chunk) == 0:
            raise StopIteration
        return self.dbobj._postprocess_results(chunk)


class DBObjectMeta(type):
    """Meta class for registering new objects.

    When any new type of object class is created, this registers it
    in a `registry` class attribute, available to all derived instance
    catalog.
    """
    def __new__(cls, name, bases, dct):
        # check if attribute objid is specified.
        # If not, create a default
        if 'registry' in dct:
            warnings.warn("registry class attribute should not be "
                          "over-ridden in InstanceCatalog classes. "
                          "Proceed with caution")
        if 'objid' not in dct:
            dct['objid'] = name
        return super(DBObjectMeta, cls).__new__(cls, name, bases, dct)

    def __init__(cls, name, bases, dct):
        # check if 'registry' is specified.
        # if not, then this is the base class: add the registry
        if not hasattr(cls, 'registry'):
            cls.registry = {}
        else:
            # add this class to the registry
            if cls.objid in cls.registry:
                warnings.warn('duplicate object id %s specified' % cls.objid)
            cls.registry[cls.objid] = cls

            if not hasattr(cls, 'columns') or cls.columns is None:
                cls.columns = cls.column_map.keys()
                cls.requirements = cls.column_map
            else:
                # build requirements dict from columns and column_map
                cls.requirements = dict([(key, cls.column_map.get(key, key))
                                         for key in cls.columns])
            
        return super(DBObjectMeta, cls).__init__(name, bases, dct)


class DBObject(object):
    """Database Object base class

    """
    __metaclass__ = DBObjectMeta
    objid = None
    tableid = None
    idColName = None
    appendint = None
    spatialModel = None
    columns = []
    column_map = {}

    @classmethod
    def from_objid(cls, objid, *args, **kwargs):
        """Given a string objid, return an instance of
        the appropriate DBObject class.

        objid should map to an entry in the objectMap.dat configuration
        file.  If objid does not match any subclass of DBObjectBase,
        then a generic DBObject will be returned.
        """
        cls = cls.registry.get(objid, DBObject)
        return cls(*args, **kwargs)

    def __init__(self, address=None):
        if (self.objid is None) or (self.tableid is None)\
           or (self.appendint is None) or (self.spatialModel is None):
            raise ValueError("DBObject must be subclassed, and "
                             "define objid, tableid, appendint, and spatialModel.")

        if address is None:
            self.address = DEFAULT_ADDRESS
        else:
            self.address = address

        self._connect_to_engine()
        self._get_table()

    def getObjectTypeId(self):
        return self.appendint

    def getSpatialModel(self):
        return self.spatialModel

    def _get_table(self):
        self.table = Table(self.tableid, self.metadata,
                           Column(self.idColName, BigInteger, primary_key=True),
                           autoload=True)

    def _connect_to_engine(self):
        """create and connect to a database engine"""
        self.engine = create_engine(self.address, echo=False)
        self.session = scoped_session(sessionmaker(autoflush=True, 
                                                   bind=self.engine))
        self.metadata = MetaData()
        self.metadata.bind = self.engine

    def _get_column_query(self, colnames=None):
        """Given a list of valid column names, return the query object"""
        if colnames is None:
            colnames = self.columns

        try:
            vals = [self.requirements[col] for col in colnames]
        except KeyError:
            raise ValueError('entries in colnames must be in self.columns')

        # Get the first query
        if self.idColName in vals:
            idLabel = colnames[vals.index(self.idColName)]
        else:
            idLabel = 'id'

        query = self.session.query(self.table.c[self.idColName].label(idLabel))

        for col, val in zip(colnames[1:], vals[1:]):
            if val == self.idColName:
                pass
            query = query.add_column(expression.literal_column(val).label(col))

        return query

    def _postprocess_results(self, results):
        """Post-process the query results.
        This can be overridden in derived classes: by default it returns
        the results un-modified.
        """
        return results

    def query_columns(self, colnames=None, chunksize=None,
                      obs_metadata=None):
        """Execute a query

        Parameters
        ----------
        colnames : list or None
            a list of valid column names, corresponding to entries in the
            `columns` class attribute.  If not specified, all columns are
            queried.
        chunksize : int (optional)
            if specified, then return an iterator object to query the database,
            each time returning the next `chunksize` elements.  If not
            specified, all matching results will be returned.
        obs_metadata : object (optional)
            an observation metadata object which has a "filter" method, which
            will add a filter string to the query.

        Returns
        -------
        result : list or iterator
            If chunksize is not specified, then result is a list of all
            items which match the specified query.  If chunksize is specified,
            then result is an iterator over lists of the given size.
        """
        query = self._get_column_query(colnames)

        if obs_metadata is not None:
            query = obs_metadata.filter(query)

        if chunksize is None:
            exec_query = self.session.execute(query)
            return self._postprocess_results(exec_query.fetchall())
        else:
            return ChunkIterator(self, chunksize)


class StarObj(DBObject):
    # XXX: this is incomplete.  We need to use all the column values from
    #      the requiredFields file.
    objid = 'msstars'
    tableid = 'starsMSRGB_forceseek'
    idColName = 'simobjid'
    appendint = 4
    spatialModel = 'POINT'
    columns = ['id', 'umag', 'gmag', 'rmag', 'imag', 'zmag',
               'raJ2000', 'decJ2000', 'sedFilename']
    column_map = {'id':'simobjid',
                  'raJ2000':'ra*PI()/180.',
                  'decJ2000':'decl*PI()/180.',
                  'sedFilename':'sedfilename'}


if __name__ == '__main__':
    #star = StarObj()
    star = DBObject.from_objid('msstars')

    obs_metadata = ObservationMetaData(circ_bounds=dict(ra=2.0,
                                                        dec=5.0,
                                                        radius=1.0))
    print star.query_columns(obs_metadata=obs_metadata)
