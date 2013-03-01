import warnings
import math
from sqlalchemy.orm import scoped_session, sessionmaker, mapper
from sqlalchemy.sql import expression
from sqlalchemy import (create_engine, ThreadLocalMetaData, MetaData,
                        Table, Column, BigInteger)


DEFAULT_ADDRESS = "mssql+pymssql://LSST-2:L$$TUser@fatboy.npl.washington.edu:1433/LSST"


# XXX: create an OpSim metadata interface class?


#------------------------------------------------------------
# Iterator for database chunks
class ChunkIterator(object):
    """Iterator for query chunks"""
    def __init__(self, exec_query, chunksize):
        self.exec_query = exec_query
        self.chunksize = chunksize

    def __iter__(self):
        return self

    def next(self):
        chunk = self.exec_query.fetchmany(self.chunksize)
        if len(chunk) == 0:
            raise StopIteration
        return chunk

#------------------------------------------------------------
# Query utilities
def box_bound(RAmin, RAmax, DECmin, DECmax,
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


def circle_bound(RA, DEC, radius,
                 RAname='ra', DECname='decl'):
    RAmax = RA + radius / math.cos(math.radians(DEC))
    RAmin = RA - radius / math.cos(math.radians(DEC))
    DECmax = DEC + radius
    DECmin = DEC - radius

    return box_bound(RAmin, RAmax, DECmin, DECmax, RAname, DECname)


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

            if not hasattr(cls, 'columns'):
                raise ValueError()

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
        # XXX: We've hard-coded simobjid here: this might not be correct
        #      in general.  We should find a better way to do this.
        self.table = Table(self.tableid, self.metadata,
                           Column('simobjid', BigInteger, primary_key=True),
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
        # XXX: this will fail if the first of the queries is a compound
        #      expression like 'ra*PI()/180.'  We should do this a better
        #      way, perhaps using the primary column (related to problem
        #      in _get_table(), above)
        query = self.session.query(self.table.c[vals[0]].label(colnames[0]))

        for col, val in zip(colnames[1:], vals[1:]):
            query = query.add_column(expression.literal_column(val).label(col))

        return query

    def query_columns(self, colnames=None, chunksize=None,
                      circ_bounds=None, box_bounds=None):
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
        circ_bounds : tuple (optional)
            if specified, then limit the query to the specified circular
            spatial range. circ_bounds = (RAcenter, DECcenter, radius),
            measured in degrees.
        box_bounds :tuple (optional)
            if specified, then limit the query to the specified quadrilateral
            spatial range.  box_bounds = (RAmin, RAmax, DECmin, DECmax),
            measured in degrees

        Returns
        -------
        result : list or iterator
            If chunksize is not specified, then result is a list of all
            items which match the specified query.  If chunksize is specified,
            then result is an iterator over lists of the given size.
        """
        query = self._get_column_query(colnames)

        if circ_bounds is not None:
            RA, DEC, radius = circ_bounds
            query = query.filter(circle_bound(RA, DEC, radius))

        if box_bounds is not None:
            RAmin, RAmax, DECmin, DECmax = box_bounds
            query = query.filter(box_bound(RAmin, RAmax, DECmin, DECmax))

        exec_query = self.session.execute(query)

        if chunksize is None:
            return exec_query.fetchall()
        else:
            return ChunkIterator(exec_query, chunksize)


class StarObj(DBObject):
    # XXX: this is incomplete.  We need to use all the column values from
    #      the requiredFields file.
    objid = 'msstars'
    tableid = 'starsMSRGB_forceseek'
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
    print star.query_columns(circ_bounds=(2.0, 5.0, 1.0))
