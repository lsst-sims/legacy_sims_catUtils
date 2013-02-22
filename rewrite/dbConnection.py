import warnings
from sqlalchemy.orm import scoped_session, sessionmaker, mapper
from sqlalchemy.sql import expression
from sqlalchemy import create_engine
from sqlalchemy import ThreadLocalMetaData
import sqlalchemy.databases as sd 
from sqlalchemy import MetaData
from sqlalchemy import Table
from sqlalchemy.ext.sqlsoup import SqlSoup
from sqlalchemy import exc as sa_exc


DEFAULT_ADDRESS = "mssql+pymssql://LSST-2:L$$TUser@fatboy.npl.washington.edu:1433/LSST"


requirementsDict = dict(
    msstars={'id':'simobjid',
             'umag':'umag'}
    )


# OpSim metadata interface class?

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

        # add this class to the registry
        if cls.objid in cls.registry:
            warnings.warn('duplicate object id %s specified' % cls.objid)
        cls.registry[cls.objid] = cls

        # build requirements dict from columns and column_map
        cls.requirements = dict([(key, cls.column_map.get(key, key))
                                 for key in cls.columns])
            
        return super(DBObjectMeta, cls).__init__(name, bases, dct)


class DBObject(object):
    __metaclass__ = DBObjectMeta
    objid = None
    tableid = None
    columns = []
    column_map = {}
    requires_tiling = False

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
        if (self.objid is None) or (self.tableid is None):
            raise ValueError("DBObject must be subclassed, and "
                             "define objid & tableid.")

        if address is None:
            self.address = DEFAULT_ADDRESS
        else:
            self.address = address
        self._connect_to_engine()
        self._get_table()

    def _get_table(self):
        table = Table(self.tableid, self.metadata, autoload=True)
        self.table = self.db.map(table,
                                 primary_key=[table.c.simobjid])

    def _connect_to_engine(self):
        self.engine = create_engine(self.address, echo=False)
        self.session = scoped_session(sessionmaker(autoflush=True, 
                                                   bind=self.engine))
        self.metadata = MetaData()
        self.metadata.bind = self.engine
        self.db = SqlSoup(self.metadata)

    def _get_column_query(self, colnames):
        """Given a list of valid column names, return the query object"""
        try:
            vals = [self.requirements[col] for col in colnames]
        except KeyError:
            raise ValueError('entries in colnames must be in self.columns')

        # Get first query
        # XXX: this will fail if the first of the queries is a compound
        #      expression like 'ra*PI()/180.'  We should do this a better
        #      way, perhaps using the primary column?
        query = self.session.query(getattr(self.table,
                                           vals[0]).label(colnames[0]))
        for col, val in zip(colnames[1:], vals[1:]):
            query = query.add_column(expression.literal_column(val).label(col))

        return query

    def add_spatial_query(self):
        # XXX: add spatial query
        pass

    def compute_tiling(self):
        # XXX: compute tiling
        pass

    def query_columns(self, colnames):
        query = self._get_column_query(colnames)
        # XXX: execute query here
        return query


class StarObj(DBObject):
    objid = 'msstars'
    tableid = 'starsMSRGB_forceseek'
    columns = ['id', 'umag', 'gmag', 'rmag', 'imag', 'zmag',
               'raJ2000', 'decJ2000']
    column_map = {'id':'simobjid',
                  'raJ2000':'ra*PI()/180.',
                  'decJ2000':'decl*PI()/180.'}
    requires_tiling = False
    

if __name__ == '__main__':
    #star = StarObj()
    star = DBObject.from_objid('msstars')
    print star.query_columns(['id', 'gmag', 'rmag', 'raJ2000', 'decJ2000'])
