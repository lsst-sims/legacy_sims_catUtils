import warnings
import math
import numpy
from collections import OrderedDict
from lsst.sims.catalogs.measures.instance import\
        InstanceCatalog

from sqlalchemy.orm import scoped_session, sessionmaker, mapper
from sqlalchemy.sql import expression
from sqlalchemy import (create_engine, ThreadLocalMetaData, MetaData,
                        Table, Column, BigInteger)

#The documentation at http://docs.sqlalchemy.org/en/rel_0_7/core/types.html#sqlalchemy.types.Numeric
#suggests using the cdecimal module.  Since it is not standard, import decimal.
#TODO: test for cdecimal and use it if it exists.
import decimal

#------------------------------------------------------------
# Iterator for database chunks
class ChunkIterator(object):
    """Iterator for query chunks"""
    def __init__(self, dbobj, query, chunk_size):
        self.dbobj = dbobj
        self.exec_query = dbobj.session.execute(query)
        self.chunk_size = chunk_size

    def __iter__(self):
        return self

    def next(self):
	if self.chunk_size is None and not self.exec_query.closed:
	    chunk = self.exec_query.fetchall()
	    if len(chunk) == 0:
                raise StopIteration
            return self.dbobj._postprocess_results(chunk)
        elif self.chunk_size is None:
	    raise StopIteration
	else:
            chunk = self.exec_query.fetchmany(self.chunk_size)
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

        return super(DBObjectMeta, cls).__init__(name, bases, dct)

    def __str__(cls):
        dbObjects = cls.registry.keys()
        outstr = "++++++++++++++++++++++++++++++++++++++++++++++\n"+\
                 "Registered object types are:\n"
        for dbObject in dbObjects:
            outstr += "%s\n"%(dbObject)
        outstr += "\n\n"
        outstr += "To query the possible column names do:\n"
        outstr += "$> DBObject.from_objid([name]).show_mapped_columns()\n"
        outstr += "+++++++++++++++++++++++++++++++++++++++++++++"
        return outstr


class DBObject(object):
    """Database Object base class

    """
    __metaclass__ = DBObjectMeta
    objid = None
    tableid = None
    idColKey = None
    appendint = None
    spatialModel = None
    columns = None
    raColName = None
    decColName = None
    #: This is the default address.  Simply change this in the class definition for other
    #: endpoints.
    dbAddress = "mssql+pymssql://LSST-2:L$$TUser@fatboy.npl.washington.edu:1433/LSST"
    #: Mapping of DDL types to python types.  Strings are assumed to be 256 characters
    #: this can be overridden by modifying the dbTypeMap or by making a custom columns
    #: list.
    dbTypeMap = {'BIGINT':(int,), 'BOOLEAN':(bool,), 'FLOAT':(float,), 'INTEGER':(int,),\
                 'NUMERIC':(decimal.Decimal,), 'SMALLINT':(int,), 'TINYINT':(int,), 'VARCHAR':(str, 256),\
                 'TEXT':(str, 256), 'CLOB':(str, 256), 'NVARCHAR':(str, 256),\
                 'NCLOB':(unicode, 256), 'NTEXT':(unicode, 256), 'CHAR':(str, 1), 'INT':(int,),\
                 'REAL':(float,)}

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
                             "define objid and tableid.")
        if (self.appendint is None) or (self.spatialModel is None):
            warnings.warn("Either appendint or spatialModel has not "
                          "been set.  Input files for phosim are not "
                          "possible.")

        if address is None:
            address = self.getDbAddress()

        self.dbAddress = address

        self._connect_to_engine()
        self._get_table()

        #Need to do this after the table is instantiated so that
        #the default columns can be filled from the table object.
        if self.columns is None:
            self._make_default_columns()
        # build column mapping and type mapping dicts from columns
        self._make_column_map()
        self._make_type_map()

    def show_mapped_columns(self):
        for col in self.columnMap.keys():
            print "%s -- %s"%(col, self.typeMap[col][0].__name__)

    def show_db_columns(self):
        for col in self.table.c.keys():
            print "%s -- %s"%(col, self.table.c[col].type.__visit_name__)


    def getCatalog(self, ftype, *args, **kwargs):
        return InstanceCatalog.new_catalog(ftype, self, *args, **kwargs)

    def getDbAddress(self):
        return self.dbAddress

    def getObjectTypeId(self):
        return self.appendint

    def getSpatialModel(self):
        return self.spatialModel

    def _get_table(self):
        '''
        self.table = Table(self.tableid, self.metadata,
                           Column(self.columnMap[self.idColKey], BigInteger, primary_key=True),
                           autoload=True)
        '''
        self.table = Table(self.tableid, self.metadata,
                           autoload=True)

    def _connect_to_engine(self):
        """create and connect to a database engine"""
        self.engine = create_engine(self.dbAddress, echo=False)
        self.session = scoped_session(sessionmaker(autoflush=True, 
                                                   bind=self.engine))
        self.metadata = MetaData()
        self.metadata.bind = self.engine
    def _make_column_map(self):
        self.columnMap = OrderedDict([(el[0], el[1] if el[1] else el[0]) 
                                     for el in self.columns])
    def _make_type_map(self):
        self.typeMap = OrderedDict([(el[0], el[2:] if len(el)> 2 else (float,))
                                   for el in self.columns])

    def _make_default_columns(self):
        self.columns = []
        for col in self.table.c.keys():
            dbtypestr = self.table.c[col].type.__visit_name__
            if dbtypestr in self.dbTypeMap.keys():
                self.columns.append((col, col)+self.dbTypeMap[dbtypestr])
            else:
                warnings.warn("Can't create default column for %s.  There is no mapping "%(col)+\
                              "for type %s.  Modify the dbTypeMap, or make a custom columns "%(dbtypestr)+\
                              "list.")

    def _get_column_query(self, colnames=None):
        """Given a list of valid column names, return the query object"""
        if colnames is None:
            colnames = [k for k in self.columnMap.keys()]
        try:
            vals = [self.columnMap[k] for k in colnames]
        except KeyError:
            for c in colnames:
                if c in self.columnMap.keys():
                    continue
                else:
                    print "%s not in columnMap"%(c)
            raise ValueError('entries in colnames must be in self.columnMap')

        # Get the first query
        idColName = self.columnMap[self.idColKey]
        if idColName in vals:
            idLabel = self.idColKey
        else:
            idLabel = idColName

        query = self.session.query(self.table.c[idColName].label(idLabel))

        for col, val in zip(colnames, vals):
            if val is idColName:
                continue
            query = query.add_column(expression.literal_column(val).label(col))

        return query

    def filter(self, query, circ_bounds=None, box_bounds=None):
        """Filter the query by the associated metadata"""
        on_clause = self.to_SQL(circ_bounds, box_bounds)
        if on_clause:
            query = query.filter(on_clause)
        return query

    def to_SQL(self, circ_bounds=None, box_bounds=None):
        if circ_bounds and box_bounds:
            raise ValueError("circ_bounds and box_bounds should not both be set")
        constraint = ""
        if box_bounds:
            bb = box_bounds
            constraint = self.box_bound_constraint(bb['ra_min'],
                                                    bb['ra_max'],
                                                    bb['dec_min'],
                                                    bb['dec_max'],
                                                    self.raColName,
                                                    self.decColName)
        if circ_bounds:
            cb = circ_bounds
            constraint = self.circle_bound_constraint(cb['ra'], cb['dec'],
                                                       cb['radius'],
                                                       self.raColName, self.decColName)
        return constraint

    @staticmethod
    def mjd_constraint(mjd_bounds, MJDname):
        raise NotImplementedError("This is better done using a constraint at run time")

    @staticmethod
    def box_bound_constraint(RAmin, RAmax, DECmin, DECmax,
                             RAname, DECname):
        #KSK:  I don't know exactly what we do here.  This is in code, but operating
        #on a database is it less confusing to work in degrees or radians?
        #(RAmin, RAmax, DECmin, DECmax) = map(math.radians,
        #                                     (RAmin, RAmax, DECmin, DECmax))

        #Special case where the whole region is selected
        if RAmin < 0 and RAmax > 360.:
            bound = "%s between %f and %f" % (DECname, DECmin, DECmax)
            return bound

        RAmin %= 360.
        RAmax %= 360.
        if RAmin > RAmax:
            # XXX is this right?  It seems strange.
            bound = ("%s not between %f and %f and %s between %f and %f"
                     % (RAname, RAmax, RAmin, 
                        DECname, DECmin, DECmax))
        else:
            bound = ("%s between %f and %f and %s between %f and %f"
                     % (RAname, RAmin, RAmax, DECname, DECmin, DECmax))

        return bound

    @staticmethod
    def circle_bound_constraint(RA, DEC, radius,
                                RAname, DECname):
        RAmax = RA + radius / math.cos(math.radians(DEC))
        RAmin = RA - radius / math.cos(math.radians(DEC))
        DECmax = DEC + radius
        DECmin = DEC - radius
        return DBObject.box_bound_constraint(RAmin, RAmax,
                                                        DECmin, DECmax,
                                                        RAname, DECname)    

    def _final_pass(self, results):
        """ Make final modifications to a set of data before returning it to the user
        
        **Parameters**
        
            * results : a structured array constructed from the result set from a query

        **Returns**
        
            * results : a potentially modified structured array.  The default is to do nothing.
        
        """
        return results

    def _postprocess_results(self, results):
        """Post-process the query results to put then
        in a structured array.
  
        **Parameters**

            * results : a result set as returned by execution of the query

        **Returns**

            * _final_pass(retresults) : the result of calling the _final_pass method on a
              structured array constructed from the query data.

        """
        if len(results) > 0:
            cols = [str(k) for k in results[0].keys()]
        else:
            return results
        dtype = numpy.dtype([(k,)+self.typeMap[k] for k in cols])
        retresults = numpy.zeros((len(results),), dtype=dtype)
        for i, result in enumerate(results):
            for k in cols:
                retresults[i][k] = result[k]
        return self._final_pass(retresults)

    def query_columns(self, colnames=None, chunk_size=None,
                      obs_metadata=None, constraint=None):
        """Execute a query

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

        **Returns**

            * result : list or iterator
              If chunk_size is not specified, then result is a list of all
              items which match the specified query.  If chunk_size is specified,
              then result is an iterator over lists of the given size.

        """
        query = self._get_column_query(colnames)

        if obs_metadata is not None:
            query = self.filter(query, circ_bounds=obs_metadata.circ_bounds, 
                    box_bounds=obs_metadata.box_bounds)

        if constraint is not None:
            query = query.filter(constraint)
        return ChunkIterator(self, query, chunk_size)
