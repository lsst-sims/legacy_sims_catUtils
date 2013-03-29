import warnings
import math
import numpy

from sqlalchemy.orm import scoped_session, sessionmaker, mapper
from sqlalchemy.sql import expression
from sqlalchemy import (create_engine, ThreadLocalMetaData, MetaData,
                        Table, Column, BigInteger)


# TODO : we need to clean up "columns", "column_map", and "requirements".
#        currently it's very confusing, and a lot of things are duplicated.
#        Initially, I thought that "columns" would specify the ordering
#        of things, and "column_map" would be used for only a handful of
#        the columns.  After seeing Simon's GalaxyObj implementations, it's
#        clear that this assumption was wrong: nearly every column needs to
#        be mapped, and additionally we need to specify type info for nearly
#        every column.
#
#        For simplicity, we should have the class define a single array with
#        column names & attributes.  It could look like this:
#
#        columns = [('galid', 'galid', str, 30),
#                   ('raJ2000', 'ra*PI()/180., float),
#                   ('decJ2000', 'dec*PI()/180.'),  # float by default
#                   ...]
#
#        The metaclass would take this column list, and turn it into two
#        ordered dictionaries, (using `from collections import OrderedDict`)
#        one consisting of the column mapping, and one consisting of the
#        type mapping.


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
    
    def __init__(self, obshistid=None, circ_bounds=None, box_bounds=None, mjd=None):
        if obshistid is not None:
            raise NotImplementedError()
        else:
            if circ_bounds is not None and box_bounds is not None:
                raise ValueError()
            self.circ_bounds = circ_bounds
            self.box_bounds = box_bounds
            self.mjd = mjd



#------------------------------------------------------------
# Iterator for database chunks
class ChunkIterator(object):
    """Iterator for query chunks"""
    def __init__(self, exec_query, chunk_size):
        self.dbobj = dbobj
        self.exec_query = dbobj.session.execute(query)
        self.chunk_size = chunk_size

    def __iter__(self):
        return self

    def next(self):
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

            # build requirements dict from columns
            cls.requirements = dict([(k, cls.columns[k][0])
                                     for k in cls.columns.keys()])
            # build dtype from columns            
            cls.dtype = numpy.dtype([(k,)+cls.columns[k][1:]
                                     for k in cls.columns.keys()])
            
        return super(DBObjectMeta, cls).__init__(name, bases, dct)


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
    mjdColName = None

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
        if self.columns is None:
            raise ValueError("DBObject must be subclasses, and define "
                             "columns.  The columns variable is a list "
                             "of tuples suitable for parsing by numpy.dtype")

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
                           Column(self.columns[self.idColKey][0], BigInteger, primary_key=True),
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
        #XXXX This is not o.k. as this is now a dict and as such is 
        #not ordered.
        if colnames is None:
            colnames = [k for k in self.columns.keys()]
        try:
            vals = [self.requirements[col] for col in colnames]
        except KeyError:
            raise ValueError('entries in colnames must be in self.columns')

        # Get the first query
        idColName = self.columns[self.idColKey][0]
        if idColName in vals:
            idLabel = colnames[vals.index(idColName)]
        else:
            idLabel = 'id'

        query = self.session.query(self.table.c[idColName].label(idLabel))

        for col, val in zip(colnames, vals):
            if val is idColName:
                continue
            query = query.add_column(expression.literal_column(val).label(col))

        return query

    def filter(self, query, circ_bounds=None, box_bounds=None, mjd=None):
        """Filter the query by the associated metadata"""
        on_clause = self.to_SQL(circ_bounds, box_bounds, mjd)
        if on_clause:
            query = query.filter(on_clause)
        return query

    def to_SQL(self, circ_bounds=None, box_bounds=None, mjd=None):
        constraint = ""
        if box_bounds is not None:
            bb = box_bounds
            constraint += self.box_bound_constraint(bb['ra_min'],
                                                    bb['ra_max'],
                                                    bb['dec_min'],
                                                    bb['dec_max'],
						    self.raColName,
						    self.decColName)
        if circ_bounds is not None:
            cb = circ_bounds
            constraint += self.circle_bound_constraint(cb['ra'], cb['dec'],
                                                       cb['radius'],
                                                       self.raColName, self.decColName)
	#KSK: Make MJD self consistent with other column labelings
        if mjd is not None:
            constraint += self.mjd_constraint(mjd, self.mjdColName)
            
        return constraint

    @staticmethod
    def mdj_constraint(MJD, MJDname):
        raise NotImplementedError("haven't implemented MJD bound yet")

    @staticmethod
    def box_bound_constraint(RAmin, RAmax, DECmin, DECmax,
                             RAname, DECname):
        #KSK:  I don't know exactly what we do here.  This is in code, but operating
        #on a database is it less confusing to work in degrees or radians?
        #(RAmin, RAmax, DECmin, DECmax) = map(math.radians,
        #                                     (RAmin, RAmax, DECmin, DECmax))

        if RAmin < 0 and RAmax > 360.:
            bound = "%s between %f and %f" % (DECname, DECmin, DECmax)

        elif RAmin < 0 and RAmax <= 360.:
            # XXX is this right?  It seems strange.
            bound = ("%s not between %f and %f and %s between %f and %f"
                     % (RAname, RAmin % (360.), RAmax,
                        DECname, DECmin, DECmax))

        elif RAmin >= 0 and RAmax > 2. * math.pi:
            bound = ("%s not between %f and %f and %s between %f and %f" 
                     % (RAname, RAmin, RAmax % (360.),
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


    def _postprocess_results(self, results):
        """Post-process the query results.
        This can be overridden in derived classes: by default it returns
        the results un-modified.
        """
        retresults = numpy.zeros((len(results),),dtype=self.dtype)
        for i, result in enumerate(results):
            for k in self.requirements.keys():
                retresults[i][k] = result[k]
        return retresults

    def query_columns(self, colnames=None, chunk_size=None,
                      obs_metadata=None, constraint=None):
        """Execute a query

        Parameters
        ----------
        colnames : list or None
            a list of valid column names, corresponding to entries in the
            `columns` class attribute.  If not specified, all columns are
            queried.
        chunk_size : int (optional)
            if specified, then return an iterator object to query the database,
            each time returning the next `chunk_size` elements.  If not
            specified, all matching results will be returned.
        obs_metadata : object (optional)
            an observation metadata object which has a "filter" method, which
            will add a filter string to the query.

        Returns
        -------
        result : list or iterator
            If chunk_size is not specified, then result is a list of all
            items which match the specified query.  If chunk_size is specified,
            then result is an iterator over lists of the given size.
        """
        query = self._get_column_query(colnames)
        if obs_metadata is not None:
            query = self.filter(query, circ_bounds=obs_metadata.circ_bounds, 
                    box_bounds=obs_metadata.box_bounds, mjd=obs_metadata.mjd)

        if constraint is not None:
            query = query.filter(constraint)
        if chunk_size is None:
            exec_query = self.session.execute(query)
            return self._postprocess_results(exec_query.fetchall())
        else:
            return ChunkIterator(self, chunk_size)


class StarObj(DBObject):
    # XXX: this is incomplete.  We need to use all the column values from
    #      the requiredFields file.
    objid = 'msstars'
    tableid = 'starsMSRGB_forceseek'
    idColKey = 'id'
    raColName = 'ra'
    decColName = 'decl'
    appendint = 4
    spatialModel = 'POINT'
    #These types should be matched to the database.
    columns = {'id':('simobjid', int),
               'umag':('umag', float),
               'gmag':('gmag', float),
               'rmag':('rmag', float),
               'imag':('imag', float),
               'zmag':('zmag', float),
               'raJ2000':('ra*PI()/180.', float),
               'decJ2000':('decl*PI()/180.', float),
               'sedFilename':('sedfilename', unicode, 40)}

if __name__ == '__main__':
    #star = StarObj()
    star = DBObject.from_objid('msstars')
    obs_metadata = ObservationMetaData(circ_bounds=dict(ra=2.0,
                                                        dec=5.0,
                                                        radius=1.0))
    print star.query_columns(obs_metadata=obs_metadata, constraint="rmag < 21.")
