import warnings
from lsst.sims.catalogs.generation.db import ChunkIterator, DBObject, ObservationMetaData

class SolarSystemObj(DBObject):
    objid = 'ssm'
    # There is no materialized table since this is a table valued function
    tableid = ''
    objectTypeId = 40

    doRunTest = True
    testObservationMetaData = ObservationMetaData(circ_bounds=dict(ra=0., dec=0., radius=1.0),
                                                  mjd=51200., bandpassName='r')

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
            ('sedFilename', 'sed_filename', unicode, 40),
            ('velRa', 'dradt*PI()/180.'),
            ('velDec', 'ddecldt*PI()/180.')]

    def _get_column_query(self, colnames=None):
        raise NotImplementedError("We are calling a stored procedure so "
                                  "no need to loop over columns")

    def _get_table(self):
        #This needs to be overridden since __init__() calls the default impl.
        return None

    def getIdColKey(self):
        return 'ssmid'

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
        circ_bounds = None
        box_bounds = None

        if obs_metadata is not None:
            if obs_metadata.circ_bounds is not None:
                 circ_bounds = obs_metadata.circ_bounds
            if obs_metadata.box_bounds is not None:
                box_bounds = obs_metadata.box_bounds
        
        if circ_bounds is not None:
            RA = circ_bounds['ra']
            DEC = circ_bounds['dec']
            radius = circ_bounds['radius']
        elif box_bounds is not None:
            raise NotImplementedError("There is no way to specify a box "
                                      "for the ssm stored procedure")
        else:
            RA, DEC, radius = (180., 0., 180.)
            warnings.warn("Searching over entire sky "
                          "since no circ_bounds specified. "
                          "This could be a very bad idea "
                          "if the database is large")

        query = "select %s from [LSST].[dbo].fSSMAll("%(mappedcolnames)+\
                "%f, %f, %f, %f)"%(obs_metadata.mjd, RA, DEC, radius*60.)

        if constraint is not None:
            query += "where %s"%(constraint)

        return ChunkIterator(self, query, chunk_size)
