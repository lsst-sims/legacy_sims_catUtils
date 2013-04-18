import warnings
import math
import numpy
import os
from dbConnection import ChunkIterator, DBObject
from MetaDataDBObject import ObservationMetaData
from sqlalchemy import Table, Column, BigInteger, MetaData

class MyGalaxyObj(DBObject):
    
    objid = 'mygalaxyBase'
    #This is the base table for the galaxies
    tableid = 'galaxies'
    idColKey = 'galid'
    raColName = 'ra'
    decColName = 'decl'
    appendint = 9
    #There is no spatial model available for coadded galaxies 
    #This will cause a warning, but it just means we can't make
    #TRIM files with this object.    
    spatialModel = None



    #The following maps column names to database schema.  The tuples
    #must be at least length 2.  If column name is the same as the name
    #in the DB the mapping element may be None.  The rest of the tuple
    #should be formatted like a numpy.dtype.  If ommitted, the dtype
    #is assumed to be float.
    columns = [('galid', None, str, 30),
            ('raJ2000', 'ra*PI()/180.'),
            ('decJ2000', 'decl*PI()/180.'),
            ('raJ2000Bulge', 'bra*PI()/180.'),
            ('decJ2000Bulge', 'bdec*PI()/180.'),
            ('raJ2000Disk', 'dra*PI()/180.'),
            ('decJ2000Disk', 'ddec*PI()/180.'),
            ('raJ2000Agn', 'agnra*PI()/180.'),
            ('decJ2000Agn', 'agndec*PI()/180.'),
            ('magNormBulge', 'magnorm_bulge'),
            ('magNormDisk', 'magnorm_disk'),
            ('magNormAgn', 'magnorm_agn'),
            ('sedFilenameBulge', 'sedname_bulge', unicode, 40),
            ('sedFilenameDisk', 'sedname_disk', unicode, 40),
            ('sedFilenameAgn', 'sedname_agn', unicode, 40),
            ('majorAxisBulge', 'a_b'),
            ('minorAxisBulge', 'b_b'),
            ('positionAngleBulge', 'pa_bulge'),
            ('sindexBulge', 'bulge_n', int),
            ('majorAxisDisk', 'a_d'),
            ('minorAxisDisk', 'b_d'),
            ('positionAngleDisk', 'pa_disk'),
            ('sindexDisk', 'disk_n', int),
            ('internalExtinctionModelBulge', 'ext_model_b', str, 3),
            ('internalAvBulge', 'av_b'),
            ('internalRvBulge', 'rv_b'),
            ('internalExtinctionModelDisk', 'ext_model_d', str, 3),
            ('internalAvDisk', 'av_d'),
            ('internalRvDisk', 'rv_d'),
            ('redshift', None),
            ('radialVelocity', 'rad_vel'),
            ('lsst_u', 'u_ab'),
            ('lsst_g', 'g_ab'),
            ('lsst_r', 'r_ab'),
            ('lsst_i', 'i_ab'),
            ('lsst_z', 'z_ab'),
            ('lsst_y', 'y_ab')]

    def getDbAddress(self):
        home_path = os.getenv("HOME")
        f=open("%s/dbLogin"%(home_path),"r")
        return (f.readline()).strip()
    

class GalaxyObj(DBObject):
    objid = 'galaxyBase'
    #This is the base table for the galaxies
    tableid = 'galaxy'
    idColKey = 'galid'
    raColName = '((CAST(ra AS NUMERIC(9,6))%360.)+360.)%360.'
    decColName = 'dec'
    appendint = 9
    #There is no spatial model available for coadded galaxies 
    #This will cause a warning, but it just means we can't make
    #TRIM files with this object.
    spatialModel = None

    #The following maps column names to database schema.  The tuples
    #must be at least length 2.  If column name is the same as the name
    #in the DB the mapping element may be None.  The rest of the tuple
    #should be formatted like a numpy.dtype.  If ommitted, the dtype
    #is assumed to be float.
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
            ('sedFilenameBulge', 'sedname_bulge', unicode, 40),
            ('sedFilenameDisk', 'sedname_disk', unicode, 40),
            ('sedFilenameAgn', 'sedname_agn', unicode, 40),
            ('majorAxisBulge', 'a_b'),
            ('minorAxisBulge', 'b_b'),
            ('positionAngleBulge', 'pa_bulge'),
            ('sindexBulge', 'bulge_n', int),
            ('majorAxisDisk', 'a_d'),
            ('minorAxisDisk', 'b_d'),
            ('positionAngleDisk', 'pa_disk'),
            ('sindexDisk', 'disk_n', int),
            ('internalExtinctionModelBulge', 'ext_model_b', str, 3),
            ('internalAvBulge', 'av_b'),
            ('internalRvBulge', 'rv_b'),
            ('internalExtinctionModelDisk', 'ext_model_d', str, 3),
            ('internalAvDisk', 'av_d'),
            ('internalRvDisk', 'rv_d'),
            ('redshift', None),
            ('radialVelocity', 'rad_vel'),
            ('lsst_u', 'u_ab'),
            ('lsst_g', 'g_ab'),
            ('lsst_r', 'r_ab'),
            ('lsst_i', 'i_ab'),
            ('lsst_z', 'z_ab'),
            ('lsst_y', 'y_ab')]

    def _final_pass(self, results):
        #This is to map the values from 0 - 2*PI() as ra goes negative currently
        results['raJ2000'] = results['raJ2000']%(numpy.pi*2.)
        results['raJ2000Bulge'] = results['raJ2000Bulge']%(numpy.pi*2.)
        results['raJ2000Disk'] = results['raJ2000Disk']%(numpy.pi*2.)
        results['raJ2000Agn'] = results['raJ2000Agn']%(numpy.pi*2.)
        return results



class GalaxyTileObj(DBObject):
    objid = 'galaxyTiled'
    #This is the base table for the galaxies
    tableid = 'galaxy'
    idColKey = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 9
    #There is no spatial model available for coadded galaxies 
    spatialModel = None

    #The following maps column names to database schema.  The tuples
    #must be at least length 2.  If column name is the same as the name
    #in the DB the mapping element may be None.  The rest of the tuple
    #should be formatted like a numpy.dtype.  If ommitted, the dtype
    #is assumed to be float.
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
            ('sedFilenameBulge', 'sedname_bulge', unicode, 40),
            ('sedFilenameDisk', 'sedname_disk', unicode, 40),
            ('sedFilenameAgn', 'sedname_agn', unicode, 40),
            ('majorAxisBulge', 'a_b'),
            ('minorAxisBulge', 'b_b'),
            ('positionAngleBulge', 'pa_bulge'),
            ('sindexBulge', 'bulge_n', int),
            ('majorAxisDisk', 'a_d'),
            ('minorAxisDisk', 'b_d'),
            ('positionAngleDisk', 'pa_disk'),
            ('sindexDisk', 'disk_n', int),
            ('internalExtinctionModelBulge', 'ext_model_b', str, 3),
            ('internalAvBulge', 'av_b'),
            ('internalRvBulge', 'rv_b'),
            ('internalExtinctionModelDisk', 'ext_model_d', str, 3),
            ('internalAvDisk', 'av_d'),
            ('internalRvDisk', 'rv_d'),
            ('redshift', None),
            ('radialVelocity', 'rad_vel'),
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
        Parameters
        ----------
        results : Structured array of results from query

        Returns
        -------
        results : Modified structured array
        """

        results['raJ2000'] = numpy.radians(results['raJ2000'])
        results['decJ2000'] = numpy.radians(results['decJ2000'])
        return results

    def query_columns(self, colnames=None, chunk_size=None, obs_metadata=None, constraint=None):
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
            object containing information on the observation including the region of the sky
            to query and time of the observation.
        constraint : string (optional)
            if specified, the predicate is added to the query verbatim using AND

        Returns
        -------
        result : structured array or iterator
            If chunksize is not specified, then result is a structured array of all
            items which match the specified query with columns named by the column
            names in the columns class attribute.  If chunksize is specified,
            then result is an iterator over structured arrays of the given size.
        """
        if colnames is None:
            colnames = [k for k in self.columnMap.keys()]
        
        #We know that galtileid comes back with the query, but we don't want 
        #to add it to the query since it's generated on the fly.
        while 'galtileid' in colnames:
            colnames.remove('galtileid')
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
                                      "for the galaxy stored procedure")
        else:
            RA, DEC, radius = (180., 0., 180.)
            warnings.warn("Searching over entire sky "
                          "since no circ_bounds specified. "
                          "This could be a very bad idea "
                          "if the database is large")


        query = "EXECUTE [LSST].[dbo].[GalaxySearchSpecColsConstraint2013]\
               @RaSearch = %f, @DecSearch = %f, @apertureRadius = %f,\
               @ColumnNames = '%s'" % (RA, DEC, radius * 60., mappedcolnames)

        if constraint is not None:
            query += ", @WhereClause = '%s'"%(constraint)

        exec_query = self.session.execute(query)
        if chunk_size is None:
            return self._postprocess_results(exec_query.fetchall())
        else:
            return ChunkIterator(self, query, chunk_size)

class GalaxyBulgeObj(GalaxyTileObj):
    objid = 'galaxyBulge'
    #This is the base table for the galaxies
    #with a bulge component
    tableid = 'galaxy_bulge'
    idColKey = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 1
    spatialModel = 'SERSIC2D'
    #The following maps column names to database schema.  The tuples
    #must be at least length 2.  If column name is the same as the name
    #in the DB the mapping element may be None.  The rest of the tuple
    #should be formatted like a numpy.dtype.  If ommitted, the dtype
    #is assumed to be float.
    columns = [('galtileid', None, numpy.int64),
            ('galid', None, str, 30),
            ('componentra','bra*PI()/180.'),
            ('componentdec', 'bdec*PI()/180.'),
            #This is actually a problem with the stored procedure.  We need to be able to map columns other than
            #just ra/dec to raJ2000/decJ2000.  This gets important when we start perturbing the three galaxy components
            ('raJ2000', 'ra'),
            ('decJ2000', 'dec'),
            ('magNorm', 'magnorm_bulge'),
            ('sedFilename', 'sedname_bulge', unicode, 40),
            ('majorAxis', 'a_b*PI()/180.'),
            ('minorAxis', 'b_b*PI()/180.'),
            ('positionAngle', 'pa_bulge*PI()/180.'),
            ('sindex', 'bulge_n', int),
            ('internalExtinctionModel', 'ext_model_b', str, 3),
            ('internalAv', 'av_b'),
            ('internalRv', 'rv_b'),
            ('redshift', None),
            ('radialVelocity', 'rad_vel'),
            ('lsst_u', 'u_ab'),
            ('lsst_g', 'g_ab'),
            ('lsst_r', 'r_ab'),
            ('lsst_i', 'i_ab'),
            ('lsst_z', 'z_ab'),
            ('lsst_y', 'y_ab')]

class GalaxyDiskObj(GalaxyTileObj):
    objid = 'galaxyDisk'
    #This is the base table for the galaxies
    #with a disk component
    tableid = 'galaxy'
    idColKey = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 2
    spatialModel = 'SERSIC2D'
    #The following maps column names to database schema.  The tuples
    #must be at least length 2.  If column name is the same as the name
    #in the DB the mapping element may be None.  The rest of the tuple
    #should be formatted like a numpy.dtype.  If ommitted, the dtype
    #is assumed to be float.
    columns = [('galtileid', None, numpy.int64),
            ('galid', None, str, 30),
            ('componentra','dra*PI()/180.'),
            ('componentdec', 'ddec*PI()/180.'),
            #This is actually a problem with the stored procedure.  We need to be able to map columns other than
            #just ra/dec to raJ2000/decJ2000.  This gets important when we start perturbing the three galaxy components
            ('raJ2000', 'ra'),
            ('decJ2000', 'dec'),
            ('magNorm', 'magnorm_disk'),
            ('sedFilename', 'sedname_disk', unicode, 40),
            ('majorAxis', 'a_d*PI()/180.'),
            ('minorAxis', 'b_d*PI()/180.'),
            ('positionAngle', 'pa_disk*PI()/180.'),
            ('sindex', 'disk_n', int),
            ('internalExtinctionModel', 'ext_model_d', str, 3),
            ('internalAv', 'av_d'),
            ('internalRv', 'rv_d'),
            ('redshift', None),
            ('radialVelocity', 'rad_vel'),
            ('lsst_u', 'u_ab'),
            ('lsst_g', 'g_ab'),
            ('lsst_r', 'r_ab'),
            ('lsst_i', 'i_ab'),
            ('lsst_z', 'z_ab'),
            ('lsst_y', 'y_ab')]

class GalaxyAgnObj(GalaxyTileObj):
    objid = 'galaxyAgn'
    #This is the base table for the galaxies
    #with an agn component
    tableid = 'galaxy_agn'
    idColKey = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 3
    spatialModel = 'ZPOINT'
    #The following maps column names to database schema.  The tuples
    #must be at least length 2.  If column name is the same as the name
    #in the DB the mapping element may be None.  The rest of the tuple
    #should be formatted like a numpy.dtype.  If ommitted, the dtype
    #is assumed to be float.
    columns = [('galtileid', None, numpy.int64),
            ('galid', None, str, 30),
            ('componentra','agnra*PI()/180.'),
            ('componentdec', 'agndec*PI()/180.'),
            #This is actually a problem with the stored procedure.  We need to be able to map columns other than
            #just ra/dec to raJ2000/decJ2000.  This gets important when we start perturbing the three galaxy components
            ('raJ2000', 'ra'),
            ('decJ2000', 'dec'),
            ('magNorm', 'magnorm_agn'),
            ('sedFilename', 'sedname_agn', unicode, 40),
            ('redshift', None),
            ('radialVelocity', 'rad_vel'),
            ('variabilityParameters', 'varParamStr', str, 256),
            ('lsst_u', 'u_ab'),
            ('lsst_g', 'g_ab'),
            ('lsst_r', 'r_ab'),
            ('lsst_i', 'i_ab'),
            ('lsst_z', 'z_ab'),
            ('lsst_y', 'y_ab')]

