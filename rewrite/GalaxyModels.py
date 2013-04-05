import warnings
import math
import numpy
from dbConnection import ChunkIterator, DBObject, ObservationMetaData
from sqlalchemy import Table, Column, BigInteger, MetaData

class GalaxyObj(DBObject):
    objid = 'galaxyBase'
    #This is the base table for the galaxies
    #with a bulge component
    tableid = 'galaxy'
    idColKey = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 9
    #There is no spatial model available for coadded galaxies 
    spatialModel = None
    #Only realy need to specify the mappings that are not 
    #"x as x" but most need to be mapped anyway.
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
    #Only realy need to specify the mappings that are not 
    #"x as x" but most need to be mapped anyway.
    columns = [('galtileid', None, numpy.int64),
            ('galid', None, str, 30),
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

    def _get_column_query(self, colnames=None):
        raise NotImplementedError("We are calling a stored procedure so "
                                  "no need to loop over columns")

    def _final_pass(self, results):
        """Modify the results of raJ2000 and decJ2000 to be in radians"""
	results['raJ2000'] = numpy.radians(results['raJ2000'])
	results['decJ2000'] = numpy.radians(results['decJ2000'])
        return results

    def query_columns(self, colnames=None, chunksize=None, obs_metadata=None, constraint=None):
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
        box_bounds : tuple (optional)
            if specified, then limit the query to the specified quadrilateral
            spatial range.  box_bounds = (RAmin, RAmax, DECmin, DECmax),
            measured in degrees
        constraint : string (optional)
            if specified, the predicate is added to the query

        Returns
        -------
        result : list or iterator
            If chunksize is not specified, then result is a list of all
            items which match the specified query.  If chunksize is specified,
            then result is an iterator over lists of the given size.
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
        
	#KSK: Should we do some checking to make sure circ_bounds and 
	#box_bounds are not both set, or if they are that they are 
	#self consistent
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
        if chunksize is None:
            return self._postprocess_results(exec_query.fetchall())
        else:
            return ChunkIterator(exec_query, chunksize)

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
    #Only realy need to specify the mappings that are not 
    #"x as x" but most need to be mapped anyway.
    columns = [('galtileid', None, numpy.int64),
            ('galid', None, str, 30),
	    ('componentra','bra'),
	    ('componentdec', 'bdec'),
            #This is actually a problem with the stored procedure.  We need to be able to map columns other than
            #just ra/dec to raJ2000/decJ2000.  This gets important when we start perturbing the three galaxy components
            ('raJ2000', 'ra'),
            ('decJ2000', 'dec'),
            ('magNorm', 'magnorm_bulge'),
            ('sedFilename', 'sedname_bulge', unicode, 40),
            ('majorAxis', 'a_b'),
            ('minorAxis', 'b_b'),
            ('positionAngle', 'pa_bulge'),
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
    #Only realy need to specify the mappings that are not 
    #"x as x" but most need to be mapped anyway.
    columns = [('galtileid', None, numpy.int64),
            ('galid', None, str, 30),
	    ('componentra','dra'),
	    ('componentdec', 'ddec'),
            #This is actually a problem with the stored procedure.  We need to be able to map columns other than
            #just ra/dec to raJ2000/decJ2000.  This gets important when we start perturbing the three galaxy components
            ('raJ2000', 'ra'),
            ('decJ2000', 'dec'),
            ('magNorm', 'magnorm_disk'),
            ('sedFilename', 'sedname_disk', unicode, 40),
            ('majorAxis', 'a_d'),
            ('minorAxis', 'b_d'),
            ('positionAngle', 'pa_disk'),
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
    #The tuples must be at least len 2
    columns = [('galtileid', None, numpy.int64),
            ('galid', None, str, 30),
	    ('componentra','agnra'),
	    ('componentdec', 'agndec'),
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

if __name__ == '__main__':
    star = DBObject.from_objid('msstars')
    galaxy = DBObject.from_objid('galaxyBase')
    galaxyTiled = DBObject.from_objid('galaxyTiled')
    galaxyBulge = DBObject.from_objid('galaxyBulge')
    galaxyDisk = DBObject.from_objid('galaxyDisk')
    galaxyAgn = DBObject.from_objid('galaxyAgn')

    objects = [star, galaxy, galaxyTiled, galaxyBulge, galaxyDisk, galaxyAgn]

    constraints = ["rmag < 21.", "gr_total_rest > 0.8", "r_ab < 20.", "mass_bulge > 1.", 
                   "DiskLSSTg < 20.", "t0_agn > 300."]
    obs_metadata = ObservationMetaData(circ_bounds=dict(ra=2.0,
                                                        dec=5.0,
                                                        radius=0.01))
    obs_metadata_gal = ObservationMetaData(circ_bounds=dict(ra=0.0,
                                                            dec=0.0,
                                                            radius=0.01))
    metadataList = [obs_metadata, obs_metadata_gal, obs_metadata, 
		    obs_metadata, obs_metadata, obs_metadata]
    for object, constraint, md in zip(objects, constraints, metadataList):
        result = object.query_columns(obs_metadata=md, constraint=constraint)
        print "Length of returned result set of %s is: %i"%(object.objid, len(result))
