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
    idColName = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 9
    #There is no spatial model available for coadded galaxies 
    spatialModel = None
    #Only realy need to specify the mappings that are not 
    #"x as x" but most need to be mapped anyway.
    columns = [('galid', str, 30),
            ('raJ2000', float),
            ('decJ2000', float),
            ('raJ2000Bulge', float),
            ('decJ2000Bulge', float),
            ('raJ2000Disk', float),
            ('decJ2000Disk', float),
            ('raJ2000Agn', float),
            ('decJ2000Agn', float),
            ('magNormBulge', float),
            ('magNormDisk', float),
            ('magNormAgn', float),
            ('sedFilenameBulge', unicode, 40),
            ('sedFilenameDisk', unicode, 40),
            ('sedFilenameAgn', unicode, 40),
            ('majorAxisBulge', float),
            ('minorAxisBulge', float),
            ('positionAngleBulge', float),
            ('sindexBulge', int),
            ('majorAxisDisk', float),
            ('minorAxisDisk', float),
            ('positionAngleDisk', float),
            ('sindexDisk', int),
            ('internalExtinctionModelBulge', str, 3),
            ('internalAvBulge', float),
            ('internalRvBulge', float),
            ('internalExtinctionModelDisk', str, 3),
            ('internalAvDisk', float),
            ('internalRvDisk', float),
            ('redshift', float),
            ('radialVelocity', float),
            ('lsst_u', float),
            ('lsst_g', float),
            ('lsst_r', float),
            ('lsst_i', float),
            ('lsst_z', float),
            ('lsst_y', float)]
    column_map = {    
            'galid':'galid',
            'raJ2000':'ra*PI()/180.',
            'decJ2000':'dec*PI()/180.',
            'raJ2000Bulge':'bra*PI()/180.',
            'decJ2000Bulge':'bdec*PI()/180.',
            'raJ2000Disk':'dra*PI()/180.',
            'decJ2000Disk':'ddec*PI()/180.',
            'raJ2000Agn':'agnra*PI()/180.',
            'decJ2000Agn':'agndec*PI()/180.',
            'magNormBulge':'magnorm_bulge',
            'magNormDisk':'magnorm_disk',
            'magNormAgn':'magnorm_agn',
            'sedFilenameBulge':'sedname_bulge',
            'sedFilenameDisk':'sedname_disk',
            'sedFilenameAgn':'sedname_agn',
            'majorAxisBulge':'a_b',
            'minorAxisBulge':'b_b',
            'positionAngleBulge':'pa_bulge',
            'sindexBulge':'bulge_n',
            'majorAxisDisk':'a_d',
            'minorAxisDisk':'b_d',
            'positionAngleDisk':'pa_disk',
            'sindexDisk':'disk_n',
            'internalExtinctionModelBulge':'ext_model_b',
            'internalAvBulge':'av_b',
            'internalRvBulge':'rv_b',
            'internalExtinctionModelDisk':'ext_model_d',
            'internalAvDisk':'av_d',
            'internalRvDisk':'rv_d',
            'redshift':'redshift',
            'radialVelocity':'rad_vel',
            'lsst_u':'u_ab',
            'lsst_g':'g_ab',
            'lsst_r':'r_ab',
            'lsst_i':'i_ab',
            'lsst_z':'z_ab',
            'lsst_y':'y_ab'}


class GalaxyTileObj(DBObject):
    objid = 'galaxyTiled'
    #This is the base table for the galaxies
    tableid = 'galaxy'
    idColName = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 9
    #There is no spatial model available for coadded galaxies 
    spatialModel = None
    #Only realy need to specify the mappings that are not 
    #"x as x" but most need to be mapped anyway.
    columns = [('galtileid', numpy.int64),
            ('galid', str, 30),
            ('raJ2000', float),
            ('decJ2000', float),
            ('raJ2000Bulge', float),
            ('decJ2000Bulge', float),
            ('raJ2000Disk', float),
            ('decJ2000Disk', float),
            ('raJ2000Agn', float),
            ('decJ2000Agn', float),
            ('magNormBulge', float),
            ('magNormDisk', float),
            ('magNormAgn', float),
            ('sedFilenameBulge', unicode, 40),
            ('sedFilenameDisk', unicode, 40),
            ('sedFilenameAgn', unicode, 40),
            ('majorAxisBulge', float),
            ('minorAxisBulge', float),
            ('positionAngleBulge', float),
            ('sindexBulge', int),
            ('majorAxisDisk', float),
            ('minorAxisDisk', float),
            ('positionAngleDisk', float),
            ('sindexDisk', int),
            ('internalExtinctionModelBulge', str, 3),
            ('internalAvBulge', float),
            ('internalRvBulge', float),
            ('internalExtinctionModelDisk', str, 3),
            ('internalAvDisk', float),
            ('internalRvDisk', float),
            ('redshift', float),
            ('radialVelocity', float),
            ('lsst_u', float),
            ('lsst_g', float),
            ('lsst_r', float),
            ('lsst_i', float),
            ('lsst_z', float),
            ('lsst_y', float)]
    column_map = {    
            'galtileid':'galtileid',
            'galid':'galid',
            'raJ2000':'ra',
            'decJ2000':'dec',
            'raJ2000Bulge':'bra*PI()/180.',
            'decJ2000Bulge':'bdec*PI()/180.',
            'raJ2000Disk':'dra*PI()/180.',
            'decJ2000Disk':'ddec*PI()/180.',
            'raJ2000Agn':'agnra*PI()/180.',
            'decJ2000Agn':'agndec*PI()/180.',
            'magNormBulge':'magnorm_bulge',
            'magNormDisk':'magnorm_disk',
            'magNormAgn':'magnorm_agn',
            'sedFilenameBulge':'sedname_bulge',
            'sedFilenameDisk':'sedname_disk',
            'sedFilenameAgn':'sedname_agn',
            'majorAxisBulge':'a_b',
            'minorAxisBulge':'b_b',
            'positionAngleBulge':'pa_bulge',
            'sindexBulge':'bulge_n',
            'majorAxisDisk':'a_d',
            'minorAxisDisk':'b_d',
            'positionAngleDisk':'pa_disk',
            'sindexDisk':'disk_n',
            'internalExtinctionModelBulge':'ext_model_b',
            'internalAvBulge':'av_b',
            'internalRvBulge':'rv_b',
            'internalExtinctionModelDisk':'ext_model_d',
            'internalAvDisk':'av_d',
            'internalRvDisk':'rv_d',
            'redshift':'redshift',
            'radialVelocity':'rad_vel',
            'lsst_u':'u_ab',
            'lsst_g':'g_ab',
            'lsst_r':'r_ab',
            'lsst_i':'i_ab',
            'lsst_z':'z_ab',
            'lsst_y':'y_ab'}

    def _get_column_query(self, colnames=None):
        raise NotImplementedError("We are calling a stored procedure so "
                                  "no need to loop over columns")

    def _postprocess_results(self, results):
        #Here we will need to do the post processing of the results to return radians in the raJ2000 and decJ2000 fields.
        #This is a result of the fact that the stored procedure does not allow the ra/dec columns to be modified in the return
        #result set.
        retresults = numpy.zeros((len(results),),dtype=self.dtype)
        for i, result in enumerate(results):
            for k in self.requirements.keys():
                retresults[i][k] = result[k]
	retresults['raJ2000'] = numpy.radians(retresults['raJ2000'])
	retresults['decJ2000'] = numpy.radians(retresults['decJ2000'])
        return retresults

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
            colnames = [el[0] for el in self.columns]
        
        #We know that galtileid comes back with the query, but we don't want 
        #to add it to the query since it's generated on the fly.
        while 'galtileid' in colnames:
            colnames.remove('galtileid')
        mappedcolnames = ["%s as %s"%(self.column_map[x], x) for x in colnames]
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


        query = "EXECUTE [LSST].[dbo].[GalaxySearchSpecColsConstraint]\
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
    idColName = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 1
    spatialModel = 'SERSIC2D'
    #Only realy need to specify the mappings that are not 
    #"x as x" but most need to be mapped anyway.
    columns = [('galtileid', numpy.int64),
            ('galid', str, 30),
	    ('componentra', float),
	    ('componentdec', float),
            ('raJ2000', float),
            ('decJ2000', float),
            ('magNorm', float),
            ('sedFilename', unicode, 40),
            ('majorAxis', float),
            ('minorAxis', float),
            ('positionAngle', float),
            ('sindex', int),
            ('internalExtinctionModel', str, 3),
            ('internalAv', float),
            ('internalRv', float),
            ('redshift', float),
            ('radialVelocity', float),
            ('lsst_u', float),
            ('lsst_g', float),
            ('lsst_r', float),
            ('lsst_i', float),
            ('lsst_z', float),
            ('lsst_y', float)]
    column_map = {    
            'galtileid':'galtileid',
            'galid':'galid',
            'componentra':'bra',
            'componentdec':'bdec',
            #This is actually a problem with the stored procedure.  We need to be able to map columns other than
            #just ra/dec to raJ2000/decJ2000.  This gets important when we start perturbing the three galaxy components
            'raJ2000':'ra',
            'decJ2000':'dec',
            'magNorm':'magnorm_bulge',
            'sedFilename':'sedname_bulge',
            'majorAxis':'a_b',
            'minorAxis':'b_b',
            'positionAngle':'pa_bulge',
            'sindex':'bulge_n',
            'internalExtinctionModel':'ext_model_b',
            'internalAv':'av_b',
            'internalRv':'rv_b',
            'redshift':'redshift',
            'radialVelocity':'rad_vel',
            'lsst_u':'u_ab',
            'lsst_g':'g_ab',
            'lsst_r':'r_ab',
            'lsst_i':'i_ab',
            'lsst_z':'z_ab',
            'lsst_y':'y_ab'}

class GalaxyDiskObj(GalaxyTileObj):
    objid = 'galaxyDisk'
    #This is the base table for the galaxies
    #with a disk component
    tableid = 'galaxy'
    idColName = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 2
    spatialModel = 'SERSIC2D'
    #Only realy need to specify the mappings that are not 
    #"x as x" but most need to be mapped anyway.
    columns = [('galtileid', numpy.int64),
            ('galid', str, 30),
	    ('componentra', float),
	    ('componentdec', float),
            ('raJ2000', float),
            ('decJ2000', float),
            ('magNorm', float),
            ('sedFilename', unicode, 40),
            ('majorAxis', float),
            ('minorAxis', float),
            ('positionAngle', float),
            ('sindex', int),
            ('internalExtinctionModel', str, 3),
            ('internalAv', float),
            ('internalRv', float),
            ('redshift', float),
            ('radialVelocity', float),
            ('lsst_u', float),
            ('lsst_g', float),
            ('lsst_r', float),
            ('lsst_i', float),
            ('lsst_z', float),
            ('lsst_y', float)]
    column_map = {    
            'galtileid':'galtileid',
            'galid':'galid',
            'componentra':'dra',
            'componentdec':'ddec',
            #This is actually a problem with the stored procedure.  We need to be able to map columns other than
            #just ra/dec to raJ2000/decJ2000.  This gets important when we start perturbing the three galaxy components
            'raJ2000':'ra',
            'decJ2000':'dec',
            'magNorm':'magnorm_disk',
            'sedFilename':'sedname_disk',
            'majorAxis':'a_d',
            'minorAxis':'b_d',
            'positionAngle':'pa_disk',
            'sindex':'disk_n',
            'internalExtinctionModel':'ext_model_d',
            'internalAv':'av_d',
            'internalRv':'rv_d',
            'redshift':'redshift',
            'radialVelocity':'rad_vel',
            'lsst_u':'u_ab',
            'lsst_g':'g_ab',
            'lsst_r':'r_ab',
            'lsst_i':'i_ab',
            'lsst_z':'z_ab',
            'lsst_y':'y_ab'}

class GalaxyAgnObj(GalaxyTileObj):
    objid = 'galaxyAgn'
    #This is the base table for the galaxies
    #with an agn component
    tableid = 'galaxy_agn'
    idColName = 'galid'
    raColName = 'ra'
    decColName = 'dec'
    appendint = 3
    spatialModel = 'ZPOINT'
    #Only realy need to specify the mappings that are not 
    #"x as x" but most need to be mapped anyway.
    columns = [('galtileid', numpy.int64),
            ('galid', str, 30),
	    ('componentra', float),
	    ('componentdec', float),
            ('raJ2000', float),
            ('decJ2000', float),
            ('magNorm', float),
            ('sedFilename', unicode, 40),
            ('redshift', float),
            ('radialVelocity', float),
	    ('variabilityParameters', str, 256),
            ('lsst_u', float),
            ('lsst_g', float),
            ('lsst_r', float),
            ('lsst_i', float),
            ('lsst_z', float),
            ('lsst_y', float)]
    column_map = {    
            'galtileid':'galtileid',
            'galid':'galid',
            'componentra':'agnra',
            'componentdec':'agndec',
            #This is actually a problem with the stored procedure.  We need to be able to map columns other than
            #just ra/dec to raJ2000/decJ2000.  This gets important when we start perturbing the three galaxy components
            'raJ2000':'ra',
            'decJ2000':'dec',
            'magNorm':'magnorm_agn',
            'sedFilename':'sedname_agn',
            'redshift':'redshift',
            'radialVelocity':'rad_vel',
            'variabilityParameters':'varParamStr',
            'lsst_u':'u_ab',
            'lsst_g':'g_ab',
            'lsst_r':'r_ab',
            'lsst_i':'i_ab',
            'lsst_z':'z_ab',
            'lsst_y':'y_ab'}

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
