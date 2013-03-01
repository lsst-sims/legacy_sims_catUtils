import warnings
from dbConnection import ChunkIterator, DBObject
from sqlalchemy import Table, Column, BigInteger, MetaData

class GalaxyBulgeObj(DBObject):
    objid = 'galaxy_bulge'
    #This is the base table for the galaxies
    #with a bulge component
    tableid = 'galaxy_bulge'
    idColName = 'galid'
    appendint = 1
    spatialModel = 'SERSIC2D'
    column_map = {    
            'galtileid':'galtileid',
            'galid':'galid',
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
    columns = column_map.keys()

    def _get_column_query(self, colnames=None):
        raise NotImplementedError("We are calling a stored procedure so "
                                  "no need to loop over columns")

    def query_columns(self, colnames=None, chunksize=None,
                      circ_bounds=None, box_bounds=None, constraint=None):
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
            colnames = self.columns
        
        #We know that galtileid comes back with the query, but we don't want 
        #to add it to the query since it's generated on the fly.
        while 'galtileid' in colnames:
            colnames.remove('galtileid')

        mappedcolnames = ["%s as %s"%(self.column_map[x], x) for x in colnames]
        mappedcolnames = ",".join(mappedcolnames)

        if circ_bounds is not None:
            RA, DEC, radius = circ_bounds
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

        if box_bounds is not None:
            raise NotImplementedError("There is no way to specify a box "
                                      "for the galaxy stored procedure")

        exec_query = self.session.execute(query)

        if chunksize is None:
            return exec_query.fetchall()
        else:
            return ChunkIterator(exec_query, chunksize)

class GalaxyDiskObj(GalaxyBulgeObj):
    #Can't use super().__init__() since that calles _get_table()
    objid = 'galaxy_disk'
    #This is the base table for the galaxies
    #with a disk component
    tableid = 'galaxy'
    appendint = 2
    spatialModel = 'SERSIC2D'
    column_map = {    
            'galtileid':'galtileid',
            'galid':'galid',
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
    columns = column_map.keys()

class GalaxyAgnObj(GalaxyBulgeObj):
    objid = 'galaxy_agn'
    #This is the base table for the galaxies
    #with a bulge component
    tableid = 'galaxy_agn'
    idColName = 'galid'
    appendint = 3
    spatialModel = 'ZPOINT'
    column_map = {    
            'galtileid':'galtileid',
            'galid':'galid',
            'raJ2000':'ra',
            'decJ2000':'dec',
            'magNorm':'magnorm_agn',
            'sedFilename':'sedname_agn',
            'redshift':'redshift',
            'internalExtinctionModel':"''none''",
            'radialVelocity':'rad_vel',
            'variabilityParameters':'varParamStr',
            'lsst_u':'u_ab',
            'lsst_g':'g_ab',
            'lsst_r':'r_ab',
            'lsst_i':'i_ab',
            'lsst_z':'z_ab',
            'lsst_y':'y_ab'}
    columns = column_map.keys()

if __name__ == '__main__':
    #star = StarObj()
    star = DBObject.from_objid('msstars')
    galaxyBulge = DBObject.from_objid('galaxy_bulge')
    galaxyDisk = DBObject.from_objid('galaxy_disk')
    galaxyAgn = DBObject.from_objid('galaxy_agn')
    print star.query_columns(circ_bounds=(2.0, 5.0, 1.0))
    print galaxyBulge.query_columns(circ_bounds=(200., -10, 0.1), constraint="sedname_bulge is not NULL")
    print galaxyDisk.query_columns(circ_bounds=(200., -10, 0.1))
    print galaxyAgn.query_columns(circ_bounds=(200., -10, 0.1), constraint="sedname_agn is not NULL")
