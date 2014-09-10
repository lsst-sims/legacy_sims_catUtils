import warnings
from lsst.sims.catalogs.generation.db import DBObject, ObservationMetaData, ChunkIterator
from sqlalchemy.sql import select, func, column

class StarBase(DBObject):
    objid = 'starbase'
    tableid = None
    idColKey = 'id'
    raColName = 'ra'
    decColName = 'decl'
    objectTypeId = -1
    #Don't run test on base class
    doRunTest = False
    #default observation metadata
    testObservationMetaData = ObservationMetaData(circ_bounds=dict(ra=210., dec=-30., radius=0.3),
                                                  mjd=52000., bandpassName='g')
    dbDefaultValues = {'varsimobjid':-1, 'runid':-1, 'ismultiple':-1, 'run':-1,
                       'runobjid':-1}
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]


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

        circ_bounds = None
        box_bounds = None
        if obs_metadata is not None:
            if obs_metadata.circ_bounds is not None:
                circ_bounds = obs_metadata.circ_bounds
            if obs_metadata.box_bounds is not None:
                box_bounds = obs_metadata.box_bounds

        if circ_bounds is not None:
            regionStr = 'REGION CIRCLE J2000 %f %f %f'%(circ_bounds['ra'], circ_bounds['dec'],
                                                        60.*circ_bounds['radius'])


        elif box_bounds is not None:
            regionStr = 'REGION RECT J2000 %f %f %f %f'%(box_bounds['ra_min'], box_bounds['dec_min'],
                                                         box_bounds['ra_max'],box_bounds['dec_max'])


        else:
            regionStr = None
            warnings.warn("Searching over entire sky "
                          "since no circ_bounds specified. "
                          "This could be a very bad idea "
                          "if the database is large")

        if obs_metadata is not None and regionStr is not None:
            #add spatial constraints to query.

            #Hint sql engine to seek on htmid
            if not self.tableid.endswith('forceseek'):
                query = query.with_hint(self.table, ' WITH(FORCESEEK)', 'mssql')

            #aliased subquery for htmid ranges covering the search region
            htmid_range_alias = select([column('htmidstart'), column('htmidend')]).\
            select_from(func.fHtmCoverBinaryAdvanced(
                    func.sph.fSimplifyString(regionStr))).alias()

            #Range join on htmid ranges
            query = query.join(htmid_range_alias,
                       self.table.c['htmID'].between(htmid_range_alias.c.htmidstart,
                                                     htmid_range_alias.c.htmidend)
                       )
            query = query.filter(func.sph.fRegionContainsXYZ(func.sph.fSimplifyString(regionStr),
                     self.table.c.cx, self.table.c.cy, self.table.c.cz) == 1)


        if constraint is not None:
            query = query.filter(constraint)

        return ChunkIterator(self, query, chunk_size)



class StarObj(StarBase):
    objid = 'allstars'
    tableid = 'starsALL_forceseek'
    objectTypeId = 4
    doRunTest = True
    testObservationMetaData = ObservationMetaData(circ_bounds=dict(ra=210., dec=-30., radius=0.1),
                                                  mjd=52000., bandpassName='g')
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

class MsStarObj(StarBase):
    objid = 'msstars'
    tableid = 'starsMSRGB_forceseek'
    objectTypeId = 5
    doRunTest = True
    testObservationMetaData = ObservationMetaData(circ_bounds=dict(ra=210., dec=-30., radius=0.1),
                                                  mjd=52000., bandpassName='g')
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

class WdStarObj(StarBase):
    objid = 'wdstars'
    tableid = 'starsWD'
    objectTypeId = 6
    doRunTest = True
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

class RRLyStarObj(StarBase):
    objid = 'rrlystars'
    tableid = 'starsRRLy'
    objectTypeId = 7
    doRunTest = True
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

class BhbStarObj(StarBase):
    objid = 'bhbstars'
    tableid = 'starsBHB'
    objectTypeId = 8
    doRunTest = True
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

class EbStarObj(StarBase):
    objid = 'ebstars'
    tableid = 'ebstars'
    objectTypeId = 9
    doRunTest = True
    testObservationMetaData = ObservationMetaData(circ_bounds=None,
                                                  mjd=52000., bandpassName='g')
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

class CepheidStarObj(StarBase):
    objid = 'cepheidstars'
    tableid = 'cepheidstars'
    objectTypeId = 10
    doRunTest = True
    testObservationMetaData = ObservationMetaData(circ_bounds=None,
                                                  mjd=52000., bandpassName='g')
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

class EasterEggStarObj(StarBase):
    objid = 'eastereggstars'
    tableid = 'AstromEasterEggs'
    objectTypeId = 11
    doRunTest = True
    testObservationMetaData = ObservationMetaData(circ_bounds=None,
                                                  mjd=52000., bandpassName='g')
    dbDefaultValues = StarBase.dbDefaultValues
    dbDefaultValues['sedid']=-1
    dbDefaultValues['especid']=-1
    dbDefaultValues['pop']=-1
    dbDefaultValues['type']=-1
    dbDefaultValues['isvar']=False
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

class DwarfGalStarObj(StarBase):
    objid = 'dwarfgalstars'
    tableid = 'dwarfGalaxies'
    objectTypeId = 12
    doRunTest = True
    testObservationMetaData = ObservationMetaData(circ_bounds={'ra':2., 'dec':-0.3, 'radius':0.1},
                                                  mjd=52000., bandpassName='g')
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

