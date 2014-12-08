import warnings
import numpy
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData, ChunkIterator
from sqlalchemy.sql import select, func, column

__all__ = ["StarBase", "StarObj", "MsStarObj", "WdStarObj", "RRLyStarObj",
           "BhbStarObj", "EbStarObj", "CepheidStarObj", "EasterEggStarObj",
           "DwarfGalStarObj"]

class StarBase(CatalogDBObject):
    objid = 'starbase'

    #: This is the default address.  Simply change this in the class definition for other
    #: endpoints.
    dbAddress = "mssql+pymssql://LSST-2:L$$TUser@fatboy.npl.washington.edu:1433/LSST"

    tableid = None
    idColKey = 'id'
    raColName = 'ra'
    decColName = 'decl'
    objectTypeId = -1
    #Don't run test on base class
    doRunTest = False
    #default observation metadata
    testObservationMetaData = ObservationMetaData(boundType='circle', unrefractedRA=210.0, unrefractedDec=-30.0,
                                                  boundLength=0.3, mjd=52000., bandpassName='g',m5=22.0)
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

        if obs_metadata is not None and obs_metadata.bounds is not None:
            if obs_metadata.bounds.boundType == 'circle':
                regionStr = 'REGION CIRCLE J2000 %f %f %f'%(obs_metadata.bounds.RAdeg,
                                                            obs_metadata.bounds.DECdeg,
                                                            60.*obs_metadata.bounds.radiusdeg)
            elif obs_metadata.bounds.boundType == 'box':
                regionStr = 'REGION RECT J2000 %f %f %f %f'%(obs_metadata.bounds.RAminDeg,
                                                             obs_metadata.bounds.DECminDeg,
                                                             obs_metadata.bounds.RAmaxDeg,
                                                             obs_metadata.bounds.DECmaxDeg)
            else:
                raise RuntimeError("StarBase does not know about boundType %s "
                                   % obs_metadata.bounds.boundType)
        else:
            regionStr = 'REGION CIRCLE J2000 180. 0. 10800.'
            warnings.warn("Searching over entire sky "
                          "since no bounds specified. "
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
    testObservationMetaData = ObservationMetaData(boundType = 'circle', unrefractedRA=210.0, unrefractedDec=-30.0,
                                                  boundLength=0.1, mjd=52000., bandpassName='g', m5=22.0)
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
    testObservationMetaData = ObservationMetaData(boundType = 'circle', unrefractedRA=210.0, unrefractedDec=-30.0,
                                                  boundLength=0.1, mjd=52000., bandpassName='g',m5=22.0)
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
    testObservationMetaData = ObservationMetaData(mjd=52000., bandpassName='g', m5=22.0)
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
    testObservationMetaData = ObservationMetaData(mjd=52000., bandpassName='g', m5=22.0)
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
    testObservationMetaData = ObservationMetaData(mjd=52000., bandpassName='g', m5=22.0)
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
    testObservationMetaData = ObservationMetaData(boundType='circle', unrefractedRA=2.0, unrefractedDec=-0.3,
                                                  boundLength=0.1, mjd=52000., bandpassName='g', m5=22.0)
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

