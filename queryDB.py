#!/usr/bin/env python
from dbModel import *
from catalogDbMap import catalogDbMap
from catalogDbMap import physicalObjectMap
import os
import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import CatalogDescription
from lsst.sims.catalogs.measures.instance import Metadata
from lsst.sims.catalogs.measures.instance import CatalogType


class queryDB(object):
  def __init__(self, objtype = 'STARS', filetypes=CatalogType.TRIM, chunksize=100000):
    setup_all()
    self._start=0
    self.filetypes = filetypes
    self.objtype = objtype
    self.chunksize=chunksize
    self.cdm = catalogDbMap()
    self.pom = physicalObjectMap()
    self.queries = None
    self.ptype = 'POINT'

  def getNextChunk(self):
    try:
      result = []
      for q in self.queries:
        if len(result) == self.chunksize:
            break
        else:
            result += q.fetchmany(self.chunksize - len(result))  
      if len(result) == 0:
        return None
      else:
        return self.makeCatalogFromQuery(result, self.ptype)
    except Exception, e:
      print "Exception of type: %s"%(type(e))
      raise Exception(e)


  def getInstanceCatalogById(self, id, radiusdeg=2.1, opsim="OPSIM361", add_columns=()):
    metadataMap = self.cdm.objectTypes[opsim]
    os =\
        eval("%s.query.filter(\"obshistid=%s\").first()"%(self.pom.objectMap[opsim][0]['table'], id))
    objtype = self.objtype
    self.metadata = Metadata()
    for k in metadataMap.keys():
      self.metadata.addMetadata(k,os.__dict__[metadataMap[k]],"")
    om = self.pom.objectMap[objtype]
    queries = []
    if objtype == 'STARS':
      self.ptype = 'POINT'
      for map in om:
        query = session.query(eval("%s.id"%(map['table'])))
        for k in self.cdm.objectTypes[map['ptype']].keys():
          if k == 'id':
	    continue
          elif k == 'magNorm':
            col = expression.literal_column(self.cdm.objectTypes[map['ptype']][k]%(os.expmjd, os.filter)).label(k)
            query = query.add_column(col)
	  else:
            col = expression.literal_column(self.cdm.objectTypes[map['ptype']][k]).label(k)
            query = query.add_column(col)
        query = query.filter("point @ scircle \'<(%fd,%fd),%fd>\'"%(os.fieldradeg, os.fielddecdeg, radiusdeg))
        if map['constraint'] is not None:
          query = query.filter(map['constraint'])
        queries.append(session.execute(query))
      self.queries = queries
      result = []
      for q in self.queries:
        if len(result) == self.chunksize:
            break
        else:
            result += q.fetchmany(self.chunksize - len(result))  

      if len(result) == 0:
        return None
      else:
        return self.makeCatalogFromQuery(result, self.ptype)
    if objtype == 'STARRAW':
      self.ptype = 'STARRAW'
      for map in om:
        query = session.query(eval("%s.id"%(map['table'])))
        for k in self.cdm.objectTypes[map['ptype']].keys():
          if k == 'id':
	    continue
          elif k == 'magNorm':
            col = expression.literal_column(self.cdm.objectTypes[map['ptype']][k]%(os.expmjd, os.filter)).label(k)
            query = query.add_column(col)
	  else:
            col = expression.literal_column(self.cdm.objectTypes[map['ptype']][k]).label(k)
            query = query.add_column(col)
        query = query.filter("point @ scircle \'<(%fd,%fd),%fd>\'"%(os.fieldradeg, os.fielddecdeg, radiusdeg))
        if map['constraint'] is not None:
          query = query.filter(map['constraint'])
        queries.append(session.execute(query))
      self.queries = queries
      result = []
      for q in self.queries:
        if len(result) == self.chunksize:
            break
        else:
            result += q.fetchmany(self.chunksize - len(result))  

      if len(result) == 0:
        return None
      else:
        print self.ptype
        return self.makeCatalogFromQuery(result, self.ptype)
    elif objtype == 'WDSTARS':
      self.ptype = 'POINT'
      for map in om:
        query = session.query(eval("%s.id"%(map['table'])))
        for k in self.cdm.objectTypes[map['ptype']].keys():
          if k == 'id':
	    continue
          elif k == 'magNorm':
            col = expression.literal_column(self.cdm.objectTypes[map['ptype']][k]%(os.expmjd, os.filter)).label(k)
            query = query.add_column(col)
	  else:
            col = expression.literal_column(self.cdm.objectTypes[map['ptype']][k]).label(k)
            query = query.add_column(col)
        query = query.filter("point @ scircle \'<(%fd,%fd),%fd>\'"%(os.fieldradeg, os.fielddecdeg, radiusdeg))
        if map['constraint'] is not None:
          query = query.filter(map['constraint'])
        queries.append(session.execute(query))
      self.queries = queries
      result = []
      for q in self.queries:
        if len(result) == self.chunksize:
            break
        else:
            result += q.fetchmany(self.chunksize - len(result))  

      if len(result) == 0:
        return None
      else:
        return self.makeCatalogFromQuery(result, self.ptype)
    elif objtype == 'ssm':
      '''Need to do query and then do the ephemeris calculation
      '''
      pass
    elif objtype == 'GALRAW':
      '''Remember that the galaxy has several components: Bulge, Disk, AGN
      '''
      self.ptype = 'GALRAW'
      tiles = self.getTiles(os.fieldradeg, os.fielddecdeg, radiusdeg)
      for tile in tiles:
        for map in om:
          query = session.query(eval("%s.id"%(map['table'])))
          if map['constraint'] is not None:
            query = query.filter(map['constraint'])
          for k in self.cdm.objectTypes[map['ptype']].keys():
            if k == 'id':
              continue
            else:
              col = expression.literal_column(self.cdm.objectTypes[map['ptype']][k]).label(k)
              query = query.add_column(col)
          query = query.filter("point @ scircle '%s' and (point + strans(0,\
              %f*PI()/180.,  %f*PI()/180., 'XYZ')) @ spoly\
              '%s'"%(tile['circ'],-tile['decmid'],tile['ramid'],tile['bbox']))
          queries.append(session.execute(query))
      self.queries = queries
      result = []
      for q in self.queries:
        if len(result) == self.chunksize:
            break
        else:
            result += q.fetchmany(self.chunksize - len(result))  

      if len(result) == 0:
        return None
      else:
        return self.makeCatalogFromQuery(result,self.ptype)
    elif objtype == 'GALAXY':
      '''Remember that the galaxy has several components: Bulge, Disk, AGN
      '''
      self.ptype = 'SERSIC2D'
      tiles = self.getTiles(os.fieldradeg, os.fielddecdeg, radiusdeg)
      for tile in tiles:
        for map in om:
          query = session.query(eval("%s.id"%(map['table'])))
          if map['constraint'] is not None:
            query = query.filter(map['constraint'])
          for k in self.cdm.objectTypes[map['ptype']].keys():
            if k == 'id':
              continue
            else:
              col = expression.literal_column(self.cdm.objectTypes[map['ptype']][k]).label(k)
              query = query.add_column(col)
          query = query.filter("point @ scircle '%s' and (point + strans(0,\
              %f*PI()/180.,  %f*PI()/180., 'XYZ')) @ spoly\
              '%s'"%(tile['circ'],-tile['decmid'],tile['ramid'],tile['bbox']))
          queries.append(session.execute(query))
      self.queries = queries
      result = []
      for q in self.queries:
        if len(result) == self.chunksize:
            break
        else:
            result += q.fetchmany(self.chunksize - len(result))  

      if len(result) == 0:
        return None
      else:
        return self.makeCatalogFromQuery(result,self.ptype)
    else:
      raise Exception('getInstanceCatalogById', 'Did not give valid object type')

  def makeCatalogFromQuery(self, result, ptype):
    if os.environ.has_key("CATALOG_DESCRIPTION_PATH"):
      catalogDescriptionPath = os.environ["CATALOG_DESCRIPTION_PATH"]
    else:
      raise Exception("Environment variable CATALOG_DESCRIPTION_PATH not set to location of the catalog description files")
    nic = InstanceCatalog()
    nic.catalogDescription = CatalogDescription.CatalogDescription(
                   catalogDescriptionPath+"requiredMetadata.dat",
                   catalogDescriptionPath+"requiredSchemaFields.dat",
                   catalogDescriptionPath+"requiredDerivedFields.dat",
                   catalogDescriptionPath+"outputFormat.dat")
    nic.metadata = self.metadata
    nic.metadata.catalogDescription =  nic.catalogDescription
    nic.catalogType = self.filetypes
    data = {}
    for k in self.cdm.objectTypes[ptype].keys():
      data[k] = []
    for s in result:
      for k in self.cdm.objectTypes[ptype].keys():
        exec("data[k].append(s.%s)"%(k))
    for k in self.cdm.objectTypes[ptype].keys():
      arr = numpy.asarray(data[k])
      nic.addColumn(arr, k)
    nic.metadata = self.metadata
    if nic == None:
        raise RuntimeError, '*** nic is None'
    if nic.metadata == None:
        raise RuntimeError, '*** nic.metadata is None'
    if len(nic.dataArray) < 1:
        raise RuntimeError, '*** nic.dataArray has len < 1'

    return nic

  def getTiles(self, raDeg, decDeg, radiusDeg):
    sq = session.query(Tiles.ramid, Tiles.decmid, Tiles.bbox).filter("scircle\
            '<(%fd, %fd), %fd>' && bbox"%(raDeg, decDeg,
            radiusDeg)).subquery().alias(name="mytiles")
    col = expression.literal_column("scircle '<(%fd, %fd), %fd>'+\
            strans(-%s.ramid*PI()/180., %s.decmid*PI()/180.,0.0,\
            'ZYX')"%(raDeg,
            decDeg, radiusDeg,"mytiles","mytiles")).label("circ")
    query = session.query(col,sq.c.ramid,sq.c.decmid,sq.c.bbox)
    result = query.all()
    tiles = []
    for r in result:
      tiles.append({'circ':r.circ, 'ramid':r.ramid, 'decmid':r.decmid,
            'bbox':r.bbox})
    return tiles
