#!/usr/bin/env python
from dbMsModel import *
import time
from copy import deepcopy
import re
import os
import math
import numpy
import warnings
from lsst.sims.catalogs.generation.config import ConfigObj
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import CatalogDescription
from lsst.sims.catalogs.measures.instance import Metadata
from lsst.sims.catalogs.measures.astrometry import Bbox
from sqlalchemy import select 
from sqlalchemy.orm import join
from sqlalchemy.orm import aliased
from sqlalchemy.sql.expression import text
from sqlalchemy.sql.expression import select
from sqlalchemy.sql.expression import between


class queryDB(object):
  def __init__(self, objtype = 'STARS', filetypes=('TRIM',), chunksize=100000,
          makeCat=True, dithered=False):
    if os.environ.has_key("CATALOG_DESCRIPTION_PATH"):
      catalogDescriptionPath = os.environ["CATALOG_DESCRIPTION_PATH"]
    else:
      raise Exception("Environment variable CATALOG_DESCRIPTION_PATH not set to location of the catalog description files")
    catalogDescriptionPath = self.getEnvironPath("CATALOG_DESCRIPTION_PATH")
    dbMapConfigFile = catalogDescriptionPath+"requiredFields.dat"
    objConfigFile = catalogDescriptionPath+"objectMap.dat"
    objEnumConfigFile = catalogDescriptionPath+"objectEnum.dat"
    metaConfigFile = catalogDescriptionPath+"requiredMetadata.dat"
    self.catDescription = CatalogDescription(catalogDescriptionPath+"/config.dat")
    self.metadata = Metadata(catalogDescriptionPath+"/config.dat")
    self.filetypes = filetypes
    self.objtype = objtype
    self.chunksize=chunksize
    self.dm = ConfigObj(dbMapConfigFile)
    self.om = ConfigObj(objConfigFile)
    self.em = ConfigObj(objEnumConfigFile)['OBJECTTYPES']
    self.mm = ConfigObj(metaConfigFile)
    self.queries = None
    self.coldesc = None
    self.ptype = 'POINT'
    self.ftype = 'POINT'
    self.component = None
    self.expmjd = None
    self.centradeg = None
    self.centdecdeg = None
    self.radiusdeg = None
    self.filter = None
    self.opsimmeta = None
    self.opsim = ""
    self.curtile = None
    self.makeCat = makeCat
    self.dithered = dithered

  def closeSession(self):
    session.close_all()

  def getNextChunk(self):
    result = self.queries.fetchmany(self.chunksize)
    if len(result) == 0:
      self.closeSession()
      return None
    else:
      if self.makeCat:
        cat = self.makeCatalogFromQuery(result, self.dithered)
      else:
	cat = True
      return cat

  def getOpsimMetadataById(self, id, opsim="OPSIM361"):
    self.opsim = opsim
    omap = self.om['METADATA'][self.opsim]
    cols = {}
    for ft in self.filetypes:
      map = self.mm[ft][self.opsim]
      for k in map:
        if k == omap['idkey']:
          continue
	if cols.has_key(k):
	  print "Replacing expression %s with %s for key %s"%(cols[k], map[k][1], k)
          cols[k] = expression.literal_column(map[k][1]).label(k)
        else:
          cols[k] = expression.literal_column(map[k][1]).label(k) 
    query =\
    session.query(eval("%s.%s.label(\"%s\")"%(omap['table'],self.mm[self.filetypes[0]][self.opsim][omap['idkey']][1],omap['idkey'])))
    query = query.filter("%s=%s"%(self.mm[self.filetypes[0]][self.opsim][omap['idkey']][1], id))
    for k in cols.keys():
      query = query.add_column(cols[k])
    self.opsimmeta = query.first()
    self.centradeg =\
    math.degrees(eval("self.opsimmeta.%s"%(omap['rakey'])))
    self.centdecdeg =\
    math.degrees(eval("self.opsimmeta.%s"%(omap['deckey'])))
    if self.centdecdeg < -90.:
        self.centdecdeg = -self.centdecdeg - 180.
        self.centradeg += 180.
        self.centradeg = self.centradeg%360.
    if self.centdecdeg > 90.:
        self.centdecdeg = -self.centdecdeg + 180.
        self.centradeg += 180.
        self.centradeg = self.centradeg%360.
    self.expmjd =\
    eval("self.opsimmeta.%s"%(omap['expmjdkey']))
    self.filter =\
    eval("self.opsimmeta.%s"%(omap['filterkey']))

  def addSpatial(self, query, map, type, ta):
    if type == "circle":
      sel = select(["htmidstart, htmidend"],
              whereclause="innerflag=0",
              from_obj=["LSST.dbo.fHtmCoverBinaryAdvanced([LSST].sph.fSimplifyString('REGION\
              CIRCLE j2000 %f %f %f' ))"%(self.centradeg, self.centdecdeg,
              self.radiusdeg*60.)])
      sela = sel.alias('c')
      onclause = between(ta.htmID, text("c.htmidstart"), text("c.htmidend"))
      query = query.join((sela, onclause))
      aperture = "[LSST].sph.fSimplifyString('REGION CIRCLE j2000 %f %f %f')"%(self.centradeg, 
                self.centdecdeg, self.radiusdeg*60.)
      filterstr = "%s = 1"%(str(func.LSST.sph.fRegioncontainsXYZ(text(aperture), ta.cx, ta.cy, ta.cz)))
      query = query.filter(filterstr)
      return query
    else:
        return query

  def getUnionQuery(self, filter=None, type="circle"):  
    queries = []
    for omkey in self.om[self.objtype].keys():
      if omkey == "formatas":
        continue
      if omkey == "component":
	self.component = self.om[self.objtype][omkey]
	continue
      appint = int(self.em[omkey]['id'])
      map = self.om[self.objtype][omkey]
      idcolstr = self.dm[self.filetypes[0]][self.component][self.ptype][map['idkey']][1]
      if idcolstr.startswith("%%"):
          idcolstr = idcolstr.lstrip("%%")
          idcolstr = eval(idcolstr)
      ta = eval("aliased(%s, name='star_alias')"%map['table'])
      query = session.query(eval("ta.%s.label(\"%s\")"%(idcolstr, map['idkey'])))
      query = query.add_column(expression.literal_column("%i"%(appint)).label("appendint"))
      query = self.addSpatial(query, map, "circle", ta)
      query = self.addUniqueCols(map, query)
      if filter is not None:
        query = query.filter(filter)
      if not len(map['constraint'].strip()) == 0:
        const = map['constraint']
        if const.startswith("%%"):
          const = const.lstrip("%%")
          const = eval(const)
        query = query.filter(const)
      queries.append(query)
    self.coldesc = query.column_descriptions
    #Workaround for adding a hint
    query = queries[0]
    for i in range(len(queries)-1):
      query=query.union_all(queries[i+1])
    #End workaround
    query = session.execute(query)
    return query

  def addUniqueCols(self, map, query):
    cols = self.getUniqueColNames(map)
    for k in cols.keys():
      cols[k] = expression.literal_column(cols[k]).label(k)
      query = query.add_column(cols[k])
    return query

  def getUniqueColNames(self, map):
    cols = {}
    for ft in self.filetypes:
      dmap = self.dm[ft][self.component][map['ptype']]
      for k in dmap.keys():
        colstr = dmap[k][1]
        if k == map['idkey']:
          continue
        else:
          if cols.has_key(k):
            print "Replacing expression %s with %s for key %s"%(cols[k], colstr, k)
            cols[k] = colstr
          else:
            cols[k] = colstr 
    return cols


  def getInstanceCatalogById(self, id, opsim="OPSIM361", radiusdeg=2.1):
    self.getOpsimMetadataById(id, opsim)
    return self.getInstanceCatalogByCirc(self.centradeg, self.centdecdeg,
            radiusdeg, expmjd = self.expmjd, filter = self.filter)

  def getInstanceCatalogByCirc(self, centradeg, centdecdeg, radiusdeg, expmjd = 0.,
          filter = 'r'): 
    self.centradeg = centradeg
    self.centdecdeg = centdecdeg
    self.radiusdeg = radiusdeg
    self.filter = filter
    self.expmjd = expmjd

    objtype = self.objtype
    om = self.om[objtype]
    const = None
    self.ftype = om['formatas']
    if re.search("GALAXY", objtype) or re.search("AGN", objtype):
      '''We need to get galaxies from every tile in the overlap region
      '''
      colstr = ""
      const = None
      for omkey in om.keys():
        if omkey == "formatas":
          continue
        if omkey == "component":
  	  self.component = om[omkey]
	  continue
        self.ptype = om[omkey]['ptype']
        appint = int(self.em[omkey]['id'])
        map = om[omkey]
        cols = self.getUniqueColNames(map)
        colarr = []
        for k in cols.keys():
          if cols[k].startswith("'") and cols[k].endswith("'"):
            cols[k] = "'"+cols[k]+"'"
          colarr.append("%s as %s"%(cols[k], k))
        colarr.append("%i as appendint"%(appint))
        colstr = ",".join(colarr)
        if not len(map['constraint'].strip()) == 0:
          const = map['constraint']
          if const.startswith("%%"):
            const = const.lstrip("%%")
            const = eval(const)
          const = "WHERE %s"%(const)
      self.queries, self.coldesc = initGalaxy(self.centradeg, self.centdecdeg,
              self.radiusdeg, colstr, constraint = const)
      return self.getNextChunk()

    elif objtype == 'SSM':
      '''Need to do query and then do the ephemeris calculation
      '''
      if self.expmjd == 0.:
        raise Exception("Expmjd cannot be None if you want to get Solar System Objects out")
      colstr = ""
      const = None
      for omkey in om.keys():
        if omkey == "formatas":
          continue
        if omkey == "component":
  	  self.component = om[omkey]
	  continue
        self.ptype = om[omkey]['ptype']
        appint = int(self.em[omkey]['id'])
        map = om[omkey]
        cols = self.getUniqueColNames(map)
        colarr = []
        for k in cols.keys():
          colarr.append("%s as %s"%(cols[k], k))
        colarr.append("%i as appendint"%(appint))
        colstr = ",".join(colarr)
        if not len(map['constraint'].strip()) == 0:
          const = map['constraint']
          if const.startswith("%%"):
            const = const.lstrip("%%")
            const = eval(const)
          const = "WHERE %s"%(const)
      if (self.expmjd > 50093.1604):
          self.queries, self.coldesc = initSSM(self.centradeg, self.centdecdeg,
                  self.radiusdeg, self.expmjd, colstr, constraint = const)
      else:
          raise Exception("There are no indexed ephemerides for this observation time")
      
      return self.getNextChunk()

    else:
      '''Deal with all other types of objects
      '''
      for omkey in om.keys():
        if omkey == "formatas":
          continue
        if omkey == "component":
  	  self.component = om[omkey]
	  continue
        self.ptype = om[omkey]['ptype']
        break
      self.queries = self.getUnionQuery()
      return self.getNextChunk()

  def makeCatalogFromQuery(self, result, dithered=False):
    nic = InstanceCatalog(cd=self.catDescription, md=deepcopy(self.metadata))

    nic.objectType = self.ftype
    nic.neighborhoodType = self.component
    nic.catalogType = self.filetypes
    if self.opsimmeta is not None:
      for k in self.opsimmeta.keys():
        nic.metadata.addMetadata(k,eval("self.opsimmeta.%s"%(k)),"")
    else:
        nic.metadata.addMetadata('filter', self.filter, "filter of observation")
        nic.metadata.addMetadata('expmjd', self.expmjd, "mjd of observation")
        nic.metadata.addMetadata('centradeg', self.centradeg, "ra of center of field")
        nic.metadata.addMetadata('centdecdeg', self.centdecdeg, "dec of center of field")
    if dithered:
        omap = self.om['METADATA'][self.opsim]        
        mra =\
        eval("self.opsimmeta.%s"%(omap['moonrakey']))
        mdec =\
        eval("self.opsimmeta.%s"%(omap['moondeckey']))
        rottel =\
        eval("self.opsimmeta.%s"%(omap['rottelkey']))
        nalt, naz, nrotsky, ndist2moon =\
        nic.recalculatePointingInfo(math.radians(self.centradeg),\
                math.radians(self.centdecdeg), mra, mdec, self.expmjd, rottel)
        nic.metadata.addMetadata(omap['altkey'], nalt, "dithered altitude",
                clobber=True)
        nic.metadata.addMetadata(omap['azkey'], naz, "dithered azimuth",
                clobber=True)
        nic.metadata.addMetadata(omap['dist2moonkey'], ndist2moon, "dithered distance to the moon",
                clobber=True)
        nic.metadata.addMetadata(omap['rotskykey'], nrotsky, "Rotator angle after dithering",
                clobber=True)
    nic.catalogType = self.filetypes
    colkeys = zip(result[0].keys(),self.coldesc)
    for i,k in enumerate(colkeys):
      if k[0] == u'variabilityParameters':
        #need to cast to string in case the result is None which happens
        #when an object is not variable
        nic.addColumn( numpy.array([eval(str(s[k[0]])) for s in result]), k[1]['name'])
      else:
        nic.addColumn( numpy.array([s[k[0]] for s in result]), k[1]['name'])
    if nic.neighborhoodType == "EXTRAGALACTIC":
        nic.dataArray['raJ2000'] *= math.pi/180.
        nic.dataArray['decJ2000'] *= math.pi/180.
    if nic == None:
        raise RuntimeError, '*** nic is None'
    if nic.metadata == None:
        raise RuntimeError, '*** nic.metadata is None'
    if len(nic.dataArray) < 0:
        raise RuntimeError, '*** nic.dataArray has len < 0'

    return nic

  def deg2rad(self, ang):
    return ang*math.pi/180.

  def rad2deg(self, ang):
    return ang*180./math.pi

  def getRaDecBounds(self, bbox):
    decmin = bbox.getDecMin()
    decmax = bbox.getDecMax()
    ramin = bbox.getRaMin()
    ramax = bbox.getRaMax()

    bound = ""
    if ramin < 0 and ramax > 360:
      bound = "decl between %f and %f"%(decmin, decmax)
    elif ramax > 360:
      bound = "(ra between %f and 360. or ra between 0. and %f) and decl between %f and %f"%(ramin,ramax-360.,decmin,decmax)
    elif ramin < 0:
      bound = "(ra between %f and 360. or ra between 0. and %f) and decl between %f and %f"%(ramin+360.,ramax,decmin,decmax)
    else:
      bound = "ra between %f and %f and decl between %f and %f"%(ramin, ramax, decmin, decmax)
    return bound

  def getRaDecBoundsCirc(self, ra, dec, radius):
    ramax = ra+radius/math.cos(self.deg2rad(dec))
    ramin = ra-radius/math.cos(self.deg2rad(dec))
    decmax = dec+radius
    decmin = dec-radius
    bbox = Bbox(ramin,ramax,decmin,decmax)
    return self.getRaDecBounds(bbox)

  def getEnvironPath(self, var):
    if os.environ.has_key(var):
      return os.environ[var]
    else:
      raise Exception("Environment variable %s not set."%(var))

