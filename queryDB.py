#!/usr/bin/env python
from dbMsModel import *
from copy import deepcopy
import re
import os
import math
import numpy
import pyoorb
from lsst.sims.catalogs.generation.config import ConfigObj
import lsst.sims.catalogs.generation.movingObjects as mo
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.astrometry import Bbox
from sqlalchemy import select 
from sqlalchemy.orm import join
from sqlalchemy.orm import aliased
from sqlalchemy.sql.expression import text
from sqlalchemy.sql.expression import select
from sqlalchemy.sql.expression import between


class queryDB(object):
  def __init__(self, objtype = 'STARS', filetypes=('TRIM',), chunksize=100000):
    if os.environ.has_key("CATALOG_DESCRIPTION_PATH"):
      catalogDescriptionPath = os.environ["CATALOG_DESCRIPTION_PATH"]
    else:
      raise Exception("Environment variable CATALOG_DESCRIPTION_PATH not set to location of the catalog description files")
    dbMapConfigFile = catalogDescriptionPath+"requiredFields.dat"
    objConfigFile = catalogDescriptionPath+"objectMap.dat"
    objEnumConfigFile = catalogDescriptionPath+"objectEnum.dat"
    metaConfigFile = catalogDescriptionPath+"requiredMetadata.dat"
    self.nictemp = InstanceCatalog(catalogDescriptionPath+"/config.dat")
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
    self.component = None
    self.expmjd = None
    self.centradeg = None
    self.centdecdeg = None
    self.radiusdeg = None
    self.filter = None
    self.opsimmeta = None
    self.opsim = ""
    self.curtile = None

  def getNextChunk(self):
    result = []
    '''
    for q in self.queries:
      if len(result) == self.chunksize:
          break
      else:
          result += q.fetchmany(self.chunksize - len(result))  
    '''
    result = self.queries.fetchmany(self.chunksize)
    if len(result) == 0:
      return None
    else:
      cat = self.makeCatalogFromQuery(result)
      return cat
    """
    except Exception, e:
      print "Exception of type: %s"%(type(e))
      raise Exception(e)
    """

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
    eval("self.opsimmeta.%s"%(omap['radegkey']))
    self.centdecdeg =\
    eval("self.opsimmeta.%s"%(omap['decdegkey']))
    self.expmjd =\
    eval("self.opsimmeta.%s"%(omap['expmjdkey']))
    self.filter =\
    eval("self.opsimmeta.%s"%(omap['filterkey']))

  def addSpacial(self, query, map, type, ta):
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
      query = self.addSpacial(query, map, "circle", ta)
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
    query = queries[0]
    self.coldesc = query.column_descriptions
    for i in range(len(queries)-1):
      query = query.union_all(queries[i+1])

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
      dmap = self.dm[ft][self.component][self.ptype]
      for k in dmap.keys():
        colstr = dmap[k][1]
        if colstr.startswith("%%"):
          colstr = colstr.lstrip("%%")
          colstr = eval(colstr)
        if k == map['idkey']:
          continue
        else:
          if cols.has_key(k):
            print "Replacing expression %s with %s for key %s"%(cols[k], colstr, k)
            cols[k] = colstr
          else:
            cols[k] = colstr 
    return cols

  def getUniqueColNamesMO(self, idkey):
    cols = {}
    for ft in self.filetypes:
      dmap = self.dm[ft]['SOLARSYSTEM']['MOVINGPOINT']
      for k in dmap.keys():
        colstr = dmap[k][1]
        if k == idkey:
          continue
        else:
          if cols.has_key(k):
            print "Replacing expression %s with %s for key %s"%(cols[k], colstr, k)
            cols[k] = colstr
          else:
            cols[k] = colstr 
    return cols

  def getInstanceCatalogByBbox(self, bbox, expmjd=0., filter='r'): 
    self.expmjd = expmjd
    self.centradeg = bbox.getCentDeg()[0]
    self.centdecdeg = bbox.getCentDeg()[1]
    self.filter = filter
    objtype = self.objtype
    om = self.om[objtype]
    if om.has_key('formatas'):
      self.ptype = om['formatas']
    else:
      self.ptype = om[om.keys()[0]]['ptype']
    if re.search("GALAXY", objtype):
      '''We need to get galaxies from every tile in the overlap region
      '''
      tiles = self.getTilesBbox(bbox)
      queries = []
      for tile in tiles:
	self.curtile = tile
        queries += self.getUnionQuery("point @ spoly '%s' and (point + strans(0,\
			              %f*PI()/180.,  %f*PI()/180., 'XYZ')) @ spoly\
				      '%s'"%(tile['tbox'],-tile['decmid'],tile['ramid'],tile['bbox']))
        #queries += self.getQueryList("point @ spoly '%s' and point @ spoly\
    	#		              '%s'"%(tile['tbox'],tile['tbbox']))
      self.queries = queries
      return self.getNextChunk()
    elif objtype == 'SSM':
      '''Need to do query and then do the ephemeris calculation
      '''
      if self.expmjd == 0.:
        raise Exception("Expmjd cannot be None if you want to get Solar System Objects out")
      oom = self.om['ORBIT']
      #We need the orbit map for the join.  I'm assuming that if the orbit map
      #is length 1 there is only one orbit table and if > 1 there is one orbit
      #table per ephemeride table.
      queries = []
      n = 0
      for k in om.keys():
        if k == "formatas":
          continue
        map = om[k]
        #equery = session.query(eval("%s.objid"%(map['table'])))
        if len(oom) == 1:
          omap = oom['ORBIT0']
        elif len(oom) == len(om):
          omap = oom['ORBIT%i'%(n)]
        else:
          raise Exception('Getting orbits...', 'The orbit map and ephemeride\
                  map do not agree in length')
        eid = eval("%s.%s.label(\"%s\")"%(map['table'],
                      self.dm[self.filetypes[0]][self.component][self.ptype][map['idkey']][1],
                      map['idkey']))
        oid = eval("%s.%s.label(\"%s\")"%(omap['table'],
                      self.dm[self.filetypes[0]][self.component][self.ptype][omap['idkey']][1],
                      omap['idkey']))
        equery = session.query(eid,oid).filter(eid == oid)
        #equery = self.addUniqueCols(map, equery)
        equery = self.addUniqueCols(omap, equery)
        if not len(map['constraint'].strip()) == 0:
          const = map['constraint'].strip()
          if const.startswith("%%"):
            const = const.lstrip("%%")
            const = eval(const)
          equery = equery.filter(const)
        if not len(omap['constraint'].strip()) == 0:
          const = omap['constraint'].strip()
          if const.startswith("%%"):
            const = const.lstrip("%%")
            const = eval(const)
          equery = equery.filter(const)
        query = equery.filter(self.getRaDecBounds(bbox))
        queries.append(session.execute(query))
        n += 1
      self.queries = queries
      return self.getNextChunk()
    else:
      filterstr = "point @ sbox \'((%fd,%fd),(%fd,%fd))\'"%(
		  bbox.getRaMin(), bbox.getDecMin(), bbox.getRaMax(), bbox.getDecMax())
      self.queries = self.getUnionQuery(filter=filterstr)
      return self.getNextChunk()

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
    if om.has_key('formatas'):
      self.ptype = om['formatas']
    else:
      self.ptype = om[om.keys()[0]]['ptype']
    if re.search("GALAXY", objtype):
      '''We need to get galaxies from every tile in the overlap region
      '''
      #THis is a hack
      self.component = "EXTRAGALACTIC"
      self.queries, self.coldesc = initGalaxy(self.centradeg, self.centdecdeg, self.radiusdeg, om.keys()[2])
      return self.getNextChunk()
    elif objtype == 'SSM':
      '''Need to do query and then do the ephemeris calculation
      '''
      if self.expmjd == 0.:
        raise Exception("Expmjd cannot be None if you want to get Solar System Objects out")
      oom = self.om['ORBIT']
      #We need the orbit map for the join.  I'm assuming that if the orbit map
      #is length 1 there is only one orbit table and if > 1 there is one orbit
      #table per ephemeride table.
      queries = []
      n = 0
      for k in om.keys():
        if k == "formatas":
          continue
        map = om[k]
        #equery = session.query(eval("%s.objid"%(map['table'])))
        if len(oom) == 1:
          omap = oom['ORBIT0']
        elif len(oom) == len(om):
          omap = oom['ORBIT%i'%(n)]
        else:
          raise Exception('Getting orbits...', 'The orbit map and ephemeride\
                  map do not agree in length')
        eid = eval("%s.%s.label(\"%s\")"%(map['table'],
                      self.dm[self.filetypes[0]][self.component][self.ptype][map['idkey']][1],
                      map['idkey']))
        oid = eval("%s.%s.label(\"%s\")"%(omap['table'],
                      self.dm[self.filetypes[0]][self.component][self.ptype][omap['idkey']][1],
                      omap['idkey']))
        equery = session.query(eid,oid).filter(eid == oid)
        #equery = self.addUniqueCols(map, equery)
        equery = self.addUniqueCols(omap, equery)
        if not len(map['constraint'].strip()) == 0:
          const = map['constraint'].strip()
          if const.startswith("%%"):
            const = const.lstrip("%%")
            const = eval(const)
          equery = equery.filter(const)
        if not len(omap['constraint'].strip()) == 0:
          const = omap['constraint'].strip()
          if const.startswith("%%"):
            const = const.lstrip("%%")
            const = eval(const)
          equery = equery.filter(const)
        query = equery.filter(self.getRaDecBoundsCirc(self.centradeg,
            self.centdecdeg, radiusdeg))
        queries.append(session.execute(query))
        n += 1
      self.queries = queries
      return self.getNextChunk()
    else:
      self.queries = self.getUnionQuery()
      return self.getNextChunk()

  def makeMovingObjectsFromOrbitList(self, results):
    objects = []
    ephem_datfile = ""
    pyoorb.pyoorb.oorb_init(ephemeris_fname=ephem_datfile)
    for r in results:
      mymo = mo.MovingObject(r['q'], r['e'], r['i'], r['node'],
                             r['argPeri'], r['timePeri'], r['epoch'],
                             magHv=r['magHv'], phaseGv=r['phaseGv'], index=r['index'],
                             n_par=r['n_par'], moid=r['moid'], 
                             objid=r['id'], objtype=r['objtype'],
                             isVar=r['isVar'], var_t0=r['var_t0'],
                             var_timescale=r['var_timescale'],
                             var_fluxmax=r['var_fluxmax'],
                             sedname=r['sedname'],
                             u_opp=r['u_opp'],g_opp=r['g_opp'], r_opp=r['r_opp'],
                             i_opp=r['i_opp'], z_opp=r['z_opp'], y_opp=r['y_opp'])      
      objects.append(mymo)
    # turn list of moving objects into movingObjectList object
    objects = mo.MovingObjectList(objects)
    # generate ephemerides for all objects at once
    objects.generateEphemeridesForAllObjects([self.expmjd], obscode=807)
    self.thismjd = mymo.mjdTaiStr(self.expmjd)
    # return the basic list of moving objects 
    objects = objects._mObjects
    return objects



  class dbMovingObject(mo.MovingObject):
    def keys(self):
      return self.keyarr
    def setkeys(keyarr):
      self.keyarr = keyarr

  def makeCatalogFromQuery(self, result):
    nic = deepcopy(self.nictemp)
    nic.objectType = self.ptype
    nic.neighborhoodType = self.component
    nic.catalogType = self.filetypes
    if self.opsimmeta is not None:
      for k in self.opsimmeta.keys():
        nic.metadata.addMetadata(k,eval("self.opsimmeta.%s"%(k)),"")
    nic.catalogType = self.filetypes
    data = {}
    if self.ptype == "MOVINGPOINT":
      result = self.makeMovingObjectsFromOrbitList(result)
      thismjd = self.thismjd
      om = self.om[self.objtype]
      colkeys = self.getUniqueColNamesMO(om[om.keys()[1]]['idkey'])
      for k in colkeys:
        data[k] = []
      for s in result:
        for k in colkeys:
          col = colkeys[k]
          if colkeys[k].startswith("%%"):
            col = col.lstrip("%%")
            col = eval(col)
            eval("data[k].append(%s)"%(col))
          else:
            eval("data[k].append(s.%s)"%(col))
      for k in colkeys:
        arr = numpy.asarray(data[k])
        nic.addColumn(arr, k)
    else:
      colkeys = zip(result[0].keys(),self.coldesc)
      for k in colkeys:
        nic.addColumn( numpy.array([s[k[0]] for s in result]), k[1]['name'])


    if nic == None:
        raise RuntimeError, '*** nic is None'
    if nic.metadata == None:
        raise RuntimeError, '*** nic.metadata is None'
    if len(nic.dataArray) < 0:
        raise RuntimeError, '*** nic.dataArray has len < 1'

    return nic

  def getTilesBbox(self, bbox):
    sq = session.query(Tiles.id, Tiles.ramid, Tiles.decmid, Tiles.bbox).filter("sbox\
            '((%fd, %fd), (%fd,%fd))' && bbox"%(bbox.getRaMin(),bbox.getDecMin(),
            bbox.getRaMax(),bbox.getDecMax())).subquery().alias(name="mytiles")
    col1 = expression.literal_column("spoly '{(%fd, %fd), (%fd, %fd), (%fd, %fd), (%fd, %fd)}'+\
            strans(-%s.ramid*PI()/180., %s.decmid*PI()/180.,0.0,\
            'ZYX')"%(\
            bbox.getRaMin(),bbox.getDecMin(), bbox.getRaMin(), bbox.getDecMax(),
	    bbox.getRaMax(),bbox.getDecMax(), bbox.getRaMax(),
            bbox.getDecMin(),"mytiles","mytiles")).label("tbox")
    col2 = expression.literal_column("bbox +\
            strans(-%s.ramid*PI()/180., %s.decmid*PI()/180.,0.0,\
            'ZYX')"%("mytiles","mytiles")).label("tbbox")
    query = session.query(col1,col2,sq.c.id,sq.c.ramid,sq.c.decmid,sq.c.bbox)
    result = query.all()
    tiles = []
    for r in result:
      tiles.append({'tid':r.id, 'tbox':r.tbox, 'tbbox':r.tbbox, 'ramid':r.ramid, 'decmid':r.decmid,
            'bbox':r.bbox})
    return tiles
  def getTilesCirc(self, raDeg, decDeg, radiusDeg):
    sq = session.query(Tiles.id, Tiles.ramid, Tiles.decmid, Tiles.bbox).filter("scircle\
            '<(%fd, %fd), %fd>' && bbox"%(raDeg, decDeg,
            radiusDeg)).subquery().alias(name="mytiles")
    col = expression.literal_column("scircle '<(%fd, %fd), %fd>'+\
            strans(-%s.ramid*PI()/180., %s.decmid*PI()/180.,0.0,\
            'ZYX')"%(raDeg,
            decDeg, radiusDeg,"mytiles","mytiles")).label("circ")
    query = session.query(col,sq.c.id,sq.c.ramid,sq.c.decmid,sq.c.bbox)
    result = query.all()
    tiles = []
    for r in result:
      tiles.append({'tid':r.id, 'circ':r.circ, 'ramid':r.ramid, 'decmid':r.decmid,
            'bbox':r.bbox})
    return tiles

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
