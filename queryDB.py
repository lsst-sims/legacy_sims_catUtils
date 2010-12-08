#!/usr/bin/env python
from dbModel import *
import os
import math
import numpy
import pyoorb
from lsst.sims.catalogs.generation.config import ConfigObj
import lsst.sims.catalogs.generation.movingObjects as mo
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.astrometry import Bbox


class queryDB(object):
  def __init__(self, objtype = 'STARS', filetypes=('TRIM',), chunksize=100000):
    if os.environ.has_key("CATALOG_DESCRIPTION_PATH"):
      catalogDescriptionPath = os.environ["CATALOG_DESCRIPTION_PATH"]
    else:
      raise Exception("Environment variable CATALOG_DESCRIPTION_PATH not set to location of the catalog description files")
    if os.environ.has_key("SED_DATA"):
      sedDataPath = os.environ["SED_DATA"]
    else:
      raise Exception("Environment variable CATALOG_DESCRIPTION_PATH not set to location of the catalog description files")
    dbMapConfigFile = catalogDescriptionPath+"requiredFields.dat"
    objConfigFile = catalogDescriptionPath+"objectMap.dat"
    metaConfigFile = catalogDescriptionPath+"requiredMetadata.dat"
    setup_all()
    self.filetypes = filetypes
    self.objtype = objtype
    self.chunksize=chunksize
    self.rootSEDdir=sedDataPath
    self.dm = ConfigObj(dbMapConfigFile)
    self.om = ConfigObj(objConfigFile)
    self.mm = ConfigObj(metaConfigFile)
    self.queries = None
    self.ptype = 'POINT'
    self.component = None
    self.expmjd = None
    self.centradeg = None
    self.centdecdeg = None
    self.filter = None
    self.opsimmeta = None
    self.opsim = ""
    self.curtile = None

  def getNextChunk(self):
    result = []
    for q in self.queries:
      if len(result) == self.chunksize:
          break
      else:
          result += q.fetchmany(self.chunksize - len(result))  
    if len(result) == 0:
      return None
    else:
      return self.makeCatalogFromQuery(result)
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

  def getQueryList(self, filter=""):  
    queries = []
    for omkey in self.om[self.objtype].keys():
      if omkey == "formatas":
        continue
      if omkey == "component":
        self.component = self.om[self.objtype][omkey]
        continue
      map = self.om[self.objtype][omkey]
      idcolstr = self.dm[self.filetypes[0]][self.component][map['ptype']][map['idkey']][1]
      if idcolstr.startswith("%%"):
          idcolstr = idcolstr.lstrip("%%")
          idcolstr = eval(idcolstr)
      query = session.query(eval("%s.%s.label(\"%s\")"%(map['table'],
          idcolstr,
          map['idkey'])))
      query = self.addUniqueCols(map, query)
      query = query.filter(filter)
      if not len(map['constraint'].strip()) == 0:
        const = map['constraint']
        if const.startswith("%%"):
          const = const.lstrip("%%")
          const = eval(const)
        query = query.filter(const)
      queries.append(session.execute(query))
    return queries

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
    if objtype == 'GALAXY' or objtype == 'ASSEMBLEDGALAXY':
      '''We need to get galaxies from every tile in the overlap region
      '''
      tiles = self.getTilesBbox(bbox)
      queries = []
      for tile in tiles:
	self.curtile = tile
        queries += self.getQueryList("point @ spoly '%s' and (point + strans(0,\
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
        if k == "component":
          self.component = om['component']
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
                      self.dm[self.filetypes[0]][self.component][map['ptype']][map['idkey']][1],
                      map['idkey']))
        oid = eval("%s.%s.label(\"%s\")"%(omap['table'],
                      self.dm[self.filetypes[0]][self.component][omap['ptype']][omap['idkey']][1],
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
      self.queries = self.getQueryList(filter=filterstr)
      return self.getNextChunk()

  def getInstanceCatalogById(self, id, opsim="OPSIM361", radiusdeg=2.1):
    self.getOpsimMetadataById(id, opsim)
    return self.getInstanceCatalogByCirc(self.centradeg, self.centdecdeg,
            radiusdeg, expmjd = self.expmjd, filter = self.filter)

  def getInstanceCatalogByCirc(self, centradeg, centdecdeg, radiusdeg, expmjd = 0.,
          filter = 'r'): 
    self.centradeg = centradeg
    self.centdecdeg = centdecdeg
    self.filter = filter
    self.expmjd = expmjd

    objtype = self.objtype
    om = self.om[objtype]
    if om.has_key('formatas'):
      self.ptype = om['formatas']
    else:
      self.ptype = om[om.keys()[0]]['ptype']
    if objtype == 'GALAXY' or objtype == 'ASSEMBLEDGALAXY':
      '''We need to get galaxies from every tile in the overlap region
      '''
      tiles = self.getTilesCirc(self.centradeg,
              self.centdecdeg, radiusdeg)
      queries = []
      for tile in tiles:
	self.curtile = tile
        queries += self.getQueryList("point @ scircle '%s' and (point + strans(0,\
			              %f*PI()/180.,  %f*PI()/180., 'XYZ')) @ spoly\
				      '%s'"%(tile['circ'],-tile['decmid'],tile['ramid'],tile['bbox']))
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
        if k == "component":
          self.component = om[k]
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
                      self.dm[self.filetypes[0]][self.component][map['ptype']][map['idkey']][1],
                      map['idkey']))
        oid = eval("%s.%s.label(\"%s\")"%(omap['table'],
                      self.dm[self.filetypes[0]][self.component][omap['ptype']][omap['idkey']][1],
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
      bbox = self.getRaDecBoundsCirc(self.centradeg, self.centdecdeg,
              radiusdeg)
      filterstr = "%s and point @ scircle \'<(%fd,%fd),%fd>\'"%(bbox,self.centradeg,
              self.centdecdeg, radiusdeg)
      self.queries = self.getQueryList(filter=filterstr)
      return self.getNextChunk()

  def makeMovingObjectsFromOrbitList(self, results):
    objects = []
    ephem_datfile = ""
    pyoorb.pyoorb.oorb_init(ephemeris_fname=ephem_datfile)
    for r in results:
      #Hack since there are no sedfilenames in the db at the moment
      mymo = mo.MovingObject(r['q'], r['e'], r['i'], r['node'],
                             r['argPeri'], r['timePeri'], r['epoch'],
                             magHv=r['magHv'], phaseGv=r['phaseGv'], index=r['index'],
                             n_par=r['n_par'], moid=r['moid'], 
                             objid=r['id'], objtype=r['objtype'],
                             isVar=r['isVar'], var_t0=r['var_t0'],
                             var_timescale=r['var_timescale'],
                             var_fluxmax=r['var_fluxmax'],
                             sedname="S.dat",
                             u_opp=r['u_opp'],g_opp=r['g_opp'], r_opp=r['r_opp'],
                             i_opp=r['i_opp'], z_opp=r['z_opp'], y_opp=r['y_opp'])  
      '''
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
      '''
      objects.append(mymo)
    # turn list of moving objects into movingObjectList object
    objects = mo.MovingObjectList(objects)
    # generate ephemerides for all objects at once
    objects.generateEphemeridesForAllObjects([self.expmjd], obscode=807)
    objects.calcAllMags(self.filter, [self.expmjd],
            self.rootSEDdir+"/ssmSed/",withErrors=False)
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
    if os.environ.has_key("CATALOG_DESCRIPTION_PATH"):
      catalogDescriptionPath = os.environ["CATALOG_DESCRIPTION_PATH"]
    else:
      raise Exception("Environment variable CATALOG_DESCRIPTION_PATH not set to location of the catalog description files")
    nic = InstanceCatalog(catalogDescriptionPath+"/config.dat")
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
      colkeys = self.getUniqueColNamesMO(om[om.keys()[2]]['idkey'])
      idkey = om[om.keys()[2]]['idkey']
      for k in colkeys:
        data[k] = []
      for s in result:
        for k in colkeys:
          if k == idkey:
            continue
          col = colkeys[k]
          if colkeys[k].startswith("%%"):
            col = col.lstrip("%%")
            col = eval(col)
            data[k].append(col)
          else:
            eval("data[k].append(s.%s)"%(col))
      for k in colkeys:
        arr = numpy.asarray(data[k])
        nic.addColumn(arr, k)
    else:
      colkeys = result[0].keys()
      for k in colkeys:
        data[k] = []
      for s in result:
        for k in colkeys:
          eval("data[k].append(s.%s)"%(k))
      for k in colkeys:
        arr = numpy.asarray(data[k])
        nic.addColumn(arr, k)


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
