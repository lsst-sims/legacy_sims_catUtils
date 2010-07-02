#!/usr/bin/env python
from dbModel import *
import os
import math
import numpy
import pyoorb
from lsst.sims.catalogs.generation.config import ConfigObj
import lsst.sims.catalogs.generation.movingObjects as mo
from lsst.sims.catalogs.measures.instance import InstanceCatalog


class queryDB(object):
  def __init__(self, objtype = 'STARS', filetypes=('TRIM',), chunksize=100000):
    if os.environ.has_key("CATALOG_DESCRIPTION_PATH"):
      catalogDescriptionPath = os.environ["CATALOG_DESCRIPTION_PATH"]
    else:
      raise Exception("Environment variable CATALOG_DESCRIPTION_PATH not set to location of the catalog description files")
    dbMapConfigFile = catalogDescriptionPath+"requiredFields.dat"
    objConfigFile = catalogDescriptionPath+"objectMap.dat"
    metaConfigFile = catalogDescriptionPath+"requiredMetadata.dat"
    setup_all()
    self.filetypes = filetypes
    self.objtype = objtype
    self.chunksize=chunksize
    self.dm = ConfigObj(dbMapConfigFile)
    self.om = ConfigObj(objConfigFile)
    self.mm = ConfigObj(metaConfigFile)
    self.queries = None
    self.ptype = 'POINT'
    self.opsimmeta = None
    self.opsim = ""

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
      return self.makeCatalogFromQuery(result, self.ptype)
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

  def getQueryList(self, filter="", joinmap=None):  
    queries = []
    for omkey in self.om[self.objtype].keys():
      if omkey == "formatas":
        continue
      map = self.om[self.objtype][omkey]
      query = session.query(eval("%s.%s.label(\"%s\")"%(map['table'],
          self.dm[self.filetypes[0]][map['component']][map['ptype']][map['idkey']][1],
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
      dmap = self.dm[ft][map['component']][map['ptype']]
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

  def getInstanceCatalogById(self, id, opsim="OPSIM361", radiusdeg=2.1, add_columns=()): 
    self.getOpsimMetadataById(id, opsim)
    objtype = self.objtype
    om = self.om[objtype]
    if om.has_key('formatas'):
      self.ptype = om['formatas']
    else:
      self.ptype = om[om.keys()[0]]['ptype']
    if objtype == 'GALAXY':
      '''We need to get galaxies from every tile in the overlap region
      '''
      tiles = self.getTiles(self.opsimmeta.Unrefracted_RA,
              self.opsimmeta.Unrefracted_Dec, radiusdeg)
      queries = []
      for tile in tiles:
        queries += self.getQueryList("point @ scircle '%s' and (point + strans(0,\
			              %f*PI()/180.,  %f*PI()/180., 'XYZ')) @ spoly\
				      '%s'"%(tile['circ'],-tile['decmid'],tile['ramid'],tile['bbox']))
      self.queries = queries
      return self.getNextChunk()
    elif objtype == 'SSM':
      '''Need to do query and then do the ephemeris calculation
      '''
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
                      self.dm[self.filetypes[0]][map['component']][map['ptype']][map['idkey']][1],
                      map['idkey']))
        oid = eval("%s.%s.label(\"%s\")"%(omap['table'],
                      self.dm[self.filetypes[0]][omap['component']][omap['ptype']][omap['idkey']][1],
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
        query = equery.filter(self.getRaDecBounds(self.opsimmeta.Unrefracted_RA,
            self.opsimmeta.Unrefracted_Dec, radiusdeg))
        queries.append(session.execute(query))
        n += 1
      self.queries = queries
      return self.getNextChunk()
    else:
      filterstr = "point @ scircle \'<(%fd,%fd),%fd>\'"%(self.opsimmeta.Unrefracted_RA, self.opsimmeta.Unrefracted_Dec, radiusdeg)
      self.queries = self.getQueryList(filter=filterstr)
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
      mymo.calcEphemeris([self.opsimmeta.Opsim_expmjd])
      objects.append(mymo)
    self.thismjd = objects[0].mjdTaiStr(self.opsimmeta.Opsim_expmjd)
    return objects

  class dbMovingObject(mo.MovingObject):
    def keys(self):
      return self.keyarr
    def setkeys(keyarr):
      self.keyarr = keyarr

  def makeCatalogFromQuery(self, result, ptype):
    if os.environ.has_key("CATALOG_DESCRIPTION_PATH"):
      catalogDescriptionPath = os.environ["CATALOG_DESCRIPTION_PATH"]
    else:
      raise Exception("Environment variable CATALOG_DESCRIPTION_PATH not set to location of the catalog description files")
    nic = InstanceCatalog(catalogDescriptionPath+"/config.dat")
    nic.objectType = ptype
    for k in self.opsimmeta.keys():
      nic.metadata.addMetadata(k,eval("self.opsimmeta.%s"%(k)),"")
    nic.catalogType = self.filetypes
    data = {}
    if ptype == "MOVINGPOINT":
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

  def deg2rad(self, ang):
    return ang*math.pi/180.

  def rad2deg(self, ang):
    return ang*180./math.pi

  def getRaDecBounds(self, ra, dec, radius):
    ramax = ra+radius/math.cos(self.deg2rad(dec))
    ramin = ra-radius/math.cos(self.deg2rad(dec))
    decmax = dec+radius
    decmin = dec-radius
    if decmin < -90:
      decmin = -90
      if((-90 - decmin) > (decmax + 90)):
        decmax = -190 - decmin
      else:
        decmax = decmax
    elif decmax > 90:
      decmax = 90
      if((decmax - 90) > (90 - decmin)):
        decmin = 180 - decmax
      else:
        decmin = decmin
    else:
      pass
        

    bound = ""
    if ramin < 0 and ramax > 360:
      bound = "decl between %f and %f"%(decmin, decmax)
    elif ramax > 360:
      bound = "(ra between %f and 0. or ra between 0. and %f) and decl between %f and %f"%(ramin,ramax-360.,decmin,decmax)
    elif ramin < 0:
      bound = "(ra between %f and 360. or ra between 0. and %f) and decl between %f and %f"%(ramin+360.,ramax,decmin,decmax)
    else:
      bound = "ra between %f and %f and decl between %f and %f"%(ramin, ramax, decmin, decmax)
    return bound
