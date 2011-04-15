#!/usr/bin/env python
from dbMsModel import *
import pickle,time
from copy import deepcopy
import re
import os
import math
import numpy
import warnings
import pyoorb
from lsst.sims.catalogs.generation.config import ConfigObj
import lsst.sims.catalogs.generation.movingObjects as mo
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
          makeCat=True, pickleRes=False):
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
    self.pickleRes = pickleRes

  def closeSession(self):
    session.close_all()

  def getNextChunk(self):
    result = self.queries.fetchmany(self.chunksize)
    if len(result) == 0:
      self.closeSession()
      return None
    else:
      if self.pickleRes:
        self.pickleResults(result)
      if self.makeCat:
        cat = self.makeCatalogFromQuery(result)
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
      dmap = self.dm[ft][self.component][map['ptype']]
      for k in dmap.keys():
        colstr = dmap[k][1]
        if colstr.startswith("%%") and (str(map['ptype']) != 'MOVINGPOINT'):
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
    if re.search("GALAXY", objtype) or re.search("AGN", objtype):
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
    const = None
    if om.has_key('formatas'):
      self.ptype = om['formatas']
    else:
      self.ptype = om[om.keys()[0]]['ptype']
    if re.search("GALAXY", objtype) or re.search("AGN", objtype):
      '''We need to get galaxies from every tile in the overlap region
      '''
      colstr = ""
      const = ""
      for omkey in self.om[self.objtype].keys():
        if omkey == "formatas":
          continue
        if omkey == "component":
  	  self.component = self.om[self.objtype][omkey]
	  continue
        appint = int(self.em[omkey]['id'])
        map = self.om[self.objtype][omkey]
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
      oom = self.om['ORBIT']['ORBITS']
      self.component = self.om[self.objtype]["component"]
      appint = int(self.em['EPHEMS']['id'])
      cols = self.getUniqueColNames(oom)
      colarr = []
      colarr.append("%i as appendint"%(appint))
      for k in cols.keys():
        if cols[k].startswith("'") and cols[k].endswith("'"):
          cols[k] = "'"+cols[k]+"'"
        colarr.append("%s as %s"%(cols[k], k))
      colstr = ",".join(colarr)
      if not len(oom['constraint'].strip()) == 0:
        const = oom['constraint']
        if const.startswith("%%"):
          const = const.lstrip("%%")
          const = eval(const)
        const = "WHERE %s"%(const)
      self.queries, self.coldesc = initSSM(self.centradeg, self.centdecdeg,
              self.radiusdeg, self.expmjd, colstr, constraint = const)
      return self.getNextChunk()
    else:
      self.queries = self.getUnionQuery()
      return self.getNextChunk()

  def makeMovingObjectsFromOrbitList(self, results):
    objects = []
    ephem_datfile = ""
    pyoorb.pyoorb.oorb_init(ephemeris_fname=ephem_datfile)
    appendint = None
    for r in results:
      appendint = r['appendint']
      mymo = mo.MovingObject(r['q'], r['e'], r['i'], r['node'],
                             r['argPeri'], r['timePeri'], r['epoch'],
                             magHv=r['magHv'], phaseGv=r['phaseGv'],
                             index=r['idx'],
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
    objects = objects.getMovingObjectsInFieldofView(self.centradeg,\
            self.centdecdeg, self.radiusdeg, self.expmjd)
    objects.calcAllMags("imsim",[self.expmjd],os.path.join(self.getEnvironPath("SED_DATA"),"ssmSED"),withErrors=False)
    self.thismjd = mymo.mjdTaiStr(self.expmjd)
    # return the basic list of moving objects 
    objects = objects._mObjects
    return objects, appendint


  def makeCatalogFromQuery(self, result):
    nic = InstanceCatalog(cd=self.catDescription, md=deepcopy(self.metadata))
    nic.objectType = self.ptype
    nic.neighborhoodType = self.component
    nic.catalogType = self.filetypes
    if self.opsimmeta is not None:
      for k in self.opsimmeta.keys():
        nic.metadata.addMetadata(k,eval("self.opsimmeta.%s"%(k)),"")
    nic.catalogType = self.filetypes
    data = {}
    if self.ptype == "MOVINGPOINT":
      result,appendint = self.makeMovingObjectsFromOrbitList(result)
      thismjd = self.thismjd
      om = self.om[self.objtype]
      colkeys = self.getUniqueColNames(om['EPHEMS'])
      for k in colkeys:
        data[k] = []
      data['appendint'] = []
      for s in result:
        for k in colkeys:
          col = colkeys[k]
          if colkeys[k].startswith("%%"):
            col = col.lstrip("%%")
            col = eval(col)
            eval("data[k].append(%s)"%(col))
          else:
            eval("data[k].append(s.%s)"%(col))
        data['appendint'].append(appendint)
      for k in colkeys:
        arr = numpy.asarray(data[k])
        nic.addColumn(arr, k)
      nic.addColumn(numpy.asarray(data['appendint']), 'appendint')
    else:
      colkeys = zip(result[0].keys(),self.coldesc)
      #arr = numpy.asarray(result)
      #arr = arr.T
      for i,k in enumerate(colkeys):
        if k[0] == u'variabilityParameters':
          #need to cast to string in case the result is None which happens
          #when an object is not variable
          nic.addColumn( numpy.array([eval(str(s[k[0]])) for s in result]), k[1]['name'])
          #nic.addColumn(arr[i], k[1]['name'])
        else:
          nic.addColumn( numpy.array([s[k[0]] for s in result]), k[1]['name'])
          #nic.addColumn(arr[i], k[1]['name'])
      # nic.addColumn(numpy.fromiter((tuple(s[k[0]] for k in colkeys) for s
      #    in result), count=len(result)), k[1]['name'])
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

  def pickleResults(self, results):
    rfh = open("/astro/net/pogo3/krughoff/results.pkl","w")
    icfh = open("/astro/net/pogo3/krughoff/cols.pkl","w")
    pickle.dump(results,rfh,protocol=2)
    pickle.dump(self.coldesc,icfh,protocol=2)
