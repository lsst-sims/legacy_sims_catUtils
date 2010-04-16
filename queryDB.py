#!/usr/bin/env python
from dbModel import *
from catalogDbMap import catalogDbMap
import os
import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import CatalogDescription
from lsst.sims.catalogs.measures.instance import Metadata
from lsst.sims.catalogs.measures.instance import CatalogType


class queryDB(object):
  def __init__(self, objtype = 'star', filetypes=CatalogType.TRIM, chunksize=100000):
    setup_all()
    self._start=0
    self.filetypes = filetypes
    self.objtype = objtype
    self.chunksize=chunksize
    self.cdm = catalogDbMap()

  def getNextChunk(self):
    try:
      result = self.query.slice(self._start, self._start+self.chunksize).all()
      self._start += self.chunksize
      if len(result) == 0:
        return None
      else:
        return self.makeCatalogFromQuery(result)
    except Exception, e:
      print "Exception of type: %s"%(type(e))
      raise Exception(e)


  def getInstanceCatalogById(self, id,  opsim="3_61", add_columns=()):
    os = OpSim3_61.query.filter("obshistid=%s"%(id)).first()
    objtype = self.objtype
    self.metadata = Metadata()
    for k in self.cdm.metadataMap.keys():
      self.metadata.addMetadata(k,os.__dict__[self.cdm.metadataMap[k]['opsim3_61']],"")
    if objtype == 'star':
      self.query = session.query(Star.id)
      for k in self.cdm.objectTypes['POINT'].keys():
	if k == 'id':
	  continue
        elif k == 'magNorm':
          col = expression.literal_column(self.cdm.objectTypes['POINT'][k]['star']%(os.expmjd, os.filter)).label(k)
	  self.query = self.query.add_column(col)
	else:
          col = expression.literal_column(self.cdm.objectTypes['POINT'][k]['star']).label(k)
          self.query = self.query.add_column(col)
      self.query = self.query.filter("point @ scircle \'<(%fd,%fd),%fd>\'"%(os.fieldradeg, os.fielddecdeg, 2.1))
      result = self.query.slice(self._start, self._start+self.chunksize).all()

      self._start += self.chunksize
      if len(result) == 0:
        return None
      else:
        return self.makeCatalogFromQuery(result)
    elif objtype == 'wd':
      avs = schema.Column('ebv').label('galacticEbv')
      mags = func.toMag(fscale).label('magNorm')
      self.query = Wd.query.add_column(mags).add_column(avs).filter("point @ scircle \'<(%fd,%fd),%fd>\'"%(os.fieldradeg, os.fielddecdeg, 2.1))
    elif objtype == 'ssm':
      pass
    elif objtype == 'galaxy':
      pass
    else:
      raise Exception('getInstanceCatalogById', 'Did not give valid object type')

  def makeCatalogFromQuery(self, result):
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
    for k in self.cdm.objectTypes['POINT'].keys():
      data[k] = []
    for s in result:
      for k in self.cdm.objectTypes['POINT'].keys():
        exec("data[k].append(s.%s)"%(k))
    for k in self.cdm.objectTypes['POINT'].keys():
      arr = numpy.asarray(data[k])
      nic.addColumn(arr, k)
    nic.metadata = self.metadata
    return nic
