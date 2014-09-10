import math
from lsst.sims.catalogs.generation.db import DBObject, ObservationMetaData
#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catUtils.exampleCatalogDefinitions import RefCatalogGalaxyBase, PhoSimCatalogPoint,\
                                                         PhoSimCatalogZPoint, PhoSimCatalogSersic2D
from lsst.sims.catUtils import makeObsParamsAzAltTel, makeObsParamsRaDecTel

def exampleReferenceCatalog():
    obs_metadata = ObservationMetaData(circ_bounds=dict(ra=0., dec=0., radius=0.01))
    dbobj = DBObject.from_objid('galaxyBase')
    t = dbobj.getCatalog('ref_catalog_galaxy', obs_metadata=obs_metadata)
    filename = 'test_reference.dat'
    t.write_catalog(filename, chunk_size=10)

def examplePhoSimCatalogs():
    obsMD = DBObject.from_objid('opsim3_61')
    obs_metadata = obsMD.getObservationMetaData(88544919, 0.1, makeCircBounds=True)
    objectDict = {}
    objectDict['testStars'] = {'dbobj':DBObject.from_objid('msstars'),
                               'constraint':None,
                               'filetype':'phoSim_catalog_POINT',
                               'obsMetadata':obs_metadata}
    objectDict['testGalaxyBulge'] = {'dbobj':DBObject.from_objid('galaxyBulge'),
                               'constraint':"mass_bulge > 1. and sedname_bulge is not NULL",
                               'filetype':'phoSim_catalog_SERSIC2D',
                               'obsMetadata':obs_metadata}
    objectDict['testGalaxyDisk'] = {'dbobj':DBObject.from_objid('galaxyDisk'),
                               'constraint':"DiskLSSTg < 20. and sedname_disk is not NULL",
                               'filetype':'phoSim_catalog_SERSIC2D',
                               'obsMetadata':obs_metadata}
    objectDict['testGalaxyAgn'] = {'dbobj':DBObject.from_objid('galaxyAgn'),
                               'constraint':"sedname_agn is not NULL",
                               'filetype':'phoSim_catalog_ZPOINT',
                               'obsMetadata':obs_metadata}

    for objKey in objectDict.keys():
        dbobj = objectDict[objKey]['dbobj']
        t = dbobj.getCatalog(objectDict[objKey]['filetype'],
                             obs_metadata=objectDict[objKey]['obsMetadata'], 
                             constraint=objectDict[objKey]['constraint'])

        print
        print "These are the required columns from the database:"
        print t.db_required_columns()
        print
        print "These are the columns that will be output to the file:"
        print t.column_outputs
        print
    
        filename = 'catalog_test_%s.dat'%(dbobj.objid)
        print "querying and writing catalog to %s:" % filename
        t.write_catalog(filename)
        filename = 'catalog_test_%s_chunked.dat'%(dbobj.objid)
        t.write_catalog(filename, chunk_size=10)
        print " - finished"

def examplePhoSimNoOpSim():
    raDeg= 15.
    decDeg = -30.
    mjd = 51999.75
    md =  makeObsParamsRaDecTel(math.radians(raDeg), math.radians(decDeg), mjd, 'r')
    obs_metadata_rd = ObservationMetaData(circ_bounds=dict(ra=raDeg,
                                                        dec=decDeg,
                                                        radius=0.1),
                                                        mjd=mjd,
                                                        bandpassName='r',
                                                        metadata=md)
    azRad = math.radians(220.)
    altRad = math.radians(79.)
    md = makeObsParamsAzAltTel(azRad, altRad, mjd, 'r')
    raDeg = math.degrees(md['Unrefracted_RA'][0])
    decDeg = math.degrees(md['Unrefracted_Dec'][0])
    obs_metadata_aa = ObservationMetaData(circ_bounds=dict(ra=raDeg,
                                                        dec=decDeg,
                                                        radius=0.1),
                                                        mjd=mjd,
                                                        bandpassName='r',
                                                        metadata=md)
    dbobj = DBObject.from_objid('msstars')
    t = dbobj.getCatalog('phoSim_catalog_POINT', obs_metadata= obs_metadata_rd)
    t.write_catalog('catalog_test_stars_rd.dat')
    t = dbobj.getCatalog('phoSim_catalog_POINT', obs_metadata= obs_metadata_aa)
    t.write_catalog('catalog_test_stars_aa.dat')

if __name__ == '__main__':
    exampleReferenceCatalog()
    examplePhoSimCatalogs()
    examplePhoSimNoOpSim()
