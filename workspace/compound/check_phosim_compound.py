import lsst.sims.catUtils.baseCatalogModels as db_models
from lsst.sims.catalogs.definitions import CompoundInstanceCatalog
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.exampleCatalogDefinitions import (PhoSimCatalogSersic2D,
PhoSimCatalogPoint, PhoSimCatalogZPoint)

class connectionMixin(object):
        database='LSSTCATSIM'
        host='fatboy.phys.washington.edu'
        driver='mssql+pymssql'
        port=1433
    
class disk_class(connectionMixin,db_models.GalaxyDiskObj):
    pass

class bulge_class(connectionMixin,db_models.GalaxyBulgeObj):
    pass
    
class agn_class(connectionMixin,db_models.GalaxyAgnObj):
    pass
    
class star_class(connectionMixin,db_models.StarObj):
    pass
    

def testCompoundCatalog():
    """
    This test writes a PhoSim input catalog and compares it, one line at a time
    to a previously written catalog that should be identical.

    is test uses CompoundInstanceCatalog
    """
    obs_metadata = ObservationMetaData(pointingRA=22.0, pointingDec=-29.0,
                                       rotSkyPos=32.1, mjd=69680.0,
                                       boundType='circle',
                                       boundLength=0.1)

    bulgeDB = bulge_class()
           
    diskDB = disk_class()
           

    agnDB = agn_class()

    starDB = star_class()
           

    # first, generate the catalog without a CompoundInstanceCatalog
    single_catName = 'old_fashioned_phosim_cat.txt'

    testBulge = PhoSimCatalogSersic2D(bulgeDB, obs_metadata = obs_metadata)
    testDisk = PhoSimCatalogSersic2D(diskDB, obs_metadata = obs_metadata)
    testAgn = PhoSimCatalogZPoint(agnDB, obs_metadata = obs_metadata)
    testStar = PhoSimCatalogPoint(starDB, obs_metadata = obs_metadata)

    testBulge.write_catalog(single_catName, write_header=False,
                            chunk_size=10000)
    testDisk.write_catalog(single_catName, write_header=False, write_mode='a',
                           chunk_size=10000)
    testAgn.write_catalog(single_catName, write_header=False, write_mode='a',
                          chunk_size=10000)
    testStar.write_catalog(single_catName, write_header=False, write_mode='a',
                           chunk_size=10000)

    # now, generate the catalog using CompoundInstanceCatalog
    #
    # because the CompoundCatalogDBObject requires that database
    # connection parameters be set in the input CatalogDBObject
    # daughter class definitions, we have to declare dummy
    # CatalogDBObject daughter classes below

    
    compoundCatalog = CompoundInstanceCatalog([PhoSimCatalogSersic2D, PhoSimCatalogSersic2D,
                                               PhoSimCatalogZPoint, PhoSimCatalogPoint],
                                              [disk_class, bulge_class,
                                               agn_class, star_class],
                                              obs_metadata=obs_metadata,
                                              compoundDBclass=db_models.GalaxyTileCompoundObj)


    compound_catName = 'compound_phosim_cat.txt'

    compoundCatalog.write_catalog(compound_catName, write_header=False,
                                  chunk_size=10000)

    # verify that the two catalogs are equivalent
    with open(single_catName, 'r') as single_file:
        with open(compound_catName, 'r') as compound_file:
            single_lines = single_file.readlines()
            compound_lines = compound_file.readlines()
            assert len(single_lines) == len(compound_lines)

            for line in single_lines:
                if line not in compound_lines:
                    raise RuntimeError('\n%s\nnot in compound' % line)

if __name__ == "__main__":
     testCompoundCatalog()
