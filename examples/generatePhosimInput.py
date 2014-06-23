from __future__ import with_statement
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import DBObject, ObservationMetaData
from lsst.sims.catalogs.measures.example_utils.exampleCatalogDefinitions.trimExamples import TrimCatalogPoint, TrimCatalogSersic2D, TrimCatalogZPoint


starObjNames = ['msstars', 'bhbstars', 'wdstars', 'rrlystars', 'cepheidstars']

obsMD = DBObject.from_objid('opsim3_61')
obs_metadata = obsMD.getObservationMetaData(88625744, 0.33, makeCircBounds = True)

doHeader= True
for starName in starObjNames:
    print "starName"
    stars = DBObject.from_objid(starName)
    sfd_trim=TrimCatalogPoint(stars,obs_metadata=obs_metadata)
    if (doHeader):
        with open("sfd_output.txt","w") as fh:
            sfd_trim.write_header(fh)
        doHeader = False
    sfd_trim.write_catalog("sfd_output.txt",write_mode='a',write_header=False,chunk_size=20000)


gals = DBObject.from_objid('galaxyBulge')
sfd_trim = TrimCatalogSersic2D(gals, obs_metadata=obs_metadata,constraint="sedname_bulge is not NULL")
sfd_trim.write_catalog("sfd_output.txt",write_mode='a',write_header=False,chunk_size=20000)

gals = DBObject.from_objid('galaxyDisk')
sfd_trim = TrimCatalogSersic2D(gals, obs_metadata=obs_metadata,constraint="sedname_disk is not NULL")
sfd_trim.write_catalog("sfd_output.txt",write_mode='a',write_header=False,chunk_size=20000)

gals = DBObject.from_objid('galaxyAgn')
sfd_trim = TrimCatalogZPoint(gals, obs_metadata=obs_metadata,constraint="sedname_agn is not NULL")
sfd_trim.write_catalog("sfd_output.txt",write_mode='a',write_header=False,chunk_size=20000)



