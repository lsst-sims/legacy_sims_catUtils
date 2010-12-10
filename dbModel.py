from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.sql import expression
from sqlalchemy import create_engine
from sqlalchemy import ThreadLocalMetaData
import sqlalchemy.databases as sd 
from sqlalchemy import func
from sqlalchemy import schema
from elixir import *


a_engine = create_engine("postgresql://cosmouser:cosmouser@172.25.79.34/cosmoDB.11.19.2009",
        echo=False, server_side_cursors=True)
a_session = scoped_session(sessionmaker(autoflush=True, 
    bind=a_engine))
a_metadata = metadata
a_metadata.bind = a_engine

b_engine = create_engine("postgresql://jobreporter:jobreporter@172.25.79.34/joblog",
        echo=False)
b_session = application_session = scoped_session(sessionmaker(autoflush=True,
     bind=b_engine))
b_metadata = ThreadLocalMetaData()
b_metadata.bind = b_engine

c_engine = create_engine("postgresql://calibuser:calibuser@128.95.99.32/calibDB.05.05.2010",
        echo=False)
c_session = application_session = scoped_session(sessionmaker(autoflush=True,
     bind=c_engine))
c_metadata = ThreadLocalMetaData()
c_metadata.bind = c_engine

class Star(Entity):
  using_options(tablename="stars", autoload=True, metadata=a_metadata,
          session=a_session)

class CalibStar(Entity):
  using_options(tablename="msrgb", autoload=True, metadata=c_metadata,
          session=c_session)
  using_mapper_options(primary_key=['simobjid'])

class Wd(Entity):
  using_options(tablename="starswd", autoload=True, metadata=a_metadata,
          session=a_session)

class Galaxy(Entity):
  using_options(tablename="galaxy_11182010", autoload=True, metadata=a_metadata,
          session=a_session)

class GalaxyRect(Entity):
  using_options(tablename="galaxy", autoload=True, metadata=a_metadata,
          session=a_session)
  using_mapper_options(primary_key=['id'])

class Ephems(Entity):
  using_options(tablename="ephems", autoload=True, metadata=a_metadata,
          session=a_session)
  using_mapper_options(primary_key=['objid','obsdate'])

class Orbits(Entity):
  using_options(tablename="orbits", autoload=True, metadata=a_metadata,
          session=a_session)

class Tiles(Entity):
  using_options(tablename="tiles", autoload=True, metadata=a_metadata,
          session=a_session)

class OpSim3_61(Entity):
  using_options(tablename="output_opsim3_61", autoload=True,
          metadata=a_metadata, session=a_session)

class CatalogEventLog (Entity):
  using_options(tablename='eventlog', metadata=b_metadata, session=b_session)
  jobid = Field(Integer, index=True)
  owner = Field(UnicodeText)
  pkey = Field(UnicodeText)
  pvalue = Field(UnicodeText)
  time = Field(DateTime(timezone=True))
  taskNumber = Field(Integer)
  ip = Field(sd.postgres.PGInet)
  description = Field(UnicodeText)
  def __repr__(self):
    return '<Log Event (%s,%s) at %s>' % (self.pkey, self.pvalue, self.time)

class JobStateLog (Entity):
  using_options(tablename='statelog', metadata=b_metadata, session=b_session)
  jobid = Field(Integer, index=True)
  owner = Field(UnicodeText)
  pkey = Field(UnicodeText)
  pvalue = Field(UnicodeText)
  time = Field(DateTime(timezone=True))
  def __repr__(self):
    return '<Log state (%s,%s) at %s>' % (self.pkey, self.pvalue, self.time)

