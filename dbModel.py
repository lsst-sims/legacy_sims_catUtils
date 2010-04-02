from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.sql import expression
from sqlalchemy import create_engine
from sqlalchemy import ThreadLocalMetaData
import sqlalchemy.databases as sd 
from sqlalchemy import func
from sqlalchemy import schema
from elixir import *


a_engine = create_engine("postgresql://cosmouser:cosmouser@172.25.79.34/cosmoDB.11.19.2009?server_side_cursors",
        echo=False)
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

class Star(Entity):
  using_options(tablename="stars", autoload=True, metadata=a_metadata,
          session=a_session)

class OpSim3_61(Entity):
  using_options(tablename="output_opsim3_61", autoload=True,
          metadata=a_metadata, session=a_session)

class CatalogEventLog (Entity):
  using_options(tablename='eventlog', metadata=b_metadata, session=b_session)
  jobid = Field(Integer, index=True)
  pkey = Field(Unicode(15))
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
  pkey = Field(Unicode(15))
  pvalue = Field(UnicodeText)
  time = Field(DateTime(timezone=True))
  def __repr__(self):
    return '<Log state (%s,%s) at %s>' % (self.pkey, self.pvalue, self.time)

