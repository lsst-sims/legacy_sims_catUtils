import warnings
from sqlalchemy.orm import scoped_session, sessionmaker, mapper
from sqlalchemy.sql import expression
from sqlalchemy import create_engine
from sqlalchemy import ThreadLocalMetaData
import sqlalchemy.databases as sd 
from sqlalchemy import MetaData
from sqlalchemy import Table
from sqlalchemy.ext.sqlsoup import SqlSoup
from sqlalchemy import exc as sa_exc


class DBMSModel(object):
    """Database Interface"""
    def __init__(self):
        warnings.simplefilter("ignore", category=sa_exc.SAWarning)

        self.sessions = {}
        self.a_engine = create_engine("mssql+pymssql://LSST-2:L$$TUser@fatboy.npl.washington.edu:1433/LSST",
                                echo=False)
        self.a_session = scoped_session(sessionmaker(autoflush=True, 
                                                     bind=self.a_engine))
        self.a_metadata = MetaData()
        self.a_metadata.bind = self.a_engine

        self.opsim = Table('output_opsim3_61',
                           self.a_metadata, autoload=True)
        self.allstar = Table('starsALL_forceseek',
                             self.a_metadata, autoload=True)
        self.star = Table('starsMSRGB_forceseek',
                          self.a_metadata, autoload=True)
        self.starcirc = Table('starsMSRGBp80m10_forceseek',
                              self.a_metadata, autoload=True)
        self.wd = Table('starsWD_forceseek',
                        self.a_metadata, autoload=True)
        self.rrly = Table('starsRRLy_forceseek',
                          self.a_metadata, autoload=True)
        self.bhb = Table('starsBHB_forceseek',
                         self.a_metadata, autoload=True)
        self.lens = Table('lens',
                          self.a_metadata, autoload=True)
        self.image = Table('image',
                           self.a_metadata, autoload=True)
        self.eb = Table('ebstars',
                        self.a_metadata, autoload=True)
        self.cepheid = Table('cepheidstars',
                             self.a_metadata, autoload=True)
        self.astromeggs = Table('AstromEasterEggs',
                                self.a_metadata, autoload=True)
        self.dwarfgal = Table('dwarfGalaxies',
                              self.a_metadata, autoload=True)
        self.db = SqlSoup(self.a_metadata)
        self.Star = self.db.map(self.star,
                                primary_key=[self.star.c.simobjid])
        self.StarCirc = self.db.map(self.starcirc,
                                    primary_key=[self.starcirc.c.simobjid])
        self.AllStar = self.db.map(self.allstar,
                                   primary_key=[self.allstar.c.simobjid])
        self.AstromEggs = self.db.map(self.astromeggs)
        self.OpSim3_61 = self.db.map(self.opsim)
        self.Wd = self.db.map(self.wd, primary_key=[self.wd.c.simobjid])
        self.BHB = self.db.map(self.bhb, primary_key=[self.bhb.c.simobjid])
        self.RRLy = self.db.map(self.rrly, primary_key=[self.rrly.c.simobjid])
        self.LENS = self.db.map(self.lens, primary_key=[self.lens.c.id])
        self.IMAGE = self.db.map(self.image, primary_key=[self.image.c.id])
        self.EBSTARS = self.db.map(self.eb, primary_key=[self.eb.c.simobjid])
        self.CEPHEIDSTARS = self.db.map(self.cepheid,
                                        primary_key=[self.cepheid.c.simobjid])
        self.DWARFGALAXYSTARS = self.db.map(self.dwarfgal)

        self.sessions['MSSQL'] = self.a_session

    def initGalaxy(self, ra, dec, radiusdeg, columns, constraint=None):
        if constraint is not None:
            query = self.a_session.execute(
                ("EXECUTE [LSST].[dbo].[GalaxySearchSpecColsConstraint]\
                   @RaSearch = %f, @DecSearch = %f, @apertureRadius = %f,\
                   @ColumnNames = '%s', @WhereClause = '%s'"
                 % (ra, dec, radiusdeg * 60., columns, constraint)))
        else:
            query = self.a_session.execute(
                ("EXECUTE [LSST].[dbo].[GalaxySearchSpecColsConstraint]\
                   @RaSearch = %f, @DecSearch = %f, @apertureRadius = %f,\
                   @ColumnNames = '%s'"
                 % (ra, dec, radiusdeg*60., columns)))
        coldesc = []
        for k in query.keys():
            coldesc.append({"name":k})
        return query, coldesc

    def initSSM(self, ra, dec, radiusdeg, expmjd, columns, constraint=None):
        if constraint is not None:
            query = self.a_session.execute(
                ("select %s from [LSST].[dbo].[fSSMAll]("
                 "%f, %f, %f,%f) where %s" 
                 % (columns, expmjd, ra, dec, radiusdeg*60., constraint)))
        else:
            query = self.a_session.execute(
                ("select %s from [LSST].[dbo].[fSSMAll]("
                 "%f, %f, %f,%f)"
                 % (columns, expmjd, ra, dec, radiusdeg*60.)))
        coldesc = []
        for k in query.keys():
            coldesc.append({"name":k})
        return query, coldesc
