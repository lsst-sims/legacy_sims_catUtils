"""configConnection.py:
This example demonstrates the different ways to set the host/port/database
for the backend.

Note that the connection to the database is made during the baseModel instantiation.
Therefore if the configuration is changed after the object is created,
there is no change in the backend connection.
"""
import os
from lsst.sims.catUtils.baseCatalogModels import StarObj, BaseCatalogConfig
from lsst.utils import getPackageDir


def demoEditDefaultFile():
    """Demonstrate changing database connection config through default config file

    This method cannot be used to make changes at run time

    Before importing baseCatalogModels, edit the default config file at:
    $SIMS_CATUTILS_DIR/config/db.py
    """
    #after editing $SIMS_CATUTILS_DIR/CONFIG/db.py
    #defaults will be applied at class definition during import
    from lsst.sims.catUtils.baseCatalogModels import SolarSystemObj
    ssmDB = SolarSystemObj()


def demoGlobalConfig():
    """Demonstrate changing database connection in script.

    This method can be used to change default connections at runtime,
    but only before instantiating a DbObject
    """
    CONFIG = BaseCatalogConfig(driver='mssql+pymssql',
                               port='5555',
                               host='localhost',
                               database='LSST')
    #Can pass to the constructor
    starDB = StarObj(**CONFIG.toDict())


def demoLoadFromFile():
    """Demonstrate loading database connection config from config file.

    This method can be used to change default connection at runtime,
    but only before instantiating a DbObject.
    """
    config = BaseCatalogConfig()
    configFilename = os.path.join(getPackageDir("sims_catUtils"), "config", "db.py")
    config.load(configFilename)
    starDB = StarObj(**config.toDict())


if __name__ == '__main__':
    demoEditDefaultFile()
    demoGlobalConfig()
    demoLoadFromFile()
