"""
This script will validate the MBA tables on fatboy by querying a random
sample of asteroids and using the orbit propagation and ephemeride
generation code in sims_movingObjects to verify that the position parametrized
on fatboy is within 5 milliarcseconds of the position calculated directly by
pyoorb.

Running this script will require installing and setting up sims_movingObjects.

You will also need to compile the coordinate transformation code coordtrans.c
by running
    gcc -lm -o coord coordtrans.c
"""

from __future__ import with_statement, print_function
from subprocess import check_output

def equatorial_from_ecliptic(lon, lat):
    """
    Convert ecliptic coordinates to RA, Dec.
    Not vectorized; must do coordinate pairs one-at-a-time.

    Parameters:
    -----------
    lon is ecliptic longitude in degrees

    lat is ecliptic latitude in degrees

    Returns:
    -------
    RA in degrees

    Dec in degrees
    """

    try:
        result = check_output(["./coord", "-e", "%.6f" % lon, "%.6f" % lat])
    except OSError:
        print('\nTrying to call ./coord resulted in an exception\n'
              'Are you sure you have compiled coordtrans.c?\n'
              '    gcc -lm -o coord coordtrans.c\n')
        raise
    result = result.split('\n')
    ra_dec = result[1]
    ra_dec = ra_dec.split(':')[1]
    ra_dec = ra_dec.strip().split()
    return 15.0*float(ra_dec[0]), float(ra_dec[1])

import numpy as np
import os
import cStringIO
from lsst.sims.movingObjects import Orbits, PyOrbEphemerides
from lsst.sims.catUtils.baseCatalogModels import MBAObj
from lsst.sims.catalogs.db import DBObject
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import _angularSeparation
from lsst.sims.utils import arcsecFromRadians

import time

def validate_orbits(obs, db, data_files=None, data_dir=None,
                    chunk_size=10000):
    """
    Take a telescope pointing, find all of the asteroids within that
    pointing, and validate the orbits as parametrized on fatboy against
    the orbits as directly evaluated by pyoorb

    Parameters
    ----------
    obs is an ObservationMetaData characterizing the telescope pointing

    db is a CatalogDBObject connecting to the Solar System object table
    we are currently validating

    data_files is a list of the data files containing the original orbit
    parameters

    data_dir is the directory containing the data files with the original
    orbit paramters

    chunk_size is an int indicating how many asteroids to process at once
    (this will balance performance versus memory usage)

    Returns
    -------
    The maximum difference between the position parametrized on fatboy
    and the position calculated by pyoorb in milliarcseconds.

    The number of objects tested.
    """

    if not hasattr(validate_orbits, 'ephemerides'):
        # construct a PyOrbEphemrides instantiation
        validate_orbits.ephemerides = PyOrbEphemerides()

    if not hasattr(validate_orbits, 'max_d'):
        validate_orbits.max_d = -1.0
        validate_orbits.n_obj = 0

    if not hasattr(validate_orbits, 'name_lookup'):
        # construct a lookup table associating objid with
        # the names of the objects as stored in the .des files
        dtype = np.dtype([('name', str, 8), ('id', int)])
        try:
            name_db = DBObject(database='LSSTSSM', host='fatboy.phys.washington.edu',
                               port=1433, driver='mssql+pymssql')
        except:
            name_db = DBObject(database='LSTSSM', host='localhost',
                               port=51433, driver='mssql+pymssql')

        query = 'SELECT name, ssmid from C14_mba_name_map'
        chunk_iter = name_db.get_arbitrary_chunk_iterator(query=query, dtype=dtype,
                                                          chunk_size=10000)

        validate_orbits.name_lookup = {}
        for chunk in chunk_iter:
            for line in chunk:
                validate_orbits.name_lookup[line['id']] = line['name']

        print('built name look up dict')

    if not hasattr(validate_orbits, 'data_dir') or validate_orbits.data_dir != data_dir:
        # create a StringIO object containing all of the data from
        # the original .des files
        validate_orbits.data_dir = data_dir
        validate_orbits.des_cache = {}
        validate_orbits.header = None
        first_file = True
        for file_name in data_files:
            if not file_name.endswith('.s3m'):
                continue
            if data_dir is not None:
                full_name = os.path.join(data_dir, file_name)
            else:
                full_name = file_name

            with open(full_name, 'r') as input_file:
                for ix, line in enumerate(input_file):
                    if ix==0:
                        if first_file:
                            validate_orbits.header = line
                            print('read header from %s' % file_name)
                            first_file = False
                        continue
                    split_line = line.strip().split()
                    validate_orbits.des_cache[split_line[0]] = line

        print('finished initialization %d' % len(validate_orbits.des_cache))

    colnames = ['objid', 'raJ2000', 'decJ2000']

    results = db.query_columns(colnames=colnames, obs_metadata=obs,
                               chunk_size=chunk_size)

    t_start = time.time()
    for chunk in results:
        validate_orbits.n_obj += len(chunk)
        orbit_obj = Orbits()
        orbit_buffer = cStringIO.StringIO()
        orbit_buffer.write(validate_orbits.header)
        for asteroid in chunk:
            ast_name = validate_orbits.name_lookup[asteroid['objid']]
            orbit_buffer.write(validate_orbits.des_cache[ast_name])
        orbit_buffer.seek(0)
        orbit_obj.readOrbits(orbit_buffer)
        validate_orbits.ephemerides.setOrbits(orbit_obj)
        eph_list = validate_orbits.ephemerides.generateEphemerides(np.array([obs.mjd.UTC]), byObject=False)
        ra_vec = np.radians(np.array([rr for rr in eph_list['ra'][0]]))
        dec_vec = np.radians(np.array([dd for dd in eph_list['dec'][0]]))
        dd = 1000.0*arcsecFromRadians(_angularSeparation(chunk['raJ2000'], chunk['decJ2000'],
                                                         ra_vec, dec_vec))
        orbit_buffer.close()
        ellapsed = time.time()-t_start
        if dd.max() > validate_orbits.max_d:
            max_dex = np.argmax(dd)
            print('object %d displacement %e TAI %.12f' % (chunk['objid'][max_dex], dd.max(), obs.mjd.TAI))
            validate_orbits.max_d = dd.max()

    return validate_orbits.max_d, validate_orbits.n_obj

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--data', type=str, default=None, nargs='+',
                        help='List of files containing the original orbit '
                             'parameters to be used in validating the model')

    parser.add_argument('--data_dir', type=str, default=None,
                        help='Directory containing the files from DATA')

    parser.add_argument('--n_fields', type=int, default=1,
                        help='Number of random fields of view to validate')

    parser.add_argument('--seed', type=int, default=99,
                        help='seed for random number generator')

    parser.add_argument('--chunk_size', type=int, default=10000,
                        help='number of asteroids to process at once; '
                             'affect performance and memory usage')

    args = parser.parse_args()
    if args.data is None:
        if args.data_dir is None:
            raise RuntimeError("Need to specify at least DATA_DIR, if not DATA")
        data_files = os.listdir(args.data_dir)
    elif isinstance(args.data, list):
        data_files = [args.data]
    else:
        data_files = args.data

    try:
        # if you are on UW campus/VPN
        mba_db = MBAObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                        port=1433, driver='mssql+pymssql')
    except:
        # if you need to use the SSH tunnel to connect to fatboy
        mba_db = MBAObj()

    rng = np.random.RandomState(args.seed)
    lon_list = rng.random_sample(args.n_fields)*360.0
    lat_list = rng.random_sample(args.n_fields)*10.0-5.0
    mjd_list = rng.random_sample(args.n_fields)*3653.0+59580.0

    max_d = -1.0
    n_obj = 0
    for lon, lat, mjd in zip(lon_list,lat_list, mjd_list):
        ra, dec = equatorial_from_ecliptic(lon, lat)

        obs = ObservationMetaData(mjd=mjd,
                                  pointingRA=ra,
                                  pointingDec=dec,
                                  boundLength=1.75,
                                  boundType='circle')

        local_max_d, local_n_obj = validate_orbits(obs, mba_db,
                                                   data_files=data_files,
                                                   data_dir=args.data_dir,
                                                   chunk_size=args.chunk_size)

        if local_max_d>max_d:
            max_d = local_max_d
        n_obj += local_n_obj

    print('max_d %e; n_obj %d ' % (max_d, n_obj))
