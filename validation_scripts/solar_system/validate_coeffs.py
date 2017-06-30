"""
This script will loop through the tables of MBA orbit coefficients, verifying
that every object in our model has a set of coefficients assigned for every
time period simulated.
"""

from __future__ import with_statement, print_function
from lsst.sims.catalogs.db import DBObject
import numpy as np
import time

def verify_month(month):
    """
    month is the starting MJD of the month to be verified (an int)

    returns a message summarizing the results
    """

    msg = ""

    if not hasattr(verify_month, 'db'):
        try:
            # if you are on the University of Washington campus or VPN
            verify_month.db = DBObject(database='LSSTSSM',
                                       host='fatboy.phys.washington.edu',
                                       port=1433, driver='mssql+pymssql')
        except:
            # try connecting through an SSH tunnel
            verify_month.db = DBObject(database='LSSTSSM',
                                       host='localhost',
                                       port=51433,
                                       driver='mssql_pymssql')

    if not hasattr(verify_month, 'n_id'):
        # first, get a list of all of the object IDs in the MBA model
        dtype = np.dtype([('id', int)])
        results = verify_month.db.execute_arbitrary('SELECT ssmid from C14_MBA_name_map',
                                                    dtype=dtype)
        id_list = results['id']
        verify_month.n_id = len(id_list)


    dtype = np.dtype([('id', int), ('mjd_start', float), ('mjd_end', float)])

    t_start = time.time()
    query = 'SELECT ssmid, mjd_start, mjd_end FROM C14_MBA_%d '% month
    query += 'ORDER BY mjd_start'

    prev_end = np.ones(verify_month.n_id+1, dtype=float)*month

    tol = 1.0e-9

    results = verify_month.db.get_arbitrary_chunk_iterator(query, dtype=dtype,
                                                           chunk_size=10000)
    row_ct = 0
    overhead = time.time()-t_start
    t_start = time.time()
    for chunk in results:
        n_chunk = len(chunk)
        sub_divisions = 0
        is_okay = False
        while not is_okay:
            sub_divisions += 1
            is_okay = True
            for i_sub in range(sub_divisions):
                i_start = n_chunk*i_sub/sub_divisions
                i_end = n_chunk*(i_sub+1)/sub_divisions
                unq = np.unique(chunk['id'][i_start:i_end])
                if len(unq) != i_end-i_start:
                    is_okay = False
                    break

        for i_sub in range(sub_divisions):
            i_start = n_chunk*i_sub/sub_divisions
            i_end = n_chunk*(i_sub+1)/sub_divisions
            id_vec = chunk['id'][i_start:i_end]
            mjd_start_vec = chunk['mjd_start'][i_start:i_end]
            mjd_end_vec = chunk['mjd_end'][i_start:i_end]

            invalid = np.where(np.abs(mjd_start_vec-prev_end[id_vec])>tol)
            try:
                assert len(invalid[0]) == 0
            except AssertionError:
                msg += "\nFound a gap in coverage\n"
                for id_val, start_val, end_val, prev_val in \
                zip(id_vec[invalid], mjd_start_vec[invalid], mjd_end_vec[invalid],
                    prev_end[id_vec][invalid]):

                    msg += "id %d\nmjd_start %.12f\nmjd_end %.12f\nprev_end %.12f\n\n" % \
                    (id_val, start_val, end_val, prev_val)

            prev_end[id_vec] = mjd_end_vec

        row_ct += n_chunk
        elapsed = time.time()-t_start
        expected = overhead + 2.09e8*elapsed/row_ct
        expected = expected/(24.0*60.0*60.0)
        #print('%d %.2f %.2f %d' % (row_ct,elapsed+overhead, expected, sub_divisions))

    invalid = np.where(np.abs(prev_end-(month+30.0))>tol)
    try:
        assert len(invalid[0]) == 1
    except AssertionError:
        msg += "\nFailed on final check of prev_end\n"
        msg += "    month %d\n" % month
        for id_val, end_val in zip(invalid[0], prev_end[invalid]):
            msg+="id %d\nprev_end %.12f\n\n" % (id_val, end_val)
        return msg

    msg += '\nMonth %d analyzed %d rows\n' % (month, row_ct)
    return msg


import argparse
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--month', type=int, default=None,
                        nargs='+',
                        help='List of MJDs denoting the database tables '
                             'to verify (5980 + n*30)')
    parser.add_argument('--out_file', type=str, default='mba_validation_log.txt',
                        help='file to write results to')

    args = parser.parse_args()
    if args.month is None:
        exit(0)
    elif not isinstance(args.month, list):
        month_list = [args.month]
    else:
        month_list = args.month

    if os.path.exists(args.out_file):
        raise RuntimeError('%s exists' % args.out_file)

    for month in month_list:
        t_start = time.time()
        msg = verify_month(month)
        elapsed = time.time()-t_start
        with open(args.out_file, 'a') as out_file:
            out_file.write('month %d completed in %.2f hours: %s\n' %
                           (month, elapsed/3600.0, msg))
