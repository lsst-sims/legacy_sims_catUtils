import sqlite3


def array2dbrecords(data):
    """
    takes an iterable such as an array, structured array or list and converts it into a a list of tuples suitable
    for insertion into a sqlite 3 database. 
    args:
        data: mandatory,
            list or input iterable
    returns:
        suitable list of tuples
        
    """
    l = list()
    for row in data:
        l.append(tuple(row))
    return l


def insertfromdata(tablename, records, multiple=True):
    """
    construct string to insert multiple records into sqlite3 database 
    args:
        tablename:
        records:
        multiple:
    returns:
    """
    if multiple:
        lst = records[0]
    else:
        lst = records
    s = 'INSERT INTO ' + str(tablename) + ' VALUES '
    s += "( " + ", ".join(["?"]*len(lst)) + ")"
    return s
