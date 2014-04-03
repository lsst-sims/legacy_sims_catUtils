import scipy
import time
import lsst.sims.photUtils.Variability as variability
import numpy

"""
As of 3 April 2014, I no longer think this script runs.
It assumes a very old variability API.
"""

def getParams(filename):
    fh = open(filename)
    parr = []
    names = fh.readline().strip().split(",")
    string_flds = ("varMethodName", "filename", "lcfilename", "sedfilename",
            "lcfile", "spectrumname", "spectrum_name")
    for l in fh:
        params = {}
        flds = l.rstrip().split(",")
        for name, fld in zip(names,flds):
            if name in string_flds:
                params[name] = fld
            else:
                params[name] = float(fld)
        parr.append(params)
    return parr

if __name__ == "__main__": 
    startmjd = 51000.
    endmjd = 51100.
    steps = 4800
    var = variability.Variability(cache=True)
    mjds = numpy.linspace(startmjd, endmjd, steps)
    arr = getParams("mflare.dat")
    """
    for a in arr:
        fhout = open("lcs/mflare_%i.out"%(a['varsimobjid']),"w")
        t0 = time.time()
        dmags = eval("var.%s(a, mjds)"%(a['varMethodName']))
        line = []
        fhout.write("#MJD,u,g,r,i,z,y\n")
        for k in a:
            line.append("%s:%s"%(k,str(a[k])))
        fhout.write("#"+",".join(line)+"\n")
        print time.time() - t0
        for i in range(steps):
            line = [mjds[i], dmags['u'][i], dmags['g'][i], dmags['r'][i], dmags['i'][i],\
                    dmags['z'][i], dmags['y'][i]]
            fhout.write(",".join([str(el) for el in line])+"\n")
        fhout.close()
    """

    arr = getParams("bh_microlens.dat")
    for a in arr:
        fhout = open("lcs/bh_microlens_%i.out"%(a['varsimobjid']),"w")
        dmags = eval("var.%s(a, mjds)"%(a['varMethodName']))
        line = []
        fhout.write("#MJD,u,g,r,i,z,y\n")
        for k in a:
            line.append("%s:%s"%(k,str(a[k])))
        fhout.write("#"+",".join(line)+"\n")
        for i in range(steps):
            line = [mjds[i], dmags['u'][i], dmags['g'][i], dmags['r'][i], dmags['i'][i],\
                    dmags['z'][i], dmags['y'][i]]
            fhout.write(",".join([str(el) for el in line])+"\n")
        fhout.close()

    var = variability.Variability(cache=True)
    arr = getParams("amcvn.dat")
    for a in arr:
        fhout = open("lcs/amcvn_%i.out"%(a['varsimobjid']),"w")
        dmags = eval("var.%s(a, mjds)"%(a['varMethodName']))
        line = []
        fhout.write("#MJD,u,g,r,i,z,y\n")
        for k in a:
            line.append("%s:%s"%(k,str(a[k])))
        fhout.write("#"+",".join(line)+"\n")
        for i in range(steps):
            line = [mjds[i], dmags['u'][i], dmags['g'][i], dmags['r'][i], dmags['i'][i],\
                    dmags['z'][i], dmags['y'][i]]
            fhout.write(",".join([str(el) for el in line])+"\n")
        fhout.close()
    arr = getParams("cepheid.dat")
    for a in arr:
        fhout = open("lcs/cepheid_%i.out"%(a['varsimobjid']),"w")
        dmags = eval("var.%s(a, mjds)"%(a['varMethodName']))
        line = []
        fhout.write("#MJD,u,g,r,i,z,y\n")
        for k in a:
            line.append("%s:%s"%(k,str(a[k])))
        fhout.write("#"+",".join(line)+"\n")
        for i in range(steps):
            line = [mjds[i], dmags['u'][i], dmags['g'][i], dmags['r'][i], dmags['i'][i],\
                    dmags['z'][i], dmags['y'][i]]
            fhout.write(",".join([str(el) for el in line])+"\n")
        fhout.close()
    arr = getParams("microlens.dat")
    for a in arr:
        fhout = open("lcs/microlens_%i.out"%(a['varsimobjid']),"w")
        dmags = eval("var.%s(a, mjds)"%(a['varMethodName']))
        line = []
        fhout.write("#MJD,u,g,r,i,z,y\n")
        for k in a:
            line.append("%s:%s"%(k,str(a[k])))
        fhout.write("#"+",".join(line)+"\n")
        for i in range(steps):
            line = [mjds[i], dmags['u'][i], dmags['g'][i], dmags['r'][i], dmags['i'][i],\
                    dmags['z'][i], dmags['y'][i]]
            fhout.write(",".join([str(el) for el in line])+"\n")
        fhout.close()

    arr = getParams("eb.dat")
    for a in arr:
        fhout = open("lcs/eb_%i.out"%(a['varsimobjid']),"w")
        dmags = eval("var.%s(a, mjds)"%(a['varMethodName']))
        line = []
        fhout.write("#MJD,u,g,r,i,z,y\n")
        for k in a:
            line.append("%s:%s"%(k,str(a[k])))
        fhout.write("#"+",".join(line)+"\n")
        for i in range(steps):
            line = [mjds[i], dmags['u'][i], dmags['g'][i], dmags['r'][i], dmags['i'][i],\
                    dmags['z'][i], dmags['y'][i]]
            fhout.write(",".join([str(el) for el in line])+"\n")
        fhout.close()

    arr = getParams("agn.dat")
    mag_o = 20.
    for a in arr: 
        fhout = open("lcs/agn_%i.out"%(a['varsimobjid']),"w")
        dmags = eval("var.%s(a, mjds)"%(a['varMethodName']))
        line = []
        fhout.write("#MJD,u,g,r,i,z,y\n")
        for k in a:
            line.append("%s:%s"%(k,str(a[k])))
        fhout.write("#"+",".join(line)+"\n")
        for i in range(steps):
            line = [mjds[i], dmags['u'][i]+mag_o, dmags['g'][i]+mag_o,\
                    dmags['r'][i]+mag_o, dmags['i'][i]+mag_o,\
                    dmags['z'][i]+mag_o, dmags['y'][i]+mag_o]
            fhout.write(",".join([str(el) for el in line])+"\n")
        fhout.close()
        mag_o += 0.04

        
    arr = getParams("rrly.dat")
    for a in arr:
        fhout = open("lcs/rrly_%i.out"%(a['varsimobjid']),"w")
        dmags = eval("var.%s(a, mjds)"%(a['varMethodName']))
        line = []
        fhout.write("#MJD,u,g,r,i,z,y\n")
        for k in a:
            line.append("%s:%s"%(k,str(a[k])))
        fhout.write("#"+",".join(line)+"\n")
        for i in range(steps):
            line = [mjds[i], dmags['u'][i], dmags['g'][i], dmags['r'][i], dmags['i'][i],\
                    dmags['z'][i], dmags['y'][i]]
            fhout.write(",".join([str(el) for el in line])+"\n")
        fhout.close()

