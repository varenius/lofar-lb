#!/usr/bin/python
import pyrap.tables as pt
import sys

if not len(sys.argv) == 2:
    print "Usage: simple2circ.py <ms name>"
    sys.exit()

msname = sys.argv[1]
t=pt.taql ("update %s/ set DATA = mscal.stokes(DATA,'circ')" % msname)
t=pt.taql ("update %s/POLARIZATION set CORR_TYPE=[5,6,7,8]" % msname)
t.close()
