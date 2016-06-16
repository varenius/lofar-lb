#!/usr/bin/python
import pyrap.tables as pt
import os,sys

if not len(sys.argv) == 2:
    print "Usage: %s <ms name>" % sys.argv[0]
    sys.exit()

msname = sys.argv[1]
if not os.path.exists(msname):
    print msname + " doesn't exist!"
    sys.exit()

antennatable = pt.table(msname + "::ANTENNA")
print "There are " + str(len(antennatable)) + " antennas"
