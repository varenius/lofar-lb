# Script to make parmdb-solutions time independent. Modified from original by Stephen Bourke
# Run e.g. using script.py *.parmdb, and it will output one .parmdb.timeext file for each parmdb file in the folder
#!/usr/bin/env python

from pyrap.tables import table
import sys
import os

def cal_expand(calname, freq=100e6, time=86400):
    with table(calname, ack=False, readonly=False) as caltab:
        caltab.putcol("STARTX", caltab.getcol("STARTX") - freq)
        caltab.putcol("ENDX", caltab.getcol("ENDX") + freq)
        caltab.putcol("STARTY", caltab.getcol("STARTY") - time)
        caltab.putcol("ENDY", caltab.getcol("ENDY") + time)

def main():
    for calin in sys.argv[1:]:
        calout = calin +'.timeext'
        os.system('cp -r '+ calin + ' ' + calout)
        cal_expand(calout)

if __name__ == "__main__":
    main()

