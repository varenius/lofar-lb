#!/usr/bin/python
# A very simple script for converting to circular polarisation using TAQL. 
# This conversion does NOT include the element beam model.
import pyrap.tables as pt

oldms = '/data/scratch/adeller/L50631/quick/L50631_SAP000_SB120_uv.MS'
newms = '/data/scratch/adeller/L50631/quick/L50631_SAP000_SB120_uv.circ.MS'

# Copy to a new MS and update to circular
t = pt.table(oldms)
t.copy (newms)
t.close()
t=pt.taql ("update %s set DATA = mscal.stokes(DATA,'circ')" % newms)

# Update the polarization definition (the values represent RR,RL,LR,LL)
t=pt.taql (" update %s/POLARIZATION set CORR_TYPE=[5,6,7,8]" % newms)
t.close()
