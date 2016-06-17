#################################################################################
#                                                                               #
# Written by Leah Morabito, 15 Mar 2016                                         #
#                                                                               #
#	Add international stations to a parmdb					#
#	Assumptions: <name>.table is the parmdb for <name>.ms			#
#                                                                               #
#################################################################################

import lofar.parmdb as lp
import pyrap.tables as pt
import numpy as np
import sys

if len(sys.argv) < 2:
        print 'Please give a parmdb name.'
        print 'Usage: python testgain.py <parmdbname>'

filename = sys.argv[1]
msfile = filename.rstrip('/').replace('.table','.ms')

## get list of antennas from the msfile
ant_table = pt.table(msfile+'/ANTENNA')
ant_names = ant_table.getcol('NAME')

## get list of international stations
int_stations = [ ant for ant in ant_names if 'CS' not in ant and 'RS' not in ant ]

## open the parmdb and read in the dictionary
parmdbmtable = lp.parmdb(filename)
dictionary = parmdbmtable.getValuesGrid('*')
dict_keys = dictionary.keys()

## get an antenna dictionary to copy
tmp_dict = dictionary[dict_keys[0]]

real_dict = tmp_dict.copy()
real_dict['values'] = np.ones(tmp_dict['values'].shape)
imag_dict = tmp_dict.copy()
imag_dict['values'] = np.zeros(tmp_dict['values'].shape)

for int_station in int_stations:
	dictionary['Gain:0:0:Real:'+int_station] = real_dict
	dictionary['Gain:1:1:Real:'+int_station] = real_dict
	dictionary['Gain:0:0:Imag:'+int_station] = imag_dict
	dictionary['Gain:1:1:Imag:'+int_station] = imag_dict
	

parmdbmtable.deleteValues('*')
parmdbmtable.addValues(dictionary)
parmdbmtable.flush()

