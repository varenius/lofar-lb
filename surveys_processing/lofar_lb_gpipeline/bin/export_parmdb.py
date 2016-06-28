#!/usr/bin/env python
# -*- coding: utf-8 -*-

# converted from locapiv3b.py by Colm Coughlan, March 2016

import os

# locapi script to transfer gain solutions from flux calibrator to phase calibrator

if __name__=='__main__':
    # Options
    print('Sorry, only available as a python-plugin currently...')

def main(inname, sample_ms, add_IS_script, correct_IS_script, wd, outsuffix='solutions.table'):

    table_plus_is = inname.split('/')[-2] + '/' + (inname.split('/')[-1]) + '_plusIS'	# get a good name for the new table (located inside the current MS)
    outname = (inname.split('/')[-2]).split('.')[0] + '.' + outsuffix			# get a good name for the final table (located in the working directory)
    outname = os.path.join(wd,outname)

    os.system('cp -r '+inname+' '+table_plus_is)
    os.system(add_IS_script + ' ' + table_plus_is + ' ' + sample_ms + ' 0')
    os.system(correct_IS_script + ' ' + table_plus_is)
    os.system('parmexportcal in='+table_plus_is+' out='+outname)	# use parmdbm to filter gains and remove time dependence



    return {'parmdbfile': outname }	# return name of final table file created to main pipeline for further use


