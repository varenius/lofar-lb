#!/usr/bin/env python

# converted from locapiv3b.py by Colm Coughlan, March 2016

import os
import sys

# short locapi script to run mscorpol on current file
# (mscorpol doesn't have a main function suitable for direct pipeline integration)
# can also use other convertors

def main(msname, cpath, do_mscorpol=0):

	do_mscorpol = int(do_mscorpol)
	if do_mscorpol == 1:
		os.system('python '+path+' -f '+msname)	# execute mscorpol command
	else:
		print('python '+cpath+' '+msname)
		os.system('python '+cpath+' '+msname)
    
if __name__ == "__main__":
    # Options
    print('Mostly a pythonplugin, but attempting a run anyway...')
    print('Input form: msname, path to converter, boolean: true if using mscorpol')
    if(len(sys.argv)==4):
		print(sys.argv[3])		
		main(sys.argv[1], sys.argv[2],sys.argv[3])
    else:
		print('Error: check your inputs.')


