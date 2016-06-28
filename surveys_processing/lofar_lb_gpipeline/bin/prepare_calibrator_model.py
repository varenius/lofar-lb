#!/usr/bin/env python

# converted from locapiv3b.py by Colm Coughlan, March 2016

import os

def main(dummy, model, wd, makesourcedb_path):

    outfile = os.path.join(wd,'fluxcal_model.mod')
    
    os.system(makesourcedb_path + ' in='+model+' out='+outfile+' format=\"<\"') 

    return {'modelfile': outfile }	# return name of final file created to main pipeline for further use
    
if __name__=='__main__':
    # Options
    print('Only a pythonplugin...')

