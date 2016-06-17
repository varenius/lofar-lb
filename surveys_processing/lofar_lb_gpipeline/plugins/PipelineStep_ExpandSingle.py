import os
import numpy as np
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


# Colm Coughlan, March 2016
# used in locapi generic pipeline implementation


def plugin_main(args, **kwargs):
    """
    Takes in mapfile_in, containing many files, and returns only one

    Parameters
    ----------
    mapfile_in : str
        Parmdbs containing phase solutions
    mapfile_dir : str
        mapfile directory
    filename : str
		output filename
    n: str
        Size to expand to
    mapfile_comp : str
		target MSs

    Returns
    -------
    result : dict
        Output datamap filename

    """
    mapfile_dir = kwargs['mapfile_dir']
    mapfile_in = kwargs['mapfile_in']
    mapfile_comp = kwargs['mapfile_comp']
    filename = kwargs['filename']
    n = int(kwargs['n_expand'])

    tables_data = DataMap.load(mapfile_in)		# these are the phase solution tables
    ms_data = DataMap.load(mapfile_comp)	# these are actual MS files
    
    iterator = DataMap.SkipIterator(ms_data) 
    ms_list = []
    for value in iterator:
		name = value.file
		name = (name.split('SB'))[1]
		name = name[0:3]
		ms_list.append(int(name))
    
    iterator = DataMap.SkipIterator(tables_data)
    sol_list = []
    for value in iterator:
		name = value.file
		try:
			name = (name.split('SBgr'))[1]	# try for SB groups first
		except:
			name = (name.split('SB'))[1]
		name = name[0:3]
		sol_list.append(int(name))

    ms_list = np.array(ms_list)
    sol_list = np.array(sol_list)
    
    if( (len(ms_list) == 0) or (len(sol_list) == 0)):
        print('ERROR: ExpandSingle: 0 solution files or MS files detected.')
        print('ERROR: ExpandSingle: # MS files = '+str(len(ms_list)))
        print('ERROR: ExpandSingle: # Solution files = '+str(len(sol_list)))
        return(1)
    
    diff_list = np.diff(sol_list)
    if(len(diff_list)>0):
        n_internal = np.min()
        if(n_internal != n):
            print('WARNING: Internal estimate of n differs from user-provided value.')
            print('WARNING: Internal estimate will be used.')
    else:
        n_internal = n
		
	
		
	# typically the sol name

    map_out = DataMap([])
    nsol = 0
    nms = 0
    iterator = DataMap.SkipIterator(tables_data)
    for value in iterator:
		for i in range(n_internal):
			if(ms_list[nms] - sol_list[nsol] > n_internal):
				print('Error in expanding phase calibration mapfile.')
				return(1)
			else:
				map_out.data.append(DataProduct(value.host,value.file, value.skip ))
			nms = nms + 1
			if(nms == len(ms_list)):
				break
		nsol = nsol + 1

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
