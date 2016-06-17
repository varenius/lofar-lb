import os
import numpy as np
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


# Colm Coughlan, March 2016
# used in locapi generic pipeline implementation

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx]):
        return array[idx-1]
    else:
        return array[idx]


def plugin_main(args, **kwargs):
    """
    Takes in mapfile_in, containing many files, and returns only one

    Parameters
    ----------
    mapfile_in : str
        Name of the input mapfile to be trimmed back.
    mapfile_dir : str
        Directory for output mapfile
    filename: str
        Name of output mapfile

    Returns
    -------
    result : dict
        Output datamap filename

    """
    mapfile_dir = kwargs['mapfile_dir']
    mapfile_in = kwargs['mapfile_ms']
    mapfile_grpd = kwargs['mapfile_grpd']
    filename = kwargs['filename']

    data = DataMap.load(mapfile_in)		# these are actual MS files
    groups = DataMap.load(mapfile_grpd)	# these are probably parmdbs
    
    name_stem = ((data[0].file.split('SB'))[0]) + 'SB'
    name_end =  '_'.join(((data[0].file.split('SB'))[1]).split('_'))[3:]

    iterator = DataMap.SkipIterator(data)
    num_list = []
    for value in iterator:
		name = value.file
		name = (name.split('SB'))[1]
		name = (name.split('_'))[0]
		name = (name.split('.'))[0]
		num_list.append(int(name))
  
    iterator = DataMap.SkipIterator(groups) 

    map_out = DataMap([])
    for value in iterator:
		name = value.file
		name = (name.split('SB'))[1]
		name = (name.split('_'))[0]
		name = (name.split('.'))[0]
		name = (name.split('gr'))[-1]	# in case of grouping!
		num = int(name)
		if num in num_list:
			map_out.data.append(DataProduct(value.host, name_stem + name + name_end, value.skip ))
		else:
			num = find_nearest(num_list, num)	# if not an exact match, use the nearest
			map_out.data.append(DataProduct(value.host, name_stem + str(num) + name_end, value.skip ))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
