import os
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

    value = DataMap.load(mapfile_in)[0]		# this the the single mapfile to be expanded
    n = len(DataMap.load(mapfile_comp))	# these are actual MS files

    map_out = DataMap([])
    for i in range(n):
        map_out.data.append(DataProduct(value.host,value.file, value.skip ))

    fileid = os.path.join(mapfile_dir, filename)
    map_out.save(fileid)
    result = {'mapfile': fileid}

    return result
