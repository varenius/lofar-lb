#!/usr/bin/env python
# This script can be used to plot the parmdb-solutions for several parmdbs as
# function of time.  This is useful e.g. to plot the phase-solutions derived
# for core stations, stored as multiple parmdb folders, usually one per
# phase-up scan. So, every 20-30 minutes there is a phase solution found for the
# CS using e.g. NDPPP. See plotparmdbsols_example.png for example output.
# This script originally written by Stephen Bourke.

from pyrap.tables import table
import sys
from glob import glob
import matplotlib.pyplot as plt
from tempfile import mkdtemp
from subprocess import check_call
from shutil import rmtree
import numpy as np

NUM_CORE_STATIONS = 24 * 2 # Dual HBA mode
FULL_RANGE = True # Plot from -180 to 180

def get_vals(cal_sets, row_num):
    vals = np.array([table(t, ack=False).getcell("VALUES", row_num)[0,0] for t in cal_sets]) / np.pi * 180
    names = [table(t+"::NAMES", ack=False).getcell("NAME", row_num) for t in cal_sets]
    assert len(set(names)) == 1, "Antennas vary among cal tables: " + str(names)
    return names[0].split(':')[-1], vals

def plot2(cal_sets, ant_num, out_dir='.'):
    plt.clf()
    p_name, p = get_vals(cal_sets, ant_num * 2)
    q_name, q = get_vals(cal_sets, ant_num * 2 + 1)
    assert p_name == q_name, "Expected two pols from the same ant, got: " + str([p_name, q_name])
    if FULL_RANGE:
        plt.ylim([-180,180])
    plt.plot(p, 'o')
    plt.plot(q, 'o')
    plt.title(p_name)
    outname = "{dir}/ant{n:02d}.png".format(dir=out_dir, n=ant_num)
    plt.savefig(outname)
    return outname

def main():
    if len(sys.argv) < 3:
        print >> sys.stderr, "Usage: plotparmdbsols.py <cal1> <cal2> .... <caln> <outplot>"
        print >> sys.stderr, "e.g. plotparmdbsols.py calscan1.parmdb calscan2.parmdb calscan3.pardb calscans.png"
        sys.exit()
    cal_sets = sys.argv[1:-1]
    out_plot = sys.argv[-1]
    tmpd = mkdtemp()
    try:
        plots = [plot2(cal_sets, i, out_dir=tmpd) for i in range(NUM_CORE_STATIONS)]
        check_call("montage -tile 4x12 -geometry 400x400".split() + plots + [out_plot])
    finally:
        rmtree(tmpd)

if __name__ == "__main__":
    main()

