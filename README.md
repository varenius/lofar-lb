lofar-lb
========

Collection of LOFAR long baseline software.

-- Installation--
Clone as ordinary git repo. "mscorpol" may need special attention, see below. 

- Getting mscorpol-
"mscorpol", a tool written by T.D. Carozzi to convert measurement sets from
linear to circular polarisation, is included as a "git submodule".  When cloning
the lofar-lb repository, you may get an empty mscorpol folder. In that case,
you need to run "git submodule init" and then "git submodule update" to fetch
all data from the mscorpol submodule as well.
