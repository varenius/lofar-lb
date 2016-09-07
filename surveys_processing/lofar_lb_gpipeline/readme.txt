This is a generic pipeline implementation of the LOFAR long baseline reduction pipline
At the moment it produces FITS files which can then be read into AIPS

A lot of code (especially plugins) has been "borrowed" from Prefactor

The pipeline can be used in two ways:

1. It can perform standard amplitude calibration using a calibrator source, similar to previous versions of locapi

2. It can use the output from the prefactor pipeline - both amplitude, clock and direction independent phase solutions can be applied to the target field.



Depending on your choice, you will need to fill out some parameters at the top of the parset file. You may also choose to provide a manual list of target sources (see the included file for the format of such a list), or you can use the LBCS (formally LOBOS) survey to return a list of good candidates.

Regardless of the means of selection, the dataset will be shifted onto each source in turn before averaging, flagging and forming a tied station (the use can choose which stations to phase up). The data are then converted to circular polarisation and written out as a FITS file.



Done:

Add flagging
Change Adam's plugin to ensure the original phase calibration instrument tables are not overwritten (write new ones to working directory)
Add in standard locapi amp calibration functions


To do:

Add Stephen's parseltongue fring function
