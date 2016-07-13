This is a generic pipeline implementation of the LOFAR long baseline reduction pipline
At the moment it produces FITS files which can then be read into AIPS

A lot of code (especially plugins) has been "borrowed" from Prefactor



Done:

Add flagging
Change Adam's plugin to ensure the original phase calibration instrument tables are not overwritten (write new ones to working directory)
Add in standard locapi amp calibration functions


To do:

Add Stephen's parseltongue fring function
