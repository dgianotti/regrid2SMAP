# regrid2SMAP
A library for regridding data sets to the SMAP 36km EASE-Grid 2.0 standard

The MATLAB scripts and functions in this library are to be used to regrid various Earth-system data products to the 36km EASE-Grid 2.0 grid used by the SMAP (Soil Moisture Active/Passive) mission.

Calculation of the sparse regridding matrices is fairly time-consuming (on the order of minutes/hours), and so this library also uses some saved MATLAB binary (.mat) files to reduce the difficult computational phase.


