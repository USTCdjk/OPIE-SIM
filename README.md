My_SIM_main. m is the main program. After opening it, you can locate the SIM dataset location by modifying the 'filepath'. The parameters that need to be modified are:
rawpixelsize = [65 65 100]; %  Camera pixel size/imaging system magnification
NA=1.49;% objective NA
Exwavelength=405;% excitation wavelength 488 561 638 668
Emwavelength=525;% emitted light wavelength 525 610 665 680
Numsteps=3;% number of sample phases
NrBands=2;% coherent light quantity, for example, 2DSIM, S (k) and S (k ± p) bands are two, 3D, and three
Numangles=3;% number of structural light directions
Readmode=2;% Data format: 1. paz 2. pza (OMX)
Among them, emwavelength will cause the algorithm to use different correction 3D-OTFs. If a special OTF is required, it can be generated by the PSF generator itself