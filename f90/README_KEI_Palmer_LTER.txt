SPECIAL KEI RUN INSTRUCTIONS - PALMER LTER
========================================================================

A lot of routines and forcing data are already compiled for the Palmer
LTER regions that facilitate running and plotting data.  Below are 
instructions for run two different types of simulations.

NOTE!!! - To get matlab to run gfortran binaries on a Mac, you
may have to rename /Applications/Matlab_XXXXX.ap/sys/os/maxi64/libgfortran.3.dylib
to something different, and then link your real copy of libgfortran.3.dylib
to this location.  The MATLAB supplied libgfortran is totally useless.

Palmer LTER mooring-forced runs
------------------------------------------------------------------------
These runs occur at three different mooring locations, 300.100, 300.160,
and 400.100, during different year during 2007-2011.  Atmospheric forcing 
data comes from ECWMF Interim at these exact locations.  Initial profiles 
come from CTD casts from the LTER , excepting 300.160 in 2010 and 2011
when no CTD was taken at these locations. The ecosystem initialization is
the same in each set of forcing data.

Available runs -------------------------
Location    Year    Forcing Data File   Run Options File
300.100     2007    kf_300_07_ECMWF.nc  run.07.300.so
300.100     2008    kf_300_08_ECMWF.nc  run.08.300.so
300.100     2010    kf_300_10_ECMWF.nc  run.10.300.so
300.100     2011    kf_300_11_ECMWF.nc  run.11.300.so
300.160     2008    kf_360_08_ECMWF.nc  run.08.360.so
300.160     2009    kf_360_09_ECMWF.nc  run.09.360.so
300.160     2010    kf_360_10_ECMWF.nc  run.10.360.so
300.160     2011    kf_360_11_ECMWF.nc  run.11.360.so
400.100     2008    kf_400_08_ECMWF.nc  run.08.400.so
400.100     2009    kf_400_09_ECMWF.nc  run.09.400.so
400.100     2010    kf_400_10_ECMWF.nc  run.10.400.so

Steps to run and plot:
1) compile KEI.run
2) create a final data directory, and change to that directory in matlab
3) open kei_moorings.m for editing
Change the code, data, and exe paths to where you are compiling, writing
temporary output, and running the executable respectively. (In many cases,
it is faster to write to a SSD drive initially [i.e. use this as you DATA
path and keep the forcing files there] The final data an plots will be 
moved to the current matlab directory you just created).
4) in matlab run the command:
  kei_moorings([2007:2011],'run_string',1,[0,1,2]);
  
  This will run and plot all the available mooring runs.  To run a subset,
  only enter the years you want in the 1st input variable, and/or change
  the last input variable to limit the number of sites ([0 = 300.100, 
  1 = 400.100, 2 = 300.160]).  Change run_string to label the run
  output files and plots.


  
200.100 yearly runs
------------------------------------------------------------------------
These runs use ECMWF Interim forcing data taken from 200.100, from the years
1997 to 2010.  These runs are intended for ecologocal exploration,
and correspond to the SEAWIFS/MODIS years for ocean chl a comparison, as
well as the Palmer LTER dataset.  Initial profiles come from CTD casts 
from the LTER, and ecosystem initialization is the same in each set of 
forcing data.


Available runs -------------------------
Location    Year    Forcing Data File   Run Options File
200.100     1997    kf_200_97_ECMWF.nc  run.97.so
200.100     1998    kf_200_98_ECMWF.nc  run.98.so
200.100     1999    kf_200_99_ECMWF.nc  run.99.so
200.100     2000    kf_200_100_2000.nc  run.00.so
200.100     2001    kf_200_01_ECMWF.nc  run.01.so
200.100     2002    kf_200_02_ECMWF.nc  run.02.so
200.100     2003    kf_200_03.nc        run.03.so
200.100     2004    kf_200_04_ECMWF.nc  run.04.so
200.100     2005    kf_200_05_ECMWF.nc  run.05.so
200.100     2006    kf_200_06_ECMWF.nc  run.06.so
200.100     2007    kf_200_07_ECMWF.nc  run.07.so
200.100     2008    kf_200_08_ECMWF.nc  run.08.so
200.100     2009    kf_200_09_ECMWF.nc  run.09.so
200.100     2010    kf_200_10_ECMWF.nc  run.10.so


Steps to run and plot:
1) compile KEI.run
2) create a final data directory, and change to that directory in matlab
3) open kei.m for editing
Change the code and data paths to where you are compiling, writing
temporary output, and running the executable respectively. (In many cases,
it is faster to write to a SSD drive initially [i.e. use this as you DATA
path and keep the forcing files there] The final data an plots will be 
moved to the current matlab directory you just created).
4) in matlab run the command:
  kei([1997:2010],'run_string',1);
  
  This will run and plot all the available 200.100 runs.  To run a subset,
  only enter the years you want in the 1st input variable.  Change run_string 
  to label the run output files and plots.
