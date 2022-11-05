!     ==================================================================
!     KPP-Ecosystem-Ice (KEI, pronounced 'key') Model
!     ==================================================================
!
!     Version History (please add a version modification line below
!     if code is edited)
!     ------------------------------------------------------------------
!     Version: 0.9 (2022-01-05, Benjamin Saenz, blsaenz@gmail.com)
!
!     This model derives from Large et al. [1994], Doney et al. [1996](KPP mixing),
!     Ukita and Martinson [2001](Mixed layer - ice interactions), Saenz and Arrigo
!     [2012 and 2014] (SIESTA sea ice model), Hunke and Lipscomb [2008] (CICE v4),
!     Moore et al. [2002,2004] (CESM ecosystem model).
!
!     DEPENDENCIES:  1) A FORTRAN 90 compiler (gfortran and ifort have been tested)
!                    2) CPP pre-processor
!                    3) NetCDF 3.6.1 or higher, with FORTRAN 90 interface built
!                    4) LAPACK (used by SIESTA for solving heat transfer)
!                    5) MATLAB (not required to run compile or run KEI, but
!                        input data generation, run scripts, and analysis tools
!                        that greatly facilitate use of KEI are MATLAB-based)

RUNNING KEI: To run KEI you need:
  1) Compiled KEI.run executable file (linked to NetCDF and LAPACK libraries)
  2) NetCDF forcing data file
  3) Text-based run options file

COMPILATION:

  1) Edit the file "makefile" to suit your system.  Setting the FC variable
  to either ifort for gfortran will give you independent control over those
  two compilers.  You will have to edit the paths to the NetCDF include and
  lib directories, and also amend the linking to some version of LAPACK.
  (Note that ifort and MacOS X ship with LAPACK - you must build the LAPACK
  fortran interfaces before using however). Windows is not supported, mostly
  because it is incredibly difficult to get the fortran NetCDF interface
  installed on windows.

  2) Type "make clean", then "make."  The executable KEI.run should be produced.
  Various warnings are produced by different compilers - they seen to be
  benign on my test systems thus far.

PREPARING FORCING DATA:

  1)  Examine the example NetCDF forcing data: kf_200_100_2000.nc
  This file corresponds to several years of hourly forcing from the
  200.100 grid station of the Palmer LTER project, on the continental
  shelf west of the Antarctic Penninsula.  The variables must be named
  exactly the same as in this file.  The resolution of the forcing and
  initialization data must match the grid (i.e. KEI does not spatially
  interpolate data).

  You can use they NetCDF tool of your choice to edit the NetCDF data,
  or (likely much easier) you can use MATLAB to load and edit the
  corresponding kf_200_100_2000.mat file, which contains the same forcing
  data.  To write the MATLAB forcing data structure to NetCDF for use
  by KEI.run, use the kei_write_forcing MATLAB function.

CREATING THE RUN OPTIONS FILE:

  The example file is called "run.so."  It is rather long, however
  much of the functionality at this point is legacy and either unused
  or untested.  (I hesitate to remove before I know it's truly useless...)
  The important lines are:

  Line 13:  start time (in days) - corresponds to the "time" forcing variable
  Line 14:  time step start (always 0), end (# desired steps), and
    time step length (seconds - 1hr (3600s) timestep is the only one tested!).
  Line 21:  The PATH where output data will be written.
  Line 24:  Toggles for switching on/off parts of the model, such as ice,
    ecosystem, KPP, etc.  I have only tested toggling the ecosystem and ice.
  Line 50:  Path/filename to the forcing NetCDF file is located

RUNNING THE KEI MODEL:

  To run showing way to much debugging information, where run.so is the run options file:
    KEI.run < run.so
  To run w/out showing debugging information, where run.so is the run options file:
    KEI.run < run.so > /dev/null

POSTPROCESSING:

    A number of matlab routines are available for examining and plotting
    data from the output NetCDF file.  By default, the single output
    file contains data from each time step and can be quite large, ~800Mb
    for a year of run time.

    PLEASE NOTE !!!:
    ------------------------------------------------------------------
    Many of the matlab commands have been hardwired to process
    data written at the hourly level.  If your time step or write step
    differers, please examine the plotting and postprocessing commands!

    DATA READING:
    kei_read:  Reads a large fraction of the output variables into memory.
    This will take a long time on a slow disk drive, or if your output
    file is big.
    kei_read_fluxes:  Reads a subset of the physical output data

    SIMPLE PLOTTING:
    kei_plot:  simple plotting of kei physical ocean data
    kei_plot_ice: lots of time series plots showing water and ice through time
    kei_plot_eco:  tons of time series plots of water, ice, biology

    NICER PLOTS:  (some of these require the forcing data structure too)
    kei_plot_ecoline
    kei_plot_fluxes
    kei_plot_ecoiceline


