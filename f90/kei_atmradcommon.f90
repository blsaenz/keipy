!  FROM: -------------------- atmrad.com

! common arguments for "init atm" and "atmstep" and "prtprof"

    integer, save :: nrow             !  latitude row index
    real, save :: clat                !  Current latitude (radians)
    real, save :: oro(plond)          !  land/ocean/sea ice flag
    real, save :: sndpth(plond)       !  snow depth (liquid water equivalent)
    real, save :: ps(plond)           !  surface pressure
    real, save :: pmid(plond,plev)    !  model level pressures
    real, save :: pint(plond,plevp)   !  model interface pressures
    real, save :: pmln(plond,plev)    !  natural log of pmid
    real, save :: piln(plond,plevp)   !  natural log of pint
    real, save :: t(plond,plev)       !  model level temperatures
    real, save :: h2ommr(plond,plev)  !  model level specific humidity
    real, save :: o3mmr(plond,plev)   !  ozone mass mixing ratio
    real, save :: cldfrc(plond,plevp) !  fractional cloud cover
    real, save :: effcld(plond,plevp) !  effective fractional cloud cover
    real, save :: clwp(plond,plev)    !  cloud liquid water path
    real, save :: plol(plond,plevp)   !  o3 pressure weighted path lengths (cm)
    real, save :: plos(plond,plevp)   !  o3 path lengths (cm)
    real, save :: tstep               !  time step in seconds
    integer, save :: numavg           ! number of iterations for averaging
    real, save :: clon(plon)          !  longitude        (radians)
    real, save :: tIN(plond,plev)     !  model level temperatures   from INPUT
    real, save :: h2ommrIN(plond,plev)! model level spec. humidity from INPUT
    character (len=80), save :: &
			label,     & 										!  labels from INPUT
			tlabel, &
			dlabel


