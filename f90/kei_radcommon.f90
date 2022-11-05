module kei_radcommon

!    common-blocks for radiative/convective model
!    includes :
!    comtim, drnl, global, prt, crdcon, crdcae, crdctl, crdalb

!     !   Basic grid point resolution parameters
		implicit none

    integer, parameter :: plev    = 18      !  number of vertical levels
    integer, parameter :: plevp   = plev+1  !  number of vertical interfaces
    integer, parameter :: plon    = 1       !  number of longitudes (T42)
    integer, parameter :: plat    = 1       !  number of latitudes (T42)
    integer, parameter :: plond   = 1       !  slt extended domain longitude
    integer, parameter :: psave   = 10      ! number of extra rad/conv variables
!                                    to be stored

!               plond= plon + 1 + 2*nxpt, !slt extended domain longitude

! model control time variables ********************************** comtim

    DOUBLE PRECISION, save :: calday ! an change calday to double precision
    integer, save :: nstep

! diurnal cycle switch; if true, does diurnal cycle *************** drnl

    logical, save :: diurnal

! global mean switch; if true, does global mean computation ***** global

    ! common /glbl/ global
    logical, save :: globl ! changes from global, which is keywork in fortran 90 - saenz 7/2011
!
! clear sky switch; if true, does diagn. clear sky comp. ****** clearsky

    logical, save :: clrsky
!
! radiation constants ******************************************* crdcon

    real, save :: gravit,   & !  gravitational acceleration
    rga,   & !  1 over gravit
    cpair,   & !  heat capacity air at constant pressure
    epsilo,   & !  ratio mmw h2o to mmw air
    sslp,   & !  standard pressure
    stebol,   & !  stephan boltzmann constant
    rgsslp,   & !  0.5 / (gravit*sslp)
    co2vmr,   & !  co2 volume mixing ratio
    dpfo3,   & !  Doppler factor for o3
    dpfco2,   & !  Doppler factor for co2
    dayspy,   & !  solar days in one year
    pie    ! pie


! water vapor narrow band constants for lw computations ********* crdcae

    real, dimension(2), save :: realk,st,a1,a2,b1,b2

! constant coefficients for water vapor absorptivity and emissivi

    real, dimension(3,4), save :: coefa, coefc, coefe
    real, dimension(4,4), save :: coefb, coefd
    real, dimension(6,2), save :: coeff, coefi
    real, dimension(2,4), save :: coefg, coefh
    real, dimension(3,2), save :: coefj, coefk
    real, dimension(4), save :: c1,c2,c3,c4,c5,c6,c7
    real, save :: c8 ,c9 ,c10,c11,c12,c13,c14,c15,c16,c17, &
    	c18,c19,c20,c21,c22,c23,c24,c25,c26,c27, &
    	c28,c29,c30,c31

! farwing correction constants for narrow-band emissivity model
! introduce farwing correction to account for the
! deficiencies in narrow-band model used to derive the
! emissivity. tuned with arkings line calculations.

        real, save :: fwcoef,fwc1,fwc2,fc1,cfa1


! radiation control variables *********************************** crdctl

! fradsw = .t. iff full shortwave computation
! fradlw = .t. iff full longwave computation

! irad = iteration frequency for radiation computation
! iradae = iteration frequency for absorptivity/
! emissivity computation


    integer, save :: iradae,irad,naclw,nacsw,fnlw,fnsw
    logical, save :: aeres


! surface albedo data ******************************************* crdalb

! vs = 0.2 - 0.7 micro-meters wavelength range
! ni = 0.7 - 5.0 micro-meters wavelength range

! s  = strong zenith angle dependent surfaces
! w  = weak   zenith angle dependent surfaces

! the albedos are computed for a model grid box by ascribing values to
! high resolution points from a vegetation dataset, then linearlly
! averaging to the grid box value; ocean and land values are averaged
! together along coastlines; the fraction of every grid box that has
! strong zenith angle dependence is included also.

    real, dimension(plond,plat), save :: & 
    albvss, & !  grid box alb for vs over strng zn srfs
    albvsw, & !  grid box alb for vs over weak  zn srfs
    albnis, & !  grid box alb for ni over strng zn srfs
    albniw, & !  grid box alb for ni over weak  zn srfs
    frctst  ! fraction of area in grid box strng zn

! surface boundary data

! rghnss is the aerodynamic roughness length for the grid box, computed
! by linear averaging of the values ascribed to high resolution
! vegetation dataset values; ocean and land values are averaged together
! at coastlines.

    real, dimension(plond,plat), save :: &
    rghnss  ! aerodynamic roughness length


! FROM: -------------------- atmrad.com

! common arguments for "init_atm" and "atmstep" and "prtprof"

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



end module kei_radcommon


