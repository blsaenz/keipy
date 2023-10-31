  module kei_icecommon

    !use kei_parameters
    ! f2py won't allow this, so specified in two places ... make sure nni = z_max_ice, nns/nnfs = z_max_snow !!
    ! use sia2_parameters, only: z_max_ice, z_max_snow

    implicit none

    integer, parameter  :: &
      nni            =  42,      &  !
      nnfs           =  26,      & !
      nns            =  26        !
    real, parameter :: fs_crit = -0.02 ! height below freeboard at which flooding occurs (m)


    !integer, parameter  :: nnis=nns+nnfs+nni
    integer, parameter  :: nnis=nns+nni

    ! ice advection
    real, parameter :: export_melt_f = 0.  ! fraction of ice exported that is assumed to melt nearby/influence salinity
    real, save :: import_melt_f ! multiplier of calculated ice melt - affects ocean salinity
    real, save :: ice_fe_local ! storage for what concentration of fe the ice has


    integer, save  :: &
      ni_cur,   &
      ns_cur

    real, save  :: &
      Tas,      &
      Ti(nni),  &
      Si(nni),  &
      Ts(nns),  &
      Tfs(nnfs),  &
      dzfs(nnfs), &
      dzs(nns), &
      dzi(nni)

    ! timice
    integer, save :: &
      ndtice,icemax,inew,iold
    real, save :: &
      runice,dti

    ! ice state
    real, save :: &
      fice,hice(0:1),TI0,qI0,uI,vI,albice

    ! ice paras
    real, save :: &
      dsice,tfscl,SWFACS,tmlt,tfrz,sice,rhoice, &
      CPice,epsi,sfrazil,fCoriolis,sslush

    ! ice force
    integer, save :: &
      NSICE
    real, save :: &
      FROCN,FRICE,R1,R0,DHDTice,dhs,dhi,dhfs,shs

    ! snow state
    real, save :: &
      hsn(0:1),rhosn,rhofsnow,CPsn,CPfsnow, &
      hfsnow(0:1)

    ! flx paras
    real, save :: &
      EL,SL,FL,FLSN,C2K,Qsfc

    real, save :: &
      lateral_freshwater_melt_flux

    real, save :: atm_flux_to_ice_surface, &
      ice_ocean_bottom_flux, &
      ice_ocean_bottom_flux_potential, &
      total_ice_melt, &
      total_ice_freeze, &
      frazil_ice_volume, &
      congelation_ice_volume, &
      snow_ice_volume, &
      snow_precip_mass



!    type, public :: flux_type
!      integer :: &
!          jptr
!      real :: &
!          o_latent, &            ! ocean latent
!          o_sensible, &          ! ocean latent
!          o_shortwave, &         ! ocean latent
!          o_longwave, &          ! ocean latent
!          o_ice, &          ! ocean latent
!          i_latent, &            ! ocean latent
!          i_sensible, &          ! ocean latent
!          i_shortwave, &         ! ocean latent
!          i_longwave, &          ! ocean latent
!          i_ocean, &          ! ocean latent
!          ice_heat
!      end type flux_type

  end module kei_icecommon

