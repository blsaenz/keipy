MODULE kei_hacks

  USE kei_ecocommon, ONLY: log_kind,int_kind,real_kind,dbl_kind

  IMPLICIT NONE

  ! overwrite in some eddy heat at depth during certain periods
  LOGICAL, PARAMETER :: assimilate_ts_hack = .true.   ! assimilate deep t profile (LTER mooring data)
  LOGICAL, PARAMETER :: assimilate_ts_stop = .false.   ! only assimilate deep t profile for 1st X steps simulation
  LOGICAL, PARAMETER :: assimilate_ts_fixed = .false.  ! restore deep t profile to a single fixed profile
  LOGICAL, PARAMETER :: eddy_hack = .false.            ! use this to do some temporary deep heating - this is sort of a old hack, need to carefully figure out what it does

  ! fraction precipitation adjustment OVER ICE ONLY
  REAL, PARAMETER :: snow_fraction = 1.00  ! multiply snow precip by this
  REAL, PARAMETER :: rain_fraction = 0.0   ! multiply rain precip by this

  ! fraction shortwave irradiance
  REAL, PARAMETER :: shortwave_multiplier = 0.80  ! multiply shortwave irradiance by this

  ! conform to ice concentration found in forcing (0=no, 1=yes)
  INTEGER, PARAMETER :: ic_conform = 0         ! conform to ic forcing observations
  INTEGER, PARAMETER :: ignore_divergence = 1  ! don't raft up ice
  INTEGER, PARAMETER :: zero_ice_before = -1

  ! lateral ice growth/melt hacks
  LOGICAL, PARAMETER :: melt_lateral_complete = .true.  ! use all extra available mixed layer heat to melt laterally (.true. is standard)
  LOGICAL, PARAMETER :: mixed_layer_freezing_in_leads_only  = .true.  ! direct all mixed layer freezing flux to new ice growth (.true. is standard)

  ! new ice thickness
  REAL, PARAMETER :: initial_ice_thickness = 0.30

  ! fixed ctd profiles for assimilation
  REAL, ALLOCATABLE, SAVE :: fixed_t(:), fixed_s(:)

  INTEGER :: assim_counter
  INTEGER :: assimilate_step_limit ! number of time steps after which ts assimilation is stopped

!-----------------------------------------------------------------------
!     initial fe/bio profile adjustments
!-----------------------------------------------------------------------
  LOGICAL, PARAMETER :: assimilate_fe_from_ts = .true.

  real (kind=dbl_kind), parameter :: &
    fe_multiplier = 1.0, &
    fe_offset = 0., &    ! negative to use pre-estimated pycnocline fe levels below, positive for m lower
    bio_offset = 0., &    ! m lower
    alk_spike_value = 3000.0, &
    meltwater_fe = 0.0      ! nm
!      meltwater_fe = 12.0e-9
  INTEGER, PARAMETER :: &
    alk_spike_timestep = 480+4344, &
    alk_spike_level_start = 0, &
    alk_spike_level_end = 8

!-----------------------------------------------------------------------
!     sea ice fe/algal profile concentrations during ice melt
!-----------------------------------------------------------------------

    LOGICAL, PARAMETER :: use_year_ice_advance = .false.

    ! Fixed ice chla & Fe concentrations, for dumping into ocean during melt
    real (kind=dbl_kind), parameter :: &
      ice_diatChl = 0.  , &      ! mgchla/m^2 (meiners et al 2012) / 0.6 m = mg/m^3
!      ice_diatChl = 25.6  / 0.6, &     ! mg chla/m^2 (meiners et al 2012) / 0.6 m = mmol/m^3
!      ice_fe = 6.              ! nm
      ice_fe = 0.              ! nm

!-----------------------------------------------------------------------
!     variables used in adding a static krill grazer (not yet implemented)
!-----------------------------------------------------------------------

    logical, parameter :: &
      use_krill_hack = .false.

    real (kind=dbl_kind), parameter :: &
      krillC = 2000.          , &! krill concentration in carbon (mg C m-3)
      krill_sp_umax = 1.0     , &
      krill_diat_umax = 4.0  , &
      krill_saturation_grazing = 2.0 ! max out krill grazing ability - mG Chl a / m^3

  CONTAINS

    SUBROUTINE hacks_init()
        assim_counter = 0
        assimilate_step_limit = find_assimilate_step_limit()
    END SUBROUTINE hacks_init

    INTEGER FUNCTION find_assimilate_step_limit()

      USE kei_common

      IMPLICIT NONE

      SELECT CASE (TRIM(title(1)))
          CASE ("run.07.300.so")
              find_assimilate_step_limit = 192
          CASE ("run.08.300.so")
              find_assimilate_step_limit = 48
          CASE ("run.08.400.so")
              find_assimilate_step_limit =  96
          CASE ("run.08.360.so")
              find_assimilate_step_limit = 192
          CASE ("run.09.360.so")
              find_assimilate_step_limit = 144
          CASE ("run.09.400.so")
              find_assimilate_step_limit = 48
          CASE ("run.10.300.so")
              find_assimilate_step_limit = 144
          CASE ("run.10.360.so")
              find_assimilate_step_limit = 96
          CASE ("run.10.400.so")
              find_assimilate_step_limit = 48
          CASE ("run.11.300.so")
              find_assimilate_step_limit = 96
          CASE ("run.11.360.so")
              find_assimilate_step_limit = 240
          CASE DEFAULT
              find_assimilate_step_limit = 48
        END SELECT
        RETURN

    END FUNCTION find_assimilate_step_limit


    SUBROUTINE do_eddy_hack(nt,X)

      USE kei_parameters
      USE kei_common
      USE kei_subs1D, ONLY: ABK80

      IMPLICIT NONE

       ! Input/Output
      REAL, INTENT (INOUT) :: X(NZP1,NSCLR)
      INTEGER, INTENT(IN) :: nt

      !IF(nt .eq. 3600 ) then !.or. nt .eq. 4200 .or. nt .eq. 5160) then
      ! write(6,*) 'Hacking in eddy heat'
        !X(20:22,1) = X(20:22,1) + 0.1
      ! X(1:27,1) = X(1:27,1) + 0.25
      ! X(28:100,1) = X(28:100,1) + 0.5
      !elseIF (nt .eq. 3696 .or. nt .eq. 4296 .or. nt .eq. 5254) then
      ! write(6,*) 'Hacking out eddy heat'
      ! !X(20:22,1) = X(20:22,1) - 0.05
      ! X(1:27,1) = X(1:27,1) - 0.15
      ! X(28:100,1) = X(28:100,1) - 0.5
      !ENDIF

    END SUBROUTINE do_eddy_hack


    subroutine assimilate_ts(X,kforce,nt)

      USE kei_parameters
      USE kei_common
      USE kei_subs1D, ONLY: ABK80

      IMPLICIT NONE

      ! Input/Output
      real, intent (INOUT) :: X(NZP1,NSCLR)
      type(kei_forcing_type), intent(IN) :: kforce   ! forcing data structure (from kei_common)
      INTEGER, INTENT(IN) :: nt

      ! Local Variables
      integer :: ii,ii_ft,ii_ud, ii_interp,ii_start, assim_24_counter
      real :: deltat,sigma_ft,sigma,sigma0,s_slope,intercept,alpha,beta,exppr

      ! rectify lower ocean, if enabled
      ii_ft = -1
      if (assimilate_ts_hack .and. (assim_counter < assimilate_step_limit)) THEN
        if (assimilate_ts_fixed  .or. kforce%f_wct == 1) THEN
          if (assimilate_ts_fixed) then
              if (.not. allocated(fixed_t)) THEN
                allocate(fixed_t(401),fixed_s(401))
                open(unit=21, file='fixed_ctd.txt', form='FORMATTED', access='sequential')
                do ii = 1,401
                  read(21,*) fixed_t(ii), fixed_s(ii)
                enddo
                close(21)
              endif

              ! overwrite every timestep
              if (nt==1) then
                ii_start = 60
              else
                ii_start = 290
              endif
              ii_ft = ii_start
              do ii=ii_start,NZP1
                 X(ii,1) = fixed_t(ii)
                !X(ii,2) = 0.1184*X(ii,1) + 34.56 - Sref - 0.25
                X(ii,2) = 0.1184*X(ii,1) + 34.56 - Sref - 0.13   ! added 0.12 from slope per Stammerjohn suggestion 09/2015
              enddo

          else

            if ((.not. assimilate_ts_stop) .or. (assim_counter < assimilate_step_limit)) THEN

                deltat = 1.0/24.0
                do ii=60,NZP1
        !          if (kforce%wct_interp(ii) .ne. -999. .and. kforce%wct_interp(ii) .gt. 0.) then
                  if (kforce%wct_interp(ii) .ne. -999.) then

                    ! record top of mooring forcing, increment assim_24_counter
                    if (ii_ft .eq. -1) then
                      ii_ft = ii
                      if (assimilate_ts_stop) THEN
                          assim_counter = assim_counter + 1
                      endif
                    endif

                    ! overwrite
                    X(ii,1) = kforce%wct_interp(ii)
                    !X(ii,2) = 0.06583*X(ii,1) + 34.56 - Sref ! t-s regression from 300-100 2007 mooring

                    ! nudge 24hr relaxation

                    !X(ii,1) = X(ii,1) + &
                    !  (kforce%wct_interp(ii) - X(ii,1)) * deltat

                    ! endpoints 33.9/-1.7, 34.35/2.1 of CDW/shelf water mixing line
                    ! slope 0.1184, intercept 34.1013

                    !X(ii,2) = 0.1184*X(ii,1) + 34.56 - Sref - 0.25
                    X(ii,2) = 0.1184*X(ii,1) + 34.56 - Sref - 0.13   ! added 0.12 from slope per Stammerjohn suggestion 09/2015

                  endif
                enddo

            endif

          endif

          if (ii_ft .gt. 0) then

            if (ii_ft .gt. last_wct_ii) then
              ! forcing is lower




              !ii_interp = last_wct_ii-1
              ii_interp = 50

              ! interp temperatures back across ii_interp to forcing depth
              !s_slope = (X(ii_ft,1)-X(ii_interp,1))/(zm(ii_ft)-zm(ii_interp))
              !intercept = X(ii_ft,1) - s_slope*zm(ii_ft)
              !do ii=last_wct_ii,ii_ft-1
              !  X(ii,1) = s_slope*zm(ii) + intercept
              !enddo

              ! interp salinities back across ii_interp to forcing depth

              s_slope = (X(ii_ft,2)-X(ii_interp,2))/(zm(ii_ft)-zm(ii_interp))
              intercept = X(ii_ft,2) - s_slope*zm(ii_ft)
              do ii=last_wct_ii,ii_ft-1
                X(ii,2) = s_slope*zm(ii) + intercept
              enddo

            endif


            ! find unstable densities above forcing
            alpha = 1.0
            beta = 1.0
            exppr = 0.0
            call ABK80(X(ii_ft,2)+Sref,X(ii_ft,1),-zm(ii_ft),alpha,beta,exppr,sigma_ft,sigma)
            ii_ud = -1
            do ii=1,ii_ft-1
              call ABK80(X(ii_ft,2)+Sref,X(ii_ft,1),-zm(ii),alpha,beta,exppr,sigma0,sigma)
              if(ii_ud .lt. 0 .and. sigma0 .gt. sigma_ft) then
                ii_ud = ii
              endif
            enddo

            ! find new stable salinities
            if (ii_ud .eq. 1) then
              ! hopefully entire water column is not unstable.
              ! set salinity to just below top forcing depth, i guess
              !X(1:ii_ft-1,2) = X(ii_ft,2) - 0.005
              ! above not anymore - mixed layer should mix anyway
              X(1:ii_ft-1,2) = X(ii_ft,2)
            elseif (ii_ud .gt. 1) then
              ! interpolate salinity between deepest stable salinity and top of forcing data
              s_slope = (X(ii_ft,2)-X(ii_ud-1,2))/(zm(ii_ft)-zm(ii_ud-1))
              intercept = X(ii_ft,2) - s_slope*zm(ii_ft)
              do ii=ii_ud,ii_ft-1
                X(ii,2) = s_slope*zm(ii) + intercept
              enddo

            endif

            ! update previous min depth (index) of assimilation data
            last_wct_ii = ii_ft

          endif

        else

          print *, 'WARNING - T/S assimilation requested but no T/S assimilation data present'

        endif

      endif

    end subroutine assimilate_ts



    SUBROUTINE find_ice_fe(start_year,ice_fe_yearly)

      IMPLICIT NONE

      INTEGER(KIND=int_kind), INTENT(IN) :: start_year
      REAL, INTENT(OUT) :: ice_fe_yearly

      IF (use_year_ice_advance) THEN

        ! actual ice advance day of year
        select case (start_year)
          case (1997)
            ice_fe_yearly = 121
          case (1998)
            ice_fe_yearly = 134
          case (1999)
            ice_fe_yearly = 165
          case (2000)
            ice_fe_yearly = 177
          case (2001)
            ice_fe_yearly = 172
          case (2002)
            ice_fe_yearly = 139
          case (2003)
            ice_fe_yearly = 162
          case (2004)
            ice_fe_yearly = 161
          case (2005)
            ice_fe_yearly = 125
          case (2006)
            ice_fe_yearly = 163
          case (2007)
            ice_fe_yearly = 162
          case (2008)
            ice_fe_yearly = 185
          case (2009)
            ice_fe_yearly = 185
          case (2010)
            ice_fe_yearly = 164
          case default
            ice_fe_yearly = 161  ! mean ice advance day 97-2010
        end select

        ! assume Apr 1 (91) is full fe inclusion
        ! assume July 1 (181) is no fe inclusion
        ice_fe_yearly = ice_fe * (181 - ice_fe_yearly)/(181-91)
        ice_fe_yearly = max(ice_fe_yearly,0.0)
        ice_fe_yearly = min(ice_fe_yearly,ice_fe)

      ELSE

        ice_fe_yearly = ice_fe

      ENDIF

    END SUBROUTINE find_ice_fe


    SUBROUTINE do_eco_hacks(X,nt,start_year)

      USE kei_parameters
      USE kei_common
      USE kei_ecocommon

      IMPLICIT NONE

      real, intent (INOUT) :: X(NZP1,NSCLR)
      integer(KIND=int_kind), INTENT(IN) :: nt, start_year

      integer(KIND=int_kind) :: fe_offset_local, bio_offset_i

      ! testing the addition of lots of alkalinity in top layers
      if (alk_spike_level_start > 0) then
        if (alk_spike_timestep == nt) then
          X(alk_spike_level_start:alk_spike_level_end,2+alk_ind) = alk_spike_value
        endif
      endif

      if (fe_offset < 0) then
        select case (start_year)
          case (1997)
            fe_offset_local = 70
          case (1998)
            fe_offset_local = 80
          case (1999)
            fe_offset_local = 65
          case (2000)
            fe_offset_local = 60
          case (2001)
            fe_offset_local = 50
          case (2002)
            fe_offset_local = 75
          case (2003)
            fe_offset_local = 60
          case (2004)
            fe_offset_local = 70
          case (2005)
            fe_offset_local = 85
          case (2006)
            fe_offset_local = 60
          case (2007)
            fe_offset_local = 65
          case (2008)
            fe_offset_local = 70
          case (2009)
            fe_offset_local = 100
          case (2010)
            fe_offset_local = 80
          case default
            fe_offset_local = 0
        end select
      elseif (fe_offset > 0) then
        fe_offset_local = fe_offset
      else
        fe_offset_local = 0
      endif

      if (nt .eq. 1) then
        X(:,2+Fe_ind) = X(:,2+Fe_ind)*fe_multiplier
      endif

      if (nt .eq. 1 .and. fe_offset_local > 0) then
        ! move Fe nutricline down by fe_offset
        X(6+fe_offset_local:401,2+Fe_ind) = X(6:401-fe_offset_local,2+Fe_ind)
        !X(6:6+fe_offset-1,2+Fe_ind) = X(5,2+Fe_ind)
        X(1:6+fe_offset_local-1,2+Fe_ind) = 3.0e-5

      endif
      if (nt .eq. 1 .and. bio_offset > 0) then
        bio_offset_i = int(bio_offset)

        ! move biology down by bio_offset_i
        X(1:bio_offset_i,2+spC_ind) = X(bio_offset_i+1,2+spC_ind)
        X(1:bio_offset_i,2+spChl_ind) = X(bio_offset_i+1,2+spChl_ind)
        X(1:bio_offset_i,2+spFe_ind) = X(bio_offset_i+1,2+spFe_ind)
        X(1:bio_offset_i,2+diatC_ind) = X(bio_offset_i+1,2+diatC_ind)
        X(1:bio_offset_i,2+diatChl_ind) = X(bio_offset_i+1,2+diatChl_ind)
        X(1:bio_offset_i,2+diatFe_ind) = X(bio_offset_i+1,2+diatFe_ind)
        X(1:bio_offset_i,2+diatSi_ind) = X(bio_offset_i+1,2+diatSi_ind)
        X(1:bio_offset_i,2+don_ind) = X(bio_offset_i+1,2+don_ind)
        X(1:bio_offset_i,2+doFe_ind) = X(bio_offset_i+1,2+doFe_ind)
        X(1:bio_offset_i,2+dop_ind) = X(bio_offset_i+1,2+dop_ind)
        X(1:bio_offset_i,2+zooC_ind) = X(bio_offset_i+1,2+zooC_ind)
      endif

      dust_flux = 0.0

    END SUBROUTINE do_eco_hacks


END MODULE kei_hacks