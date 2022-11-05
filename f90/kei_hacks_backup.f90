MODULE kei_hacks

  USE kei_ecocommon, ONLY: log_kind,int_kind,real_kind,dbl_kind

  IMPLICIT NONE
  
  ! overwrite in some eddy heat at depth during certain periods
  LOGICAL, PARAMETER :: eddy_hack = .false.
  LOGICAL, PARAMETER :: assimilate_ts_hack = .true.
  
  ! fraction precipitation adjustment OVER ICE ONLY
  REAL, PARAMETER :: snow_fraction = 0.70    
  REAL, PARAMETER :: rain_fraction = 0.0
  
  ! fraction shortwave irradiance
  REAL, PARAMETER :: shortwave_multiplier = 0.80

  ! conform to ice concentration found in forcing (0=no, 1=yes)
  INTEGER, PARAMETER :: ic_conform = 0

!-----------------------------------------------------------------------
!     initial fe/bio profile adjustments
!-----------------------------------------------------------------------
  
  real (kind=dbl_kind), parameter :: &
    fe_multiplier = 1.0, &
    fe_offset = -100., &
    bio_offset = 100., &
    meltwater_fe = 0.0
!      meltwater_fe = 12.0e-9

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
  
  
    subroutine assimilate_ts(X,kforce)
      
      USE kei_parameters
      USE kei_common
      USE kei_subs1D, ONLY: ABK80

      IMPLICIT NONE
     
      ! Input/Output
      real, intent (INOUT) :: X(NZP1,NSCLR)
      type(kei_forcing_type), intent(IN) :: kforce   ! forcing data structure (from kei_common)
      
      ! Local Variables
      integer :: ii,ii_ft,ii_ud, ii_interp
      real :: deltat,sigma_ft,sigma,sigma0,s_slope,intercept,alpha,beta,exppr

      ! rectify lower ocean, if enabled
      if (assimilate_ts_hack) THEN
 
        deltat = 1.0/24.0
        if (kforce%f_wct == 1) then
          ii_ft = -1
          do ii=60,NZP1
  !          if (kforce%wct_interp(ii) .ne. -999. .and. kforce%wct_interp(ii) .gt. 0.) then
            if (kforce%wct_interp(ii) .ne. -999.) then

            ! record top of mooring forcing
              if (ii_ft .eq. -1) then
                ii_ft = ii
              endif
            
              ! overwrite
              X(ii,1) = kforce%wct_interp(ii)
              !X(ii,2) = 0.06583*X(ii,1) + 34.56 - Sref ! t-s regression from 300-100 2007 mooring 

              ! nudge 24hr relaxation
            
              !X(ii,1) = X(ii,1) + &
              !  (kforce%wct_interp(ii) - X(ii,1)) * deltat    
              
              ! endpoints 33.9/-1.7, 34.35/2.1 of CDW/shelf water mixing line
              ! slope 0.1184, intercept 34.1013
              X(ii,2) = 0.1184*X(ii,1) + 34.56 - Sref - 0.25
            endif
          enddo
  
      
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
              X(1:ii_ft-1,2) = X(ii_ft,2) - 0.005
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
  
    SUBROUTINE do_eco_hacks(X,nt,start_year)

      USE kei_parameters
      USE kei_common
      USE kei_ecocommon
    
      IMPLICIT NONE
    
      real, intent (INOUT) :: X(NZP1,NSCLR)
      integer(KIND=int_kind), INTENT(IN) :: nt, start_year
      
      integer(KIND=int_kind) :: fe_offset_local, bio_offset_i
        
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