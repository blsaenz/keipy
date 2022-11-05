
    SUBROUTINE atm_step(kforce)

      use kei_parameters
      use kei_common

      implicit none

      ! Input
      real :: kforce(forcing_var_cnt)
      DOUBLE PRECISION :: &
        day_of_year,ros,sza,rod,romean,ws,d_lat,d_lon

      ! determine shortwave surface reflectance
      day_of_year = time
      d_lat = dble(dlat)
      d_lon = dble(dlon)
      do while(day_of_year .gt. 365.)
        day_of_year = day_of_year - 365.
      enddo

      call zenith(day_of_year,d_lat,d_lon,sza)
      ws = sqrt(kforce(taux_f_ind)**2 + &
                kforce(tauy_f_ind)**2)
      call atm_sfcrfl(ws,sza,ros,rod)

      romean = (ros + rod)/2.0

      ! load into vaf variable used by atm routines
      vaf(1) = kforce(taux_f_ind)
      vaf(2) = kforce(tauy_f_ind)
      !vaf(3) = kforce%f_interp(qswins_f_ind)*(1.0-romean)*sw_scale_factor
      vaf(3) = kforce(qswins_f_ind)*sw_scale_factor
      vaf(4) = kforce(qlwdwn_f_ind)
      vaf(5) = kforce(tz_f_ind)
      vaf(6) = kforce(qz_f_ind)
      vaf(7) = kforce(prain_f_ind)
      vaf(8) = kforce(psnow_f_ind)
      vaf(9) = kforce(divu_f_ind)
      vaf(11) = kforce(msl_f_ind)
      vaf(12) = kforce(h_f_ind)

    end SUBROUTINE atm_step


    ! not calling -
!     SUBROUTINE init_atm(timed,X,atmtemp,atmshum,kforce)
!
! !     initialize the atmosphere
!
! !     "atm.f" is used instead of "atmrad.f" when there is no rad/conv
! !     code included.
!
!     use kei_parameters
!     use kei_common
!     use kei_icecommon
!     use kei_unit
!
!     implicit none
!
! ! input
!     DOUBLE PRECISION :: timed ! time (in days), at which initial forcing
! !                              array VAF gets computed
!       real :: &
!         X(NZP1,NSCLR)! needed for salinity correction
! ! output
!     real :: atmtemp(nplev),atmshum(nplev)
!     type(kei_forcing_type) :: kforce   ! forcing data structure
!
! !     local variables
!     integer :: j,i
!     real :: umag, freshwater, weight, mean_sal
!
!     call kei_init_forcing(kforce)
!
!     ! find salinity addition from total freshwater input, hourly time step:
!     freshwater = 3600 * &
!         (sum(kforce%f_data(:,prain_f_ind)) &
!          + sum(kforce%f_data(:,psnow_f_ind)))
!     ! find mean initial salinity
!     mean_sal = 0.
!     weight = 0.
!     do i=1,nzp1
!         mean_sal = mean_sal + (X(i,2)+Sref) * hm(i)
!         weight = weight + hm(i)
!     enddo
!     mean_sal = mean_sal/weight
!     ! total salt deficit (kg) ~= fresh_mass (kg) * mean ptt (kg/kg/1000)
!     ! salinity correction rate = total salt deficit (kg) / total days (day) / (s/day)
!     sal_correction_rate = &
!         freshwater * mean_sal / 1000 &
!         / maxval(kforce%f_data(:,date_f_ind)) / 86400.  ! (kg/s)
!
!     ! too high b/c of evap, potential rain problem in ice model
!     sal_correction_rate = sal_correction_rate / 2.
!
!     write(6,*) 'sal_correction_rate: ', sal_correction_rate
!
! !    sal_correction_rate = 0.0
!
! ! Initialize variables for repeat data reads: set offset for reads
!
!     if (lrepeatdat) then
!         offdatday   = 0.
!     ! restart
!     !         offdatday   = 18250. + 2996.
!         startdatday = startt
!     endif
!
! ! Initialize noise for harmonics
!
!     IF( (LAFLX == -1) .OR. (LAFLX == +1) .OR. (LAFLX == +2) ) THEN
!         do 40 j=1,NSFLXSM1
!             xnoise(j) = 0.0
!         40 END DO
!     ENDIF
!
! ! Initialize atm parameters (VAF) at startt
!
!     call preatm(timed,kforce)
!
! ! If RAD/CONV model is used:
! !     set atmospheric variables that are needed for calculating
! !     the fluxes in atmflx when called by init flx.
!
! !     IF((LAFLX.eq.2).or.(LAFLX.eq.3).or.(LAFLX.eq.5)) THEN
! ! tm         tz = t(1,plev) - tk0
! ! tm         qz = 0.77 * qsat(0, t(1,plev) )
! !     ENDIF
!
! ! If forcing via state variables:
! !     fudge fluxes  needed for initial "ldflx"
!
!     IF(LAFLX >= 1) THEN
!         umag = sqrt(vaf(1)**2 + vaf(2)**2)
!         sflux(1,1,0) = 1.2*umag*vaf(1) * (0.0012 * focn + 0.0030 * fice)
!         sflux(2,1,0) = 1.2*umag*vaf(2) * (0.0012 * focn + 0.0030 * fice)
!         sflux(5,1,0) = 0.0
!
!     ENDIF
!
!     return
!     end SUBROUTINE init_atm

!***********************************************************************

!     SUBROUTINE preatm(timed,kforce)
!
! !     prepare for an atm timestep
! !     load the atm flux parameter array VAF(nsflxs)
! ! NB   set up Equivalences in COMMON according to LAFLX
!
! !      VAF (1)  (2)  (3)  (4) (5) (6)  (7) (8) (9) (10) (11) (12)
! !      -4 tauU tauV QSWo QLWn  Hs  E   Prn Psn     QSWi QSWl  :data file
! !      -1 tauU tauV QSWo QLWn  Hs  E   Prn Psn     QSWi QSWl  :harmonics
! ! AFLX= 0 tauU tauV QSWo QLWn  Hs  E   Prn Psn     QSWi QSWl  :constant
! !       1  uZ   vZ  QSWo QLWn  tZ qZ   Prn Psn     QSWi QSWl  :harmonics
! !       2  uZ   vZ                     Prn Psn                :harmonics
! !                   QSWo QLWn  tZ qZ               QSWi QSWl  :atm model
! !       3  uZ   vZ  QSWo QLWn  tZ qZ   Prn Psn     QSWi QSWl  :atm model
! !       4  uZ   vZ Cloudf SST  tZ Tdew Prn Psn     ***  ***   :data file
!
! !       5  uZ   vZ                     Prn Psn                :data file
! !                   QSWo QLWn  tZ qZ               QSWi QLWl  :atm model
! !       6  uZ   vZ  SWdn LWdn  tZ qZ   Prn Psn DivU SST SWup LWdn
!
!     use kei_parameters
!     use kei_common
!     use kei_unit
!
!     implicit none
!
! ! Input
!     DOUBLE PRECISION :: timed        ! time for getting forcing
!     type(kei_forcing_type) :: kforce   ! forcing data structure
! ! Local
!     DOUBLE PRECISION :: datday
!     real :: dayfrac
!
! ! Diagnostic
!
! ! Constant flux forcing ?
!
!     if(LAFLX == 0) return
!
! ! Get time variable forcing
!
! ! If forcing specified by harmonics, calculate VAF
!
!     if ( iabs(laflx) == 1 .OR. laflx == 2)  call FCOMP(timed)
!
! ! If forcing specified in datafile, read VAF
!
!     if ( iabs(laflx) == 4 .OR. laflx == 5 .OR. laflx == 6 ) then
!
!     !        Repeat data input
!
!         if (lrepeatdat) then
!             datday = timed - offdatday
!         !            write(6,*) 'forcing data: timed=',timed,' datday=',datday,
!         !     $                 ' offdatday=',offdatday
!             if (datday > (startdatday+repeatdatdays)) then
!             !              adjust offdatday, but keeping same fraction of day as
!             !              model time "timed", set new startdatday,
!             !              and start reading from beginning of data file again
!                 11 write(6,*) 'repeat data: model time = timed =',timed
!                 write(6,*) '       old datday=',datday,' offdatday=', &
!                 offdatday,' startdatday=',startdatday
!                 dayfrac   = datday - int(datday)
!                 datday    = datday - repeatdatdays
!                 datday    = int(datday) + dayfrac
!                 offdatday = timed - datday
!             !              check necessary for restarts
!                 if(datday > (startdatday+repeatdatdays)) goto 11
!                 rewind(nuforc)
!                 !call initfread(datday)
!                 startdatday = fdata(1,idataold)
!                 write(6,*) '       new datday=',datday,' offdatday=', &
!                 offdatday,' startdatday=',startdatday
!             endif
!             call FREAD(datday,kforce)
!         else
!             call FREAD(timed,kforce)
!         endif
!     endif
!
!     return
!     end SUBROUTINE preatm
!
! !***********************************************************************
!
!     SUBROUTINE fcomp(tday)
!
! !     Compute forcing with harmonics. Note: parameter on function
! !     gausd with 0 argument returns next random number in sequence
! !           with 1          restarts the geneartor
! !           with any other uses idum as seed for generator.
! !     If harm0(3) for QSW is specified by a negative number, solar
! !     radiation comes from geometry and simple cloud model.
! !     If mean stresses harm0(1) and harm0(2) are 999.9 then
! !              compute analytic winds.
!
!     use kei_parameters
!     use kei_common
!
!     implicit none
!
!     interface
!      FUNCTION gausd(idum)
!      implicit none
!       integer, intent(in) :: idum
!       integer :: iset = 0
!       real :: v1,v2,r,fac,gset
!       real :: gausd
!      end function gausd
!     end interface
!
!     ! input
!     DOUBLE PRECISION :: tday
!
!     ! local
!     integer :: j,i
!     real :: at, dt, pifac,Tidays
!
!     do 100 j=1,nsflxsm1
!     !     xnoise(j) = areg(j) * xnoise(j) + ampharm(0,j) * gausd(0)
!         VAF(j) = harm0(j) + xnoise(j)
!         if(nharm(j) >= 1) then
!             do 110 i=1,nharm(j)
!                 VAF(j) = VAF(j)+ampharm(i,j)*cos(freharm(i,j)*tday-phharm(i,j))
!             !      write(6,*) 'flux: vaf=',vaf(j),' fre=',freharm(i,j),' tday=',
!             !    +           tday,' phharm=',phharm(i,j)
!             110 END DO
!         endif
!     100 END DO
!
! ! SW solar forcing from geometry and simple cloud model
!     if(harm0(3) < 0.) then
!     !         jokta = 4              ! = cloudiness (oktas)
!         call solar( -int( harm0(3) ),tday)
!     endif
!
! ! *********************  analytic winds for inertial runs **********
!     IF((harm0(1) > 999.) .AND. (harm0(2) > 999.)) THEN
!         if(tday <= (16./24.)) then
!             at    = 0.6 * sin( 3.1417 * tday *24. / 16. )
!             pifac =  2. * 3.1417 * 24. / 16.
!             vaf(1) =-at * sin( pifac * tday )
!             vaf(2) =-at * cos( pifac * tday )
!         else
!             vaf(1) = 0.0
!             vaf(2) = 0.0
!         endif
!     ENDIF
!     IF((harm0(1) < -999.) .AND. (harm0(2) < -999.)) THEN
!         Tidays  = 314.17 /18./24.
!         pifac =  2. * 3.1417 / Tidays
!         vaf(1) = 0.2 * sin( pifac * tday )
!         vaf(2) = 0.2 * cos( pifac * tday )
!     ENDIF
! !      write(6,*) 'exit analytic',tday,vaf(1),vaf(2)
!
!     return
!     end SUBROUTINE fcomp
!
! !***********************************************************************



!***********************************************************************

!     REAL FUNCTION gausd(idum)
! !     for changing random number with uniform deviates (0-1)
! !     to a normal or gaussian distribution using the Box-Muller Method
! !     see Numerical Recipes p202-203
! !     (Note: for double precision "rand" needs to be changed to"drand".)
!
! #ifdef ifort
!     use IFPORT
! #endif
!
!     implicit none
!
!     integer, intent(in) :: idum
!     integer :: iset = 0
!     real :: v1,v2,r,fac,gset
!
!     if (iset == 0) then
!         10 v1=2.*rand(idum)-1.
!         v2=2.*rand(idum)-1.
!         r=v1**2+v2**2
!         if(r >= 1) goto 10
!         fac = sqrt(-2.*log(r)/r)
!         gset=v1*fac
!         gausd =v2*fac
!         iset=1
!     else
!         gausd =gset
!         iset=0
!     endif
!     return
!     end FUNCTION gausd

!************************************************************

!     subroutine atmstep(timed  ,rlwup  ,qsen   ,qlat  ,tg    , &
!     focn   ,fice   ,flnd   ,snhice,ti0   ,tmlt, &
!     laflx  ,icount ,rhdat  ,clddat, &
!     lrhdat ,lclddat, &
!     modatm ,atmad , &
!     rswdno ,rswdni ,rswdnl ,rlwnet,tz    ,rhz , &
!     atmt   ,atmh2ommr,atmsave)
!
! !     dummy subroutine, from rad/conv code
!
!     return
!     end subroutine atmstep




!************************************************************

!     subroutine prtprof
!
! !     dummy subroutine, from rad/conv code
!
!     return
!
!
!     end subroutine prtprof
!
! !***************************************************************
!     subroutine getdat(nstart , &
!     gravx  ,cpairx ,epsilox,stebolx,   nrow, &
!     clat   ,oro    ,sndpth ,ts     ,     tg, &
!     ps     ,pmid   ,pint   ,pmln   ,   piln, &
!     t      ,h2ommr ,o3mmr  ,plol   ,plos   , &
!     cldfrc ,clwp   ,effcld , &
!     label  ,tlabel ,dlabel)
! !---------------------------------------------------------------
!
!     return
!     end subroutine getdat
!
! !*********************************************************************
!
!     subroutine radini(gravx,cpairx,epsilox,stebolx)
! !---------------------------------------------------------------------
!
! ! modified version of initialization for radiation scheme; done
! ! for the column radiation model
! !c sets constants in crdcon and crdcae
! !c incudes former stparm-subroutine
!
!     return
!     end subroutine radini

!***********************************************************************

!     SUBROUTINE solar(jokta,timed)
!
! !     compute solar shortwave flux (W m-2)
! !     lat and long are in deg., time is in days
! !     see Stull (1988)
!
! !     change(3-17-93): include eccentricity of the earth's orbit,
! !                      (r/r0)**2
!
!
!     use kei_parameters
!     use kei_common
!
!     implicit none
!
!     ! input
!     integer :: jokta
!     DOUBLE PRECISION :: timed
!
!     ! local
!     real :: s0,day,hour,ror0,ds,sinphi,tr
!
!     s0 = 1368.0
! !     deg. to radians
! !         ztime = startt + (ntime + ndtatm/2.)*dt/spd
!     day = aint(timed)
!     hour = timed-day
!
! !     compute sine of solar elevation (sinphi) as a function of time
!
!     ror0   = 1. + 0.033*cos(twopi*day/365.25)
!     ds = 0.409 * cos(twopi*(day-173)/365.25)
!     sinphi = sin(rlat)*sin(ds) - cos(rlat)*cos(ds)* &
!     cos(twopi*hour+ rlon)
!
! !     calculate transmitance (tr) and surface short wave
!
!     if(sinphi >= 0) then
!     !           call cloud(jokta,sinphi,tr)
!         call cloudpapa(jokta,sinphi,tr)  ! cloud correction for papa
!         vaf(3) = s0 * ror0 * sinphi * tr
!     else
!         vaf(3) = 0.
!     endif
!
!     RETURN
!     end SUBROUTINE solar
!
! !**********************************************************************
!
!     SUBROUTINE cloudpapa(JCLOU,S,CFAC)
!
! !     Calculate cloud transmittance using method of Dobson&Smith (1988).
! !     These are regression fits for weather station PAPA (1959-1975),
! !     Modified: clear sky transmittance is too low, so jclou >= 1.
!
! !     Input variables:
! !     jclou = cloudiness (oktas)
! !     s     = sine of the solar elavation
! !     Output variable:
! !     cfac  = cloud transmittance
!
!     implicit none
!
!     ! inputs
!     integer :: JCLOU
!     real :: S,CFAC
!
! !     jokta   =  0    1    2    3    4    5    6    7    8    9
!     integer :: jokta         ! = max(1,jclou)
!     real :: a(0:9) = (/ .400,.517,.474,.421,.380,.350,.304,.230,.106, &
!       .134 /)    ! okta 9 = sky obscured
!     real :: b(0:9) = (/ .386,.317,.381,.413,.468,.457,.438,.384,.285, &
!       .295 /)
!
!     jokta = max(1,jclou)
!
!     CFAC  = a(jokta) + b(jokta) * S
!
!     RETURN
!     end SUBROUTINE cloudpapa
!
! !**********************************************************************
!
!     SUBROUTINE cloud(JCLOU,S,CFAC)
!
! !     Calculate cloud transmittance using method of Smith & Dobson
! !     (1984).
!
!       implicit none
!
! !     Input variables:
!       integer :: jclou ! cloudiness (oktas)
!       real :: s     ! sine of the solar elavation
! !     Output variable:
!       real :: cfac  ! cloud transmittance
!
!       ! local
!       real :: C
!       REAL, DIMENSION(0:8) :: A,B,D,E
!
!       A(6:8) = (/0.310,0.235,0.103/)
!       B(6:8) = (/0.439,0.388,0.296/)
!       D(0:5) = (/0.240,0.070,-0.010,0.055,0.070,0.090/)
!       E(0:5) = (/0.0520,0.0525,0.0430,0.0395,0.0375,0.0345/)
!
!       if (JCLOU <= 5) then
!           C=JCLOU/8.
!           CFAC=(E(JCLOU)+S*EXP(-D(0)/S)*(C*EXP(-D(JCLOU)/S)+1.-C))
!       else
!           CFAC=S*(A(JCLOU)+B(JCLOU)*S)
!       endif
!
!       RETURN
!     end SUBROUTINE cloud
!
! !**********************************************************************


   subroutine zenith(day_of_year,lat,lon,zenith_angle)

      real(8), intent(in) :: &
        day_of_year, &
        lat, lon
      real(8), intent(out) :: &
        zenith_angle

      real(8) :: &
        solar_dec,rad_date,rad_time

      real(8), parameter :: &
        pi = 3.1415926535898D0, &
        c2 = 2.d0, &
        pi2 = pi*c2


      ! convert julian date to radians
      rad_date = (day_of_year/365)*pi2

      ! Calculate solar declination (between -23deg27' and 23deg27') for use in
      ! Solar elevation equation
      solar_dec = 0.39637 - 22.9133*cos(rad_date) + 4.02543*sin(rad_date)  &
          - 0.3872*cos(2.0*rad_date) + 0.052*sin(2.0*rad_date)

      ! Convert degrees to radians
      solar_dec = pi*solar_dec/180.

      ! Convert time to radians - these equations seem to give midnight as noon,
      ! so I added pi to the time to correct -
      rad_time = (day_of_year - int(day_of_year))*c2*pi + pi

      zenith_angle = acos(sin(solar_dec)*sin(lat/360*pi2) &
      + cos(solar_dec)*cos(lat/360*pi2)*cos(rad_time) )   ! radians

      ! Convert radians to degrees
      zenith_angle = 180.*zenith_angle/pi

  end subroutine zenith


!*******************************************************************

       subroutine atm_sfcrfl(ws,sza,ros,rod)

!---------------------------------------------------------------------
! declare local variables
!---------------------------------------------------------------------
!     cdrag: drag coefficient
!    drspec: direct specular reflectance
!    dfspec: diffuse specular reflectance
!      foam: sea-foam reflectance
!   dd1,dd2
!   dd3,dd4: coefficients relating wind stress to foam reflectance
!
!      sinr: sine of refracted angle
!    refrac: refracted angle of direct specular reflectance
!
!     lrho : rhoa in units of g/m3

      real(8), intent(in) :: ws, sza !(degrees)
      real(8), intent(out) :: ros,rod


      real(8), parameter :: &
        radian = 180.0D0/3.1415926535898D0


      real(8) ::nlrho
      real(8) ::bcdrag, drspec, dfspec, foam
      real(8) ::bdd1, dd2, dd3, dd4
      real(8) ::btheta2, sinr, refrac, diff, sum, sin_term, tan_term
      real(8) ::ba, b

      dd1 = 2.2d-5
      dd2 = 4.0d-4
      dd3 = 4.5d-5
      dd4 = 4.0d-5

      rhoa = 1.270
      rn = 1.341d0

!---------------------------------------------------------------------
! Foam and diffuse specular reflectance
!---------------------------------------------------------------------

      lrho = rhoa*1000.  ! rhoa: reference air density (kg/m3)
      if (ws .gt. 4.0)then  ! ws: wind speed (m/s)
       if (ws .le. 7.0)then
           cdrag = (0.62d0 + 1.56d0/ws)*1.0d-3
           foam = dd1*lrho*cdrag*ws*ws - dd2  ! dd1, dd2: coefficients relating wind stress to foam reflectance
       else
           cdrag = (0.49d0 + 0.065d0*ws)*1.0e-3
           foam = (dd3*lrho*cdrag - dd4)*ws*ws! dd3, dd4: coefficients relating wind stress to foam reflectance
       endif
       dfspec = 0.057
      else
         foam = 0.0d0
       dfspec = 0.066
      endif

!---------------------------------------------------------------------
! Direct reflectance. For theta < 50 and ws < 2 m/s calculate using
! Fresnel's law, otherwise use an empirical fit (which is different
! than that shown in paper)
!---------------------------------------------------------------------

      if (sza .lt. 50.0 .or. ws .lt. 2.0)then
       if (sza .eq. 0.0)then
        drspec = 0.0211
       else
        theta2 = sza/radian
          sinr = sin(theta2)/rn
        refrac = asin(sinr)
        diff = theta2 - refrac
         sum = theta2 + refrac
         sin_term = (sin(diff)*sin(diff))/(sin(sum)*sin(sum))
         tan_term = (tan(diff)*tan(diff))/(tan(sum)*sin(sum))
         drspec = 0.5*(sin_term + tan_term)
       endif
      elseif (sza .lt. 90.0) then
       a =  5.25E-4*ws + 0.065
       b = -1.67E-3*ws + 0.074
       drspec = a*exp(b*(sza-60.0))

       if (drspec .lt. 0.) then
        drspec = dfspec
       endif
      else
       drspec = dfspec  ! no direct light anyhow
      endif

!---------------------------------------------------------------------
! Surface specular reflectance totals
!---------------------------------------------------------------------

       ros = dfspec + foam
       rod = drspec + foam

      return
  end

