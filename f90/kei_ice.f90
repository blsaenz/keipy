module kei_ice

    use kei_parameters, only: NZP1,NSCLR
    use kei_icecommon
    use kei_ecocommon, only: ice_to_ocean_eflux,ecosys_tracer_cnt, &
      Fe_ind,diatC_ind,diatFe_ind,diatChl_ind,diatSi_ind
    use sia2_types
    use kei_hacks, ONLY: snow_fraction, rain_fraction, ic_conform, &
      ice_diatChl, ice_fe, find_ice_fe

    implicit none

    type(ice_type), save, private :: &
      ice           ! where all ice data is held

    ! COMMON statements from ice.f
    ! ------------------------------------------------------------------
    ! Tarray
    real, save, private  :: Tair

    ! Farray
    real, save, private  :: Flws,Fsws,Flhs,Fml,Fis,Fm,Fshs,Flwo,dhdt

    ! harray
    integer, save, private  :: ni,ns
    real, save, private  :: hi,hs,hfs,Aic

    type(ice_pack_type), save :: ice_pack

contains

!***********************************************************************
! subroutine: init_ice
! Initializes ice at start of run (if any ice is specified), sets up
! common variables, reports back the heat, fresh water, and salt content
! of the ice
! ----------------------------------------------------------------------
    SUBROUTINE init_ice (t_ocn,s_ocn,hc_ice,fc_ice,sc_ice,start_year)

    use kei_icecommon
    use sia2_constants
    use sia2_parameters
    use sia2_grid
    use sia2_types

    implicit none

    integer :: i,start_year
    real :: tmp1,tmp2,tmp3,t_ocn,s_ocn
    double precision :: hc_ice,fc_ice,sc_ice


    ! initialize hc_ice (heat content) and sc_ice (salt content)
    Qsfc = 0.  ! net solar surface flux

    ! assign icecommon block names to local variable names
    hs     = hsn(0)
    hi     = hice(0)
    hfs    = hfsnow(0)
    Aic    = fice
    Tas    = TI0  + kelvin0
    dhdt   = DHDTice

    ! assign melt factor multiplier/fudge
    select case (start_year)
      case (1997)
        import_melt_f = 1.0
      case (1998)
        import_melt_f = 1.0
      case (1999)
        import_melt_f = 1.0
      case (2000)
        import_melt_f = 1.0
      case (2001)
        import_melt_f = 1.0
      case (2002)
        import_melt_f = 1.0
      case (2003)
        import_melt_f = 1.0
      case (2004)
        import_melt_f = 1.0
      case (2005)
        import_melt_f = 1.0
      case (2006)
        import_melt_f = 1.0
      case (2007)
        import_melt_f = 1.0
      case default
        import_melt_f = 1.0
    end select

    CALL find_ice_fe(start_year,ice_fe_local)

    call sia2_null_ice(ice)

    if (hi .gt. z_th_min) then

      ! create ice with 1st year ice profile
      call sia2_create_ice(ice,1,hi,hs,fice,Tair,t_ocn,s_ocn)

      ! find initial ice cotents
      call ice_hfs(ice,tmp1,tmp2,tmp3)
      hc_ice = tmp1
      fc_ice = tmp2
      sc_ice = tmp3

      print *, ' '
      print *, 'Initial call to ice_hfs'
      print *, 'hc [J/m2]      = ', hc_ice
      print *, 'fc_ice [kg/m2] = ', fc_ice
      print *, 'sc_ice [kg/m2] = ', sc_ice
      print *, ' '

    else

      ! nuke ice/snow
      hc_ice = 0.
      fc_ice = 0.
      sc_ice = 0.

    endif

    ! update common ice variables from sia2 ice structure, find albedo, shortwave stuff
    call ice_report(ice,0)

    print *, ' '
    print *, 'Initial T Profile:'
    write (6,700) Tas-kelvin0
    do i=1,ns+ni
        if (i <= ns) then
            write (6,701) i, Ts(i)-kelvin0
        else
            write (6,702) i-ns, Ti(i-ns)-kelvin0
        end if
    enddo
    print *, ' '

    700 format (/,'TI0',f9.3)
    701 format ('snow layer',i5,f9.3)
    702 format ('ice layer ',i5,f9.3)

  end SUBROUTINE init_ice


!***********************************************************************

  SUBROUTINE icestep (nt,nisteps,timed,tmix,smix,rhosw,Ta,Qz,shum,v10, &
    msl,DivU,ic,ain,aout,sflux,NSFLXS_local,hc_ice,fc_ice,sc_ice,flx,X)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_types
    use sia2_flux_heat
    use sia2_parameters
    use sia2_state

    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    integer, intent(in) :: &
      nt,             & !
      nisteps
    double precision, intent (in) :: &
      timed
    real, intent(in) :: &
      tmix,           & !
      smix,           & !
      rhosw,          & !
      Ta,             & !
      Qz,             & !
      shum,           & !
      v10,            & !
      msl,            & !
      DivU,           & !
      ic,             & !
      ain,            & !
      aout
    real, dimension(NSFLXS_local,5), intent(inout) :: &
      sflux
    integer, intent(in) :: &
      NSFLXS_local      !
    double precision, intent (inout) :: &
      hc_ice,         & !
      fc_ice,         & !
      sc_ice            !
    real :: &
      flx(11)   ! diagnostic flux data structure
    real, intent(in) :: &
      X(NZP1,NSCLR)


  ! Internal Variables
  ! --------------------------------------------------------------------
    integer :: &
      n_sia2_steps,   & !
      jjj               !
    real :: &
      s_ocn,          & !
      t_ocn,          & !
      sps,            & !
      rps,            & !
      rhow,           & !
      hc_old,         & !
      hc_new,         & !
      temp,           & !
      temp2,          & !
      temp3,          & !
      div,            & !
      ain1,           & !
      aout1,          & !
      ic1,            & !
      regrid_trigger, & !
      dtt_s,          & !
      snow_gl,        & !
      s_gl,           & !
      c_gl              !

    type(forcing_type) :: ff          ! forcing structure required by sia2 routines
    real, dimension (n_dh) :: dh  ! array for tracking changes in ice/snow
    type (heat_pre_calc) :: &
      pc                      ! pre-calculated surface heat flux parameters


  ! Function Code
  ! --------------------------------------------------------------------

    ! make local copies of external parameters and fluxes
    Tair = Ta+kelvin0   ! ice model expects Tair in Kelvin
    rhow = rhosw
    t_ocn = tmix
    s_ocn = smix
    Flws = sflux(4,3)                ! incoming longwave
    Fsws = (1.-albice)*sflux(3,3)    ! incoming shortwave
    Flhs = sflux(6,3)*SL             ! latent
    Fshs = sflux(5,3)                ! sensible
    Flwo = sflux(9,3)                ! outgoing longwave
    dhdt   = DHDTice                 ! set sensible/latent turbulent flux term -  added Saenz 10/2011
    Fml  = sflux(8,4)        ! melt potential [W/m2], <=0
    Fm   = -sflux(7,4)*Fl    ! ocean frazil flux, [kg/m2/s]*[J/kg]=[W/m2], >=0

    ! precipitation partitioning
    sps  = sflux(8,3)*snow_fraction        ! snow precip [kg/m2/s]
    rps  = sflux(7,3)*rain_fraction        ! rain precip [kg/m2/s]
    IF (rain_fraction > c0) THEN
      if (Ta < 0.) then
          sps = sps + rps
          rps = 0.;
      elseif (Ta > 2.) then
          rps = rps + sps
          sps = 0.
      endif
    ENDIF

    div = DivU        ! SSM/I convergence/magical divergence or melt

    ain1 = ain
    aout1 = aout
    ic1 = ic

    !if (fml .lt. 0.) then
    !  div = 0.01/86400. ! = DivU
    !else
    ! div = 0.
    !endif

    ! initalizations
    hc_old = 0.               ! pre main loop heat content
    hc_new = 0.               ! post main loop heat content
    hc_ice = 0.               ! post main loop heat content
    fc_ice = 0.               ! pre main loop salt content
    sc_ice = 0.               ! post main loop salt content
    lateral_freshwater_melt_flux = 0.  ! output variable
    atm_flux_to_ice_surface = 0.     ! output variable
    ice_ocean_bottom_flux_potential = 0.
    ice_ocean_bottom_flux = 0.
    total_ice_melt  = 0.  ! Joules
    total_ice_freeze  = 0. ! Joules
    frazil_ice_volume = 0.
    congelation_ice_volume = 0.
    snow_ice_volume = 0.
    snow_precip_mass = 0.
    ice%flux = 0.
    dh = 0. ! vector assignment
    regrid_trigger = z_th_min / 3.0
    ff%at = Tair ! K
    ff%t = t_ocn ! degC
    ff%s = s_ocn
    ff%d = rhow*1000. ! g/m^3
    ff%no3 = 1.
    ff%nh4 = 1.
    ff%po4 = 1.
    ff%sioh4 = 1.
    ff%poc = 1.

    ice_to_ocean_eflux = 0.0D0 ! vector assignment

    flx(6:11) = 0.

    if(nt .eq. 2961 .or. nt .eq. 2962 ) then
        write(6,*) 'zero ice break'
    endif


    print *,' '
    print *,'Flws, Fsws, Flhs, Fshs, Flwo = ', &
    Flws, Fsws, Flhs, Fshs, Flwo
    print *, 'Fml (8,4), Fm (7,4), sps = ', Fml, Fm, sps
    print *, ' '
    print *, 'hi, hs, hfs   = ', hi, hs, hfs
    print *, 'ni, ns        = ', ni, ns
    print *, 'dzi = ', dzi
    print *, 'dzs = ', dzs
    print *, 'nt = ', nt
    print *, 'Tair = ', Tair-kelvin0,Tair
    print *, 'Ts = ', Tas-kelvin0,Tas

    ! set timing (specified in run.xxxx)
    n_sia2_steps = dti/720. ! 12 minute dt_ice steps - this is messed up if dti not divisible by 720...
    dtt_s = 720.

    ! calculate starting Heat Content (to be used in iceflx/sflux(4,4) calc)
    if (hi > 0) then
        print *,' '
        print *,'second call to calc_hc (pre main loop)'
        call sia2_hc (ice,temp,hc_old)
        print *, 'hc_old = ', hc_old
    end if


    ! CATEGORIES:::::::::::::::::::::::::::::::::::::::::::::::::::
    ! Dice, new ice call outside of loop. lateral melt inside loop? maybe not, apply to all categories, since it doesn't matter?

    ! CATEGORIES:::::::::::::::::::::::::::::::::::::::::::::::::::
    ! iterate over ice categories

    call sia2_heat_pre_calc(Tair,v10,msl,shum,-1.,-1., &
      ice%snow%depth,-1.,Fsws,Flws,pc)

    do jjj=1,n_sia2_steps

      print *, ' '
      print *,'sub-ice step in   ',jjj

      if ((ic_conform .eq. 1 .and. (ic .gt. c0 .or. ice%af .gt. c0)) .or. &
        ic_conform .ne. 1) then

        ! call timestep
        call ice_sub_step(ice,ff,div,ic1,ain1,aout1,sps,rps,Tair,dtt_s, &
          ff%t,ff%s,ff%d,pc,dh,flx)

      endif

      ! summarize ice boundary changes
      if (ice%z .gt. 0) then
        snow_gl = dh(sn_precip) + dh(sn_melt) + dh(sn_subl)
        s_gl = dh(ice_s_melt) + dh(ice_s_subl)
        c_gl = dh(ice_b_melt_ml) + dh(ice_b_melt_con) + dh(ice_b_grow_ml) &
          + dh(ice_b_grow_con)

        if (nt .gt. 9000 .and. c_gl .lt. 0.) then
          c_gl = c_gl + 0.
        endif


        if ((abs(s_gl) .ge. regrid_trigger) .or. &
            (abs(c_gl) .ge. regrid_trigger) .or. &
            (dh(sn_flood) .gt. 0.) .or. &
            (dh(sn_rain) .gt. 0.) .or. &
            (abs(snow_gl) .ge. regrid_trigger) .or. &
            (jjj .eq. n_sia2_steps)) then

          ! regridding & mass/energy flux updates
          congelation_ice_volume = congelation_ice_volume &
            + ice%af * max(0.,c_gl)
          snow_ice_volume = snow_ice_volume &
            + ice%af * dh(sn_flood)*ice%snow%d(1)/max(IceD/2.d0,ice%snow%d(1))  ! potentially compress snow
          call boundary_adjust(ice,ff,dh)

        endif
      endif

    enddo

    print *, ' '
    print *, 'End of Time integration'


    ! CATEGORIES:::::::::::::::::::::::::::::::::::::::::::::::::::
    ! potential ice merge, reduce categories


    ! CATEGORIES:::::::::::::::::::::::::::::::::::::::::::::::::::
    ! Find mean ice from categories, update sums for reporting


    ! Update old ice model variables for coupling
    ! hs,hi,ni,ns,hfs,ts,ti,dzi/dzs
    ! -------------------------------------------------------------
    call ice_report(ice,nt)

    700 format (/,'old and new Tsf',2f9.3)
    701 format ('snow layer',i5,2f9.3)
    702 format ('ice layer ',i5,2f9.3)
    703 format ('layer #',i5,f9.3)

!     Calculate ending Heat, Freshwater, and Salt Content
!     (hc also used in iceflx/sflux(4,4) calc)

    if (hi > 0) then
        call ice_hfs (ice,hc_new,temp2,temp3)        ! heat content of snow-ice

        hc_ice = dble(hc_new)
        fc_ice = dble(temp2)
        sc_ice = dble(temp3)

        call sia2_hc(ice,flx(11),temp2) ! find total ice/snow heat content

        print *, ' '
        print *, 'hc [J/m2]      = ', hc_ice
        print *, 'fc_ice [kg/m2] = ', fc_ice
        print *, 'sc_ice [kg/m2] = ', sc_ice
    end if

!     Return ice-ocean fluxes

    print *, ' '
    print *, 'calling ice_flux'
    call ice_flux(ice,sflux,NSFLXS_local,hc_old,hc_new,dti,flx,X)
    print *, ' '
    print *, 'End of icestep'
    print *, ' '

  end subroutine icestep                            ! End Main Driver

!***********************************************************************


!***********************************************************************
  SUBROUTINE ice_sub_step(ice,ff,div,ic,ain,aout,sps,rps,Tair,dtt_s, &
    t_ocn,s_ocn,d_ocn,pc,dh,flx)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_constants
    use sia2_parameters
    use sia2_types
    use sia2_flux_heat
    use sia2_desalination
    use sia2_state

    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    type(ice_type), intent(inout) :: &
      ice
    type(forcing_type), intent(in) :: &
      ff
    real, intent(in) :: &
      div,            & ! SSM/I ice divergence (fraction)
      ic,             & ! SSM/I ice fraction (fraction)
      ain,            & ! SSM/I ice advection in (fraction)
      aout,           & ! SSM/I ice advection out (fraction)
      sps,            & ! snow precipitation (kg m-2 s-1)
      rps,            & ! rain precipitation (kg m-2 s-1)
      Tair,           & ! 2m air temperature
      dtt_s,          & ! ice sub-timestep length (s)
      t_ocn,          & ! ocean surface mixed layer T (degC)
      s_ocn,          & ! ocean surface mixed layer Salinity (psu)
      d_ocn             ! ocean surface mixed layer density (g m-3)
    type (heat_pre_calc), intent(in) :: &
      pc                      ! pre-calculated surface heat flux parameters
    real, intent(inout) :: &
      dh(n_dh)          ! ice thickness changes (m)
    real :: &
      flx(11)   ! diagnostic flux data structure

  ! Internal Variables
  ! --------------------------------------------------------------------
    integer :: &
      i,                  & !
      z_last,             & !
      z_ice,              & !
      sldm(z_max_ice),    & ! if == 1, SLDM desal occurred in layer
      jjj                   !

    real :: &
      AicChg,                 & ! change in total ice fraction (fraction)
      Ts_prev,                & ! current surface temperature (degC)
      Ts_next,                & ! current surface temperature (degC)
      T_diff,                 & ! running largest difference in temperature between iterated solutions (degC)
      T_step,                 & ! difference in temperature across timestep (degC)
      Fe,                     & ! latent heat flux/sublimation  (W/m^2)
      Fsens,                  & ! sensible turbulent heat flux (W/m^2)
      Flongo,                 & ! outgoing longwave (W/m^2)
      Fio,                    & ! ocean heat flux to ice  (W/m^2)
      F0,                     & ! surface heat flux, minus conductive (W/m^2)
      dF0,                    & ! derivative (w/ respect to Ts) surface heat flux, minus conductive (W/m^2/K)
      Fc_top,                 & ! surface conductive heat flux (W/m^2)
      Fc_bot,                 & ! basal conductive heat flux (W/m^2)
      Fm_bot,                 & ! mixed layer frazil production (allocated to bottom ice)
      lhcoef,                 & !
      shcoef,                 & !
      ed_w_ice(z_max_ice),    & !
      ed_w_snow(z_max_snow),  & !
      dsdt_out(z_max_ice),    & !
      lfheat(z_max_ice),      & !
      dth(z_max_pack+1),      & !
      ki(z_max_pack+1),       & !
      t_next(z_max_pack+1),   & !
      t_last(z_max_pack+1),   & !
      t_prev(z_max_pack+1)      !


  ! Function Code
  ! ------------------------ --------------------------------------------

    ! initializations
    ed_w_snow = c0    ! not used b/c no internal absorption
    ed_w_ice = c0     ! not used b/c no internal absorption

    z_last = ice%z!-z_sk  ! record z to see if ice was created
    z_last = max(0,z_last)

    ! find ocean heat flux to ice
    call sia2_ohf(t_ocn,s_ocn,d_ocn,0.,0.,0.,0.,.true.,Fio)

    ! Limit ocean heat flux to >= 2W. !f ocean model uses a different
    ! method to calculate the freezing temp (or it is fixed) ohf can go
    ! negative which is unphysical

    ! Fio = max(2.,Fio)
    Fio = max(0.,Fio)


    ice_ocean_bottom_flux_potential = Fio ! maybe not enough heat in boundary layer
    print *,'Calculated Ocean Heat Flux (Fio): ',Fio
    if (Fio .gt. 0.) then
      Fio = min(Fio,abs(Fml)) ! limit bottom melt to the mixed melting potential
    endif

    ! assume that available mixed layer heat potential is ocean heat flux,
    ! ignoring boundary heat transfer of sia2_ohf calc
    !Fio = -Fml
    !Fio = 0.

    ice_ocean_bottom_flux = Fio


    ! find potential change in ice area, make new ice, destroy ice
!    call dAic (ice,ff,div,ic,ain,aout,dtt_s,c1,c0,t_ocn,s_ocn,dh, &
!      AicChg,Fio,Fm_bot)

    if (z_last .ne. 0 .and. ice%z .gt. 0) then

      z_ice = ice%z!-z_sk

      ! prepare contiguous conductivity & midpoint thickness arrays
      call sia2_conductivity(ice%t,ice%s,ice%th,ice%d, &
        ice%snow%t,ice%snow%d,ice%snow%th,z_ice,ice%snow%z,ki,dth)

      ! desalination
      call sia2_desal(dtt_s,s_ocn,t_ocn,ice,ice%s,ed_w_ice,dth,ki, &
        sldm,dsdt_out,lfheat)
      ! lfheat = 0.  ! these line neutralize desalination for testing
      ! sldm=0
      ! dsdt_out = 0.
      ice%flux(q_totalfreeze) = ice%flux(q_totalfreeze) + sum(lfheat)*dtt_s*ice%af


      ! BEGIN CONVERGENCE ITERATIONS
      ! ----------------------------------------------------------------
      jjj = 0 ! set iteration counter to zero
      T_diff = 1.e8 ! set iteration tolerance to something ridiculously large to enter the loop
      ! find initial temperature vectors
      if (ice%snow%z .gt. 0) then
        do i=ice%snow%z,1,-1
          t_prev(ice%snow%z-i+1) = ice%snow%t(i)
        enddo
        t_prev(ice%snow%z+1:ice%snow%z+z_ice) = ice%t(1:z_ice)
      else
        t_prev(1:z_ice) = ice%t(1:z_ice)
      endif
      t_next = t_prev
      ts_prev = ice%snow%ts
      ts_next = ts_prev

      do while((abs(T_diff) .gt. nr_tol) .and. (jjj .lt. max_it))

        call sia2_F0(ts_next+kelvin0,pc,F0,dF0,Fe,Fsens,Flongo)

        ! find new temps/flux heat
        call sia2_heat_solver(ts_next,t_next,ts_prev,t_prev,ki,dth, &
          ice%s,ice%d,ice%bd,ice%bv,ice%th,ice%snow%d,ice%snow%th, &
          z_ice,ice%snow%z,ed_w_ice,ed_w_snow,lfheat,F0,dF0,t_ocn, &
          Tair-kelvin0,T_diff,T_step,Fc_top,Fc_bot,dtt_s)

         ! increment interation counter
         jjj = jjj+1

      enddo
      ! END CONVERGENCE ITERATIONS
      ! -------------------------------------------------------------

      ! update ice state, calc mass fluxes associate with desalination
      call sia2_update_ice_t_s(ice,dsdt_out,t_next,ts_next, &
        ice%flux(w_desal),ice%flux(s_desal))

      ! perform boundary and flux calculations
      call boundary_calcs(ice,Fe,dtt_s,Tair-kelvin0,sps,F0,Fio,Fc_top, &
        Fc_bot,Fm_bot,t_ocn,s_ocn,dh)

      ! flooding
      call flood(ice,d_ocn,dh(sn_flood))

      ! rain water flooding
      dh(sn_rain) = dh(sn_rain) + &
        dtt_s * rps / 1000.    ! surface snow from precipitation: [s]*[kg/m2/s]/[kg/m3] = [m]

      ! report fluxes
      call sia2_F0(ts_next+kelvin0,pc,F0,dF0,Fe,Fsens,Flongo)
      atm_flux_to_ice_surface = F0
      flx(7) = flx(7) + (pc%Fl + Flongo)*dtt_s*ice%af
      flx(6) = flx(6) + pc%Fr*dtt_s*ice%af
      flx(8) = flx(8) + Fe*dtt_s*ice%af
      flx(9) = flx(9) + Fsens*dtt_s*ice%af

    endif

    ! find potential change in ice area, make new ice, destroy ice
    call dAic (ice,ff,div,ic,ain,aout,dtt_s,c1,c0,t_ocn,s_ocn,dh, &
      AicChg,Fio,Fm_bot)


  END SUBROUTINE ice_sub_step

!***********************************************************************


!***********************************************************************
! Compute the change in ice concentration (Aic)
!
! ----------------------------------------------------------------------
  SUBROUTINE dAic (ice,ff,div,ic,ain,aout,dtt_s,Amax,Amin,t_ocn,s_ocn, &
    dh,AicChg,Fio,Fm_bot)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_constants, only: kelvin0,pi
    use sia2_parameters, only: z_sk,z_th_min,z_max_ice,z_th_fr, &
      Fm_a_switch,melt_f_denom
    use sia2_grid
    use sia2_types
    use kei_hacks, only: initial_ice_thickness,melt_lateral_complete, &
      mixed_layer_freezing_in_leads_only,ignore_divergence


    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    type(ice_type), intent(inout) :: &
      ice
    type(forcing_type), intent(in) :: &
      ff
    real :: &
      div,            & ! SSM/I ice divergence (fraction)
      ic,             & ! SSM/I ice fraction (fraction)
      ain,            & ! SSM/I ice advection in (fraction)
      aout,           & ! SSM/I ice advection out (fraction)
      dtt_s,          & !
      Amax,           & !
      Amin,           & !
      t_ocn,          & !
      s_ocn             !
    real, intent(inout) :: &
      Fio               ! ocean heat flux calculated from boundary layer theory
    real, intent(inout) :: &
      dh(n_dh)                  ! ice thickness changes
    real, intent(out) :: &
      AicChg,                 & !
      Fm_bot                    ! mixed layer frazil production (allocated to bottom ice)

  ! Internal Variables
  ! --------------------------------------------------------------------
    type(ice_type)  :: ice_new
    type(meta_type) :: m_dummy
    type(snow_type) :: new_snow_dummy
    integer :: z_snow,int_z,i,j,sk_1,sk_z,z_new
    real, dimension(z_max_ice) :: th_new, id_new
    real :: lamdam,lamdaf,dA,dAn,dAdc,denN,cc1,cc2,hc_i,hc_total,Flat, &
      Flat_used,dAe,rAAm,tmp1,hi_new,dh20,dsalt,dAlgm,gsize,dA_adv_melt, &
      dmy1,dmy2,dmy3,dmy4,Ta,hs_new,Alm,Alf,dA_melt,ic_diff,daadv,fm_frac, &
      ice_af_old


  ! Function Code
  ! ------------------------ --------------------------------------------
    !dhiem = 0.      ! equiv ice melt used to bring T to Tmelt
    !dhsem = 0.      ! equiv snow melt used to bring T to Tmelt
    !dhief = 0.      ! equiv ice freeze used to bring T to Tf
    !dhirg = 0.      ! basal ice rafting growth (3:3)
    !dhsrg = 0.      ! nonflooded-snow rafting growth (3:5)
    !dhfsrg = 0.     ! flooded-snow rafting growth (3:5)

    Aic = ice%af ! working copy of aice area

    dA = 0.        ! total change in ice conc
    dAn = 0.        ! total change in ice conc
    dAdc = 0.       ! dynamic change in Aic (2:4)
    dAlgm = 0.      ! de/increase in Aic by lateral melt/growth (2:8), unitless
    AicChg = 0.     ! Total Aic change:  Aic + dA
    Fm_bot = 0.
    dmy1 = 0.
    dmy2 = 0.
    dmy3 = 0.
    dmy4 = 0.
    Ta = ff%at - kelvin0
    new_snow_dummy%z = 0
    new_snow_dummy%depth = c0
    new_snow_dummy%th = 0
    ic_diff = 0.

! Compute the floe #-density and lat. area, denN and Al
!    DGM's eqn 14 (OI Interaction, JGR-sub)

    denN=0.         ! init number density
    cc1=2.8147e+04   ! constant for number density
    cc2=10.3487      ! constant for number density
    gsize=1         ! grid size [m]

    denN=cc1*exp(-cc2*Aic)                ! number density

    print *, 'old Aic = ', Aic

    if (Aic > 0) then

      !hi_new = ice%id(ice%z-z_sk)
      hi_new = ice%id(ice%z)
      hs_new = ice%snow%depth
      Alm=4*Aic*(hi_new+hs_new)*sqrt(denN)/gsize ! normalized Aic*lateral snow/ice area

      if ( (Fml+Fm) < 0) then            ! melt w/ sflux(8,4)

        IF (melt_lateral_complete) THEN
          Flat = -Fml
          if (Fio .gt. 0.) then
            Flat = max(0.,Flat-Fio)/Flat*Aic + Flat*(1.-Aic)
          endif

          !IF (maykut_perovich_lateral_melt) THEN
          !  Flat_used = 16.0e-7 * (S_ocn*mu - T_ocn) ** 1.36  ! lateral decrease in m s-1
          !  Flat = min(Flat,Flat_used)
          !ENDIF

        ELSE
          ! Steele 1992 edge melt
          Flat = -(Fml * (1.-Aic))* pi * ice%id(ice%z) / &
            (mean_floe_diameter * alpha_lateral_melt)
        ENDIF

        !Flat = -(Fml * (1.-Aic)) ! * Alm) / (Aic + Alm)  ! make sflux(8,4) a positive lateral heat
        !Flat = -Fml
        !Flat = max(-Fml*(1-Aic),-Fml*0.1)  ! inverse availability of mixed layer heat to melt ice, with minimun
        !Flat_used = -Flat*dtt_s  ! < 0 - heat passed back ice->ocn
        !ice%flux(q_latmelt) = ice%flux(q_latmelt) + Flat_used ! lateral flux returned to sflux(4,4), <=0

        call sia2_hc_melt(ice,t_ocn,hc_i,hc_total)      ! heat content of ice, total pack
        dAlgm = dtt_s*Flat/hc_total     ! fractional lateral decrease (hcm<0)

      elseif (Fm > 0) THEN                                ! freeze w/ -sflux(7,4)*Fl

        IF  (mixed_layer_freezing_in_leads_only) THEN
          if (Aic .ge. Fm_a_switch) then
            !Flat = (Fm * (1.-Aic))                  ! use Fm = -sflux(7,4)*Fl, >=0
            fm_frac = (1.0-Aic)/(1.0-Fm_a_switch)
            Flat = Fm * fm_frac
            fm_bot = Fm * (1.0-fm_frac) / Aic ! redirect extra cooling into bottom accretion
          else
            Flat = Fm   ! all cooling makes new ice
            ! fm_bot = 0.
          endif
        ELSE
          Flat = (Fm * (1.-Aic))                  ! use Fm = -sflux(7,4)*Fl, >=0
          fm_bot = Fm
        ENDIF

        ! -------------- new method - merging new ice
        call sia2_create_ice(ice_new,0,hi_new=initial_ice_thickness, &
          hs_new=0.,af_new=1.,Tair=Ta,t_ocn=t_ocn,s_ocn=s_ocn)
        call sia2_hc_melt(ice_new,t_ocn,hc_i,hc_total)
        dAlgm = -dtt_s*Flat/hc_i     ! fractional lateral increase (hcm<0)

        if ((dAlgm + Aic) .gt. 0.95) then ! protect against > 100% ice
          fm_frac = -dtt_s*Fm/hc_i
          dAlgm = max(0.95 - Aic,0.)
          fm_bot = Fm * (1.0 + (fm_frac-dAlgm)/fm_frac*(1.0-Aic))  ! redirect extra cooling into bottom accretion
        endif
        ice%flux(q_totalfreeze) = ice%flux(q_totalfreeze) - dAlgm*hc_i
        ! -------------- end new method - merging new ice


        ! ----------------old method
        !call sia2_hc_melt(ice,t_ocn,hc_i,hc_total)      ! heat content of ice, total pack
        !dAlgm = dtt_s*Flat/hc_total     ! fractional lateral decrease (hcm<0)
        ! ----------------old method


        !dAlgm = -dtt_s*Flat/hc_i     ! fractional lateral increase, >=0

      endif

      ! heat content of ice, total pack - used in varius locations below
      call sia2_hc_melt(ice,t_ocn,hc_i,hc_total)

      ! calculate change in ice area due to divergence/convergence
      !dAdc = div*dtt_s                ! dynamic change in Aic

      ! change in area due to inputs/outputs
      daadv = dtt_s*(ain - aout)/86400.  ! convert from area change/day to area change/timestep
      IF (ignore_divergence == 1) THEN
        daadv = 0.
      ENDIF


      ! total area change = divergence + thermal change + advection
      dAn = dAlgm + daadv
      AicChg = Aic + dAn

      ! deal with area > 1, should only happen with ic_diff = 0.
      !if (AicChg .gt. Amax) then
      if (AicChg .gt. 0.95) then
        dA = Amax - AicChg
        if (dAlgm .lt. 0.) then
          ! melt extra area
          !dAlgm = dAlgm + dA

          ! export extra area
          daadv = daadv + dA  ! dA negative
        else
          ! converge extra area
          dAdc = dAdc + dA  ! dA negative
        endif
      endif

      ! deal with conformity
      if (ic_conform .eq. 1) then

        ! find change in ice area so that is conforms to specified ice concentration
        ic_diff = ic - (Aic + dAlgm + daadv + dAdc)
        if (ic_diff .gt. 0.) then
!         if (dAlgm .lt. 0.) then
            ! if freezing, freeze more?
!           dAlgm = dAlgm + ic_diff
!         else
            ! if melting, bring in laterally
            daadv = daadv + ic_diff
!         endif
        else
!         if (dAlgm .lt. 0.) then
            ! if freezing, converge I guess
!           dAdc = dAdc + ic_diff
            ! leave laterally
            daadv = daadv + ic_diff
!         else
            ! if melting, melt
!           dAlgm = dAlgm + ic_diff
            ! leave laterally
!           daadv = daadv + ic_diff
!         endif
        endif
      endif


      ! total area change = thermal change + advection + convergence (+ correction)
      dAn = dAlgm + daadv + dAdc
      AicChg = Aic + dAn

      ! add to rafting if too much convergence
      ! Perform Rafting
      if (dAdc .lt. 0.) then

        dA = Aic + dAlgm + daadv

        if (dA .gt. 0.5) then

          rAAm = dA/(dA + dAdc)  ! height multiplier for convergence

          ! inflate layers of ice during rafting, then regrid
          sk_z = ice%z
          int_z = sk_z - z_sk
          sk_1 = int_z+1
          !hi_new = ice%id(int_z)*rAAm
          hi_new = ice%id(ice%z)*rAAm

  !       dhirg = hi_new            ! rafted ice thickness

          !tmp1 = ice%th(ice%z) ! preserve skeletal thickness
          ice%th = ice%th*rAAm
          ice%id = ice%id*rAAm
          !ice%th(sk_1:sk_z) = tmp1
          !do i=sk_1,sk_z
          ! ice%id(i) = ice%id(i-1)+tmp1
          !enddo
          call sia2_new_grid(hi_new,z_th_min,10.,2,1,z_max_ice-z_sk,.false., &
            th_new,id_new,z_new)
          ! copy skeletal layers
          !j=z_new+1
          !do i=ice%z-z_sk+1,ice%z
          ! th_new(j) = ice%th(i)
          ! id_new(j) = id_new(j-1) + th_new(j)
          ! j = j + 1
          !enddo
          call sia2_ice_remap(ice,ff,m_dummy,z_new,th_new,id_new, &
            t_flood=-1.86,s_flood=18.,flooded=0.,melt_flood=0., &
            s_gl=dmy1,c_gl=dmy2,c_gl_sal=dmy3,dhdt_cm_s=0.,f_depth=dmy4)

          ! loose an equivalent snow area to melting on contact w/ water
          call sia2_snow_mass(ice,dh20)
          ice%flux(w_snowmelt) = ice%flux(w_snowmelt) + dh20*dAdc ! proportioned by daadv

          ! take away snow melt energy from mixed layer
          ice%flux(q_latmelt) = ice%flux(q_latmelt) - (hc_total - hc_i)*dAdc
          ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + (hc_total - hc_i)*dAdc

        else
          ! convergence/ridging is not likely - assume area is exiting ice
          daadv = daadv + dAdc
          dAdc = 0.

        endif

      endif

      if (AicChg <= Amin) then  ! all snow/ice melted/advected away

        if (dAlgm .lt. 0) then
          ! melt extra below minimum - why? taken out 3/6/2014
          !dA_melt = -1.*(min(0.,AicChg)) + dAlgm
          dA_melt = dAlgm

          ! get rid off all ice and snow, accounting for mass and energy
          !sk_z = ice%z
          !int_z = sk_z+1
          !sk_1 = int_z+1

          ! add some more melting from ice advected away
          if (daadv .lt. 0.) then
            dA_adv_melt = abs(daadv)*export_melt_f
          else
            dA_adv_melt = 0.
          endif

          call sia2_ice_mass(ice,dh20,dsalt)
          ice%flux(w_latmelt) = ice%flux(w_latmelt) &
            - dh20*(dA_melt + dA_adv_melt)
          ice%flux(s_latmelt) = ice%flux(s_latmelt) &
            - dsalt*(dA_melt  + dA_adv_melt)
          call sia2_snow_mass(ice,dh20)
          ice%flux(w_snowmelt) = ice%flux(w_snowmelt) &
            - dh20*(dA_melt + dA_adv_melt) ! proportioned dy area

          ! lateral flux returned to sflux(4,4), <=0
          ice%flux(q_latmelt) = ice%flux(q_latmelt) &
             - hc_total*dA_melt  ! hc_total < 0
          ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + hc_total*dA_melt ! > 0

        endif ! end of melting

        call sia2_null_ice(ice)

      else ! entire concentration change allowed

        if (dAn .ne. 0.) then

          ! no change in snow mass ... decrease thickness accordingly & remap
          if (dAlgm .gt. 0.) then

! ------------ old method - just increase area of current ice
!            i = w_latfreeze
!            j = s_latfreeze
!            ice%snow%depth = ice%snow%depth * Aic / (dA + Aic)
!            ice%snow%th = ice%snow%th * Aic / (dA + Aic)
!            if (ice%snow%depth .lt. z_th_min) then
!              if (ice%snow%z .gt. 0) then
!                ! take away snow grid, not enough snow left to conduct heat properly
!                ice%snow%d_small = ice%snow%d(1)
!                ice%snow%t = c0
!                ice%snow%th = c0
!                ice%snow%d = c0
!                ice%snow%heat = c0
!                ice%snow%melt = c0
!                ice%snow%z = c0
!              endif
!            else
!            ! regrid reduced snow thickness
!              hs_new = ice%snow%depth
!              call sia2_new_grid(hs_new,z_th_min,5.e0,2,1,z_max_snow,.false., &
!                th_new,id_new,z_new)
!              dmy1 = c0
!              dmy2 = c0
!              call sia2_snow_remap(ice,ff,new_snow_dummy, &
!                z_new,th_new,dmy1,dmy2)
!            endif
! ------------ end old method - just increase area of current ice

! ------------ new method - merge new ice
            ! operating under assumption of no convergence... so dAn = dAlgm + daadv
            i = w_latfreeze
            j = s_latfreeze
            ice_new%af = dAlgm + dAlgm * daadv/ice%af*0.5
            ice_af_old = ice%af  ! record old ice af
            ice%af = ice%af + daadv*0.5

            IF (ice_new%af > 1.0 .or. ice_new%af < 0.0) THEN
              ice_new%af = dAlgm
              ice%af = ice_af_old
            ENDIF

            IF (ice%af > 1.0 .or. ice%af < 0.0) THEN
              ice_new%af = dAlgm
              ice%af = ice_af_old
            ENDIF

            !IF (daadv < 0) THEN
            !  ice_new%af = max(ice_new%af + daadv,c0)
            !ENDIF

            IF(ice_new%af > c0) THEN
              call sia2_merge_ice(ice_new,ice,t_ocn,s_ocn)
              call sia2_new_grid(ice%id(ice%z),z_th_min,10.,2,1,z_max_ice, &
                .false.,th_new,id_new,z_new)
              call sia2_ice_remap(ice,ff,m_dummy,z_new,th_new,id_new, &
                t_flood=-1.86,s_flood=18.,flooded=0.,melt_flood=0., &
                s_gl=dmy1,c_gl=dmy2,c_gl_sal=dmy3,dhdt_cm_s=0.,f_depth=dmy4)
              if (ice%snow%depth .lt. z_th_min) then
                if (ice%snow%z .gt. 0) then
                  ! take away snow grid, not enough snow left to conduct heat properly
                  ice%snow%d_small = ice%snow%d(1)
                  ice%snow%t = c0
                  ice%snow%th = c0
                  ice%snow%d = c0
                  ice%snow%heat = c0
                  ice%snow%melt = c0
                  ice%snow%z = c0
                endif
              ENDIF
            else
            ! regrid reduced snow thickness
              hs_new = ice%snow%depth
              call sia2_new_grid(hs_new,z_th_min,5.e0,2,1,z_max_snow,.false., &
                th_new,id_new,z_new)
              dmy1 = c0
              dmy2 = c0
              call sia2_snow_remap(ice,ff,new_snow_dummy, &
                z_new,th_new,dmy1,dmy2)
            endif

            frazil_ice_volume = frazil_ice_volume + dAlgm*initial_ice_thickness

            ice%af = ice_af_old   ! restore old area in case of further operations
! ------------ end new method - merge new ice

          else
            ! record heat flux used from melting for coupling
            Flat_used = dAlgm*hc_total
            ice%flux(q_latmelt) = ice%flux(q_latmelt) - Flat_used ! lateral flux returned to sflux(4,4), <=0
            ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + Flat_used ! > 0

            i = w_latmelt
            j = s_latmelt
          endif

          ! compute mass fluxes resulting from changes in ice area
          call sia2_ice_mass(ice,dh20,dsalt)

          if (dAlgm .lt. 0.) then
            ! record eco ice->ocn fluxes
            ice_to_ocean_eflux(Fe_ind) = ice_to_ocean_eflux(Fe_ind) &
              - ice_Fe_local*ice%id(ice%z)*dAlgm*import_melt_f
            ice_to_ocean_eflux(diatChl_ind) = ice_to_ocean_eflux(diatChl_ind) &
              - ice_diatChl*ice%id(ice%z)*dAlgm*import_melt_f
            ! account for only lateral melt/growth mass transfer
            ice%flux(i) = ice%flux(i) - dh20*dAlgm*import_melt_f
            ice%flux(j) = ice%flux(j) - dsalt*dAlgm*import_melt_f

            lateral_freshwater_melt_flux = lateral_freshwater_melt_flux -  &
                dh20*dAlgm

          else
            ! account for only lateral melt/growth mass transfer
            ice%flux(i) = ice%flux(i) - dh20*dAlgm
            ice%flux(j) = ice%flux(j) - dsalt*dAlgm
          endif

        endif

      endif

      ! update ice concentration
      Aic = Aic + dAn
      ice%af = Aic

    else                             ! Aic=0

    ! Ice growth onset via ocean frazil
      if (Fm > 0 .or. ((ic_conform .eq. 1) .and. ic .gt. c0)) then

        call sia2_create_ice(ice,0,hi_new=initial_ice_thickness, &
          hs_new=0.,af_new=1.,Tair=Ta,t_ocn=t_ocn,s_ocn=s_ocn)
        call sia2_hc_melt(ice,t_ocn,hc_i,hc_total)

        if ( Fm > 0 ) then
          ! form new ice using heat - ic will be corrected next time step if needed
          Aic = abs(Fm * dtt_s / hc_total)
          if (Aic > Amax) then
            tmp1 = z_th_fr*Aic/Amax
            tmp1 = tmp1/z_th_fr
            Aic = Amax
            !ice%th(1:ice%z-z_sk) = ice%th(1:ice%z-z_sk)*tmp1
            !ice%id(1:ice%z-z_sk) = ice%id(1:ice%z-z_sk)*tmp1
            ice%th(1:ice%z) = ice%th(1:ice%z)*tmp1
            ice%id(1:ice%z) = ice%id(1:ice%z)*tmp1
            !do i=ice%z-z_sk+1,ice%z
            ! ice%id(i) = ice%id(i-1)+ice%th(ice%z)
            !enddo
          endif

          ice%af = Aic
          frazil_ice_volume = frazil_ice_volume + Aic*initial_ice_thickness

          ! fluxes ice -> ocn (kg/m^2) (no snow yet...)
          call sia2_ice_mass(ice,dh20,dsalt)
          ice%flux(w_latmelt) = ice%flux(w_latmelt) - dh20*Aic
          ice%flux(s_latmelt) = ice%flux(s_latmelt) - dsalt*Aic

          ! correct ice concentration of created ice
          if (ic_conform .eq. 1) then
            ice%af = ic
          endif

        else

          ! ic increase due to drift - initialize to ic w/ no heat transfer
          Aic = ic
          ice%af = ic

          ! make sure ice advected in is not created with unrealistic temps
          if (Tair .gt. ice%s(1)*mu) then
            ! fix ice temp/sal for summertime advective ice
            ice%snow%ts = min(-1.,t_ocn)
            do i=1,ice%z
                ice%s(i) = 5.
                ice%t(i) = min(-1.,t_ocn)
                call sia2_ice_state(ice%t(i),ice%s(i),ice%bs(i),ice%bd(i), &
                  ice%d(i),ice%bv(i),ice%heat(i))
            enddo
            do i=1,ice%snow%z
                ice%snow%t(i) = min(-1.,t_ocn)
                ice%snow%heat(i) = sia2_snow_heat(ice%snow%d(i),ice%snow%t(i))
            enddo
          endif

        endif

      endif                  ! if Fm > 0

    endif                   ! Aic >/= 0

    print *, 'new Aic = ', Aic

  end SUBROUTINE dAic

!***************************************************************


!***********************************************************************

  subroutine sia2_update_ice_t_s(ice,dsdt,t_next,ts_next,fh2o,fsalt)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_parameters, only: den_s_switch
    use sia2_state
    use sia2_types

    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    type(ice_type), intent(inout) :: &
      ice
    real, intent(in) :: &
      dsdt(z_max_ice),        & !
      t_next(z_max_pack+1),   & !
      ts_next                   !
    real, intent(inout) :: &
      fh2o,                   & !
      fsalt                     !

  ! Internal Variables
  ! --------------------------------------------------------------------
    integer :: i
    real :: dh2o,dsalt

  ! Function Code
  ! ------------------------ --------------------------------------------

    ! update ice state
    ice%snow%ts = Ts_next
    if (ice%snow%z .gt. 0) then
      do i=1,ice%snow%z
        ice%snow%t(i) = t_next(ice%snow%z-i+1)
        if(ice%snow%t(i) .gt. den_s_switch) then
          ice%snow%melt(i) = 1.
        endif
        ice%snow%heat(i) = sia2_snow_heat(ice%snow%d(i),ice%snow%t(i))
      enddo
    endif
    do i=1,ice%z !-z_sk
      dsalt = sia2_layer_salt(ice%s(i),ice%d(i))
!     dh2o = sia2_layer_h2o(ice%s(i),ice%d(i))
      ice%t(i) = t_next(i+ice%snow%z)
      ice%s(i) = ice%s(i) + dsdt(i)
      call sia2_ice_state(ice%t(i),ice%s(i),ice%bs(i),ice%bd(i), &
        ice%d(i),ice%bv(i),ice%heat(i))
      ! ignoring small freshwater flux - overshadowed by density changes that
      ! are difficult to deal with
!     dh2o = (sia2_layer_h2o(ice%t(ii),ice%t(ii)) - dh2o)*ice%th(i)
!     fh2o = fh2o + dh2o*ice%af
      ! salt flux from salinity change  - density changes here will
      ! cause some inaccuracies .... FLAG
      dsalt = (sia2_layer_salt(ice%s(i),ice%d(i)) - dsalt)*ice%th(i)
      fsalt = fsalt + dsalt*ice%af

    enddo

  end subroutine sia2_update_ice_t_s

!***********************************************************************


!***********************************************************************

  SUBROUTINE boundary_calcs(ice,Fe,dtt_s,Tair,sps,F0,F_ocn,Fc_top, &
    Fc_bot,Fm_bot,t_ocn,s_ocn,dh)


  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_constants, only: c0,nc1,Lf,Lv,ci0,c_1e3
    use sia2_parameters, only: den_s_dry,den_s_wet,th_min
    use sia2_flux_heat
    use sia2_grid
    use sia2_state
    use sia2_types

    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    type(ice_type), intent(inout) :: &
      ice
    real, intent(in) :: &
      Fe,                       & ! latent turbulent heat flux
      dtt_s,                    & ! time step in seconds
      Tair,                     & ! degC
      sps,                      & ! precipitation (kg m-2 s-1)
      F0,                       & ! surface heat balance less conductive flux (W/m^2)
      F_ocn,                    & ! ocean heat flux (W/m^2)
      Fc_top,                   & ! surface conductive heat flux
      Fc_bot,                   & ! basal conductive heat flux
      Fm_bot,                   & ! mixed layer frazil production (allocated to bottom ice)
      t_ocn,                    & ! ocean mixed layer temperature (degC)
      s_ocn                       ! ocean mixed layer salinity (psu)
    real, intent(inout) :: &
      dh(n_dh)                  ! ice thickness changes

  ! Internal Variables
  ! --------------------------------------------------------------------
    integer :: &
      z_snow,                 & !
      int_z                     !
    real :: &
      f_salt,                 & !
      h1,                     & !
      h2,                     & !
      th_ice,                 & ! total current ice thickness
      s_new,                  & !
      F_temp,                 & !
      F_prev,                 & !
      F_mlm,                  & !
      Fbmelt,                 & !
      F_surf,                 & !
      F_bot,                  & !
      Fe_tot,                 & !
      q_melt,                 & !
      dh_melt,                & !
      dh_subl,                & !
      dh_grow,                & !
      d_snow                    !

  ! Function Code
  ! --------------------------------------------------------------------
    f_salt = c0
    z_snow = ice%snow%z
    th_ice = ice%id(ice%z)
    if(z_snow .gt. 0) then
      d_snow = ice%snow%d(z_snow)
    else
      if (Tair .lt. den_s_switch) then
         d_snow = den_s_dry*c1e6
      else
         d_snow = den_s_wet*c1e6
      endif
    endif
    !int_z = ice%z-z_sk
    !int_z = max(0,int_z)
    !ice%z = int_z  ! temporarily assign bottom to int_z - KEI does not use skeletal layers yet
    dh(ice_b_dhdt) = 0.  ! cm/s

    ! surface latent heat
    Fe_tot = Fe*dtt_s ! (J m-2)

    ! surface flux balance
    F_surf = (F0 - Fc_top)*dtt_s ! (J m-2)

    ! bottom flux balance (F_mlm > 0 is melting, Fc_bot > 0 is growing,
    !                      Fm_bot > 0 is growing)

    ! let Fm (mixed freezing layer potential) take care of non-conductive growth,
    ! otherwise heat mixed layer heat accounting can be messed up
    F_mlm = max(0.,F_ocn)
    F_bot = (Fc_bot - F_mlm + Fm_bot)*dtt_s ! (J m-2)

    print *,'========================================================='
    print *,'ICEICEICEICEICEICEICEICEICEICEICEICEICEICEICEICEICEICEICE'
    print *,'========================================================='
    print *,'Fc_bot: ',Fc_bot
    print *,'F_ocn: ',F_ocn
    print *,'Fm_bot: ',Fm_bot
    print *,'========================================================='
    print *,'ICEICEICEICEICEICEICEICEICEICEICEICEICEICEICEICEICEICEICE'
    print *,'========================================================='


    ! ------------------------------------------------------------------
    ! Snow Calcs - need to set:
    ! ------------------------------------------------------------------
    !dhssm = c0   ! snow surface melt/sublimation
    dh(sn_precip) = dh(sn_precip) + &
      dtt_s * sps / (d_snow*c_1e3)    ! surface snow from precipitation: [s]*[kg/m2/s]/[kg/m3] = [m]
    snow_precip_mass = snow_precip_mass + dtt_s * sps * ice%af
    ! ------------------------------------------------------------------
    if (ice%snow%depth .gt. c0) then
      ! surface latent heat flux
      if (Fe_tot .lt. c0) then

        F_temp = abs(Fe_tot)

        ! First try to sublime new snow - dhssps
        if (dh(sn_precip) .gt. c0) then
          q_melt = d_snow*(Lf + Lv)
          dh_subl = F_temp/q_melt
          if (dh_subl .gt. dh(sn_precip)) then
            F_temp = q_melt*(dh_subl - dh(sn_precip))
            dh_subl = dh(sn_precip)
          else
            F_temp = c0
          endif
          dh(sn_precip) = dh(sn_precip) - dh_subl
          ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + dh_subl*q_melt*ice%af
        endif

        ! Next try to sublime older snow in grid - dhssps
        h1 = dh(sn_subl) + dh(sn_melt)
        Fe_tot = F_temp
        call sia2_melt_subl_snow(ice,F_temp,h1,t_ocn,.false., &
          .true.,dh_subl)
        ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + (Fe_tot - F_temp)*ice%af
        Fe_tot = nc1*F_temp       ! setup for ice sublimation later
        dh(sn_subl) = dh(sn_subl) - dh_subl
      endif

      ! surface melt
      !if (F_surf .gt. c0) then
      if (F_surf .gt. c0 .and. ice%snow%ts .ge. c0) then

        ! First try to melt new snow - dhssps
        if (dh(sn_precip) .gt. c0) then
          q_melt = d_snow*(Lf + Lv)
          dh_melt = F_surf/q_melt
          if (dh_melt .gt. dh(sn_precip)) then
            F_surf = q_melt*(dh_melt - dh(sn_precip))
            dh_melt = dh(sn_precip)
          else
            F_surf = c0
          endif
          dh(sn_precip) = dh(sn_precip) - dh_melt
          ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + dh_melt*q_melt*ice%af
          ice%flux(w_snowmelt) = ice%flux(w_snowmelt) + &
            d_snow*dh_melt*c_001*ice%af    ! report kg melt water -> ocn
        endif
        h1 = dh(sn_subl) + dh(sn_melt)
        F_temp = F_surf
        call sia2_melt_subl_snow(ice,F_surf,h1,t_ocn,.false., &
          .false.,dh_melt)
        ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + (F_temp-F_surf)*ice%af
        dh(sn_melt) = dh(sn_melt) - dh_melt

      endif

    endif ! end of is there any snow? test


    ! below is not needed with sia2 model
    ! ------------------------------------------------------------------
    ! Flooded Snow Calcs - need to set:
    ! ------------------------------------------------------------------
    !dhfsm  = c0    ! flooded snow surface melt
    !dhfsf  = c0    ! total new flooded snow
    !dhfsbf = c0    ! flooded snow freeze / loss to ice
    ! ------------------------------------------------------------------

    ! ------------------------------------------------------------------
    ! Ice Calcs - need to set:
    ! ------------------------------------------------------------------
    !dhism  = c0        ! ice surface melt
    !dhisg  = nc1*dhfsbf    ! ice surface growth from re-freezing flooded snow
    !dhbc = c0          ! basal growth/melt from conductive heat flux
    !dhbf = c0          ! basal growth from ocean frazil
    !dhbm = c0          ! basal melt from ocean melt potential
    !dhieb = c0        ! NOT USED BY SIA2 CODE - equiv basal melt due to T adjustment
    ! ------------------------------------------------------------------

    h1 = max(c0,abs(dh(ice_s_melt) + dh(ice_s_subl)))
    h2 = max(c0,abs(dh(ice_b_melt_ml) + dh(ice_b_melt_con) + dh(ice_b_grow_ml) &
      + dh(ice_b_grow_con)))

    ! latent heat flux (sublimation of ice)
    if (Fe_tot .lt. c0 .and. th_ice .gt. th_min) then

      F_temp = abs(Fe_tot)
      F_prev = F_temp
      call sia2_melt_subl_ice(ice,F_temp,h1,h2,t_ocn,.false., &
        .false.,.true.,dh_subl,f_salt)
      ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + (F_prev-F_temp)*ice%af
      Fe_tot = nc1*F_temp
      dh(ice_s_subl) = dh(ice_s_subl) - dh_subl
      dh(ice_s_sal) = dh(ice_s_sal) + f_salt
      h1 = max(c0,abs(dh(ice_s_melt) + dh(ice_s_subl)))

    endif

    ! surface melt
    if (F_surf .gt. c0 .and. th_ice .gt. th_min) then

      F_prev = F_surf
      call sia2_melt_subl_ice(ice,F_surf,h1,h2,t_ocn,.false., &
        .false.,.false.,dh_melt,f_salt)
      ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + (F_prev-F_surf)*ice%af
      dh(ice_s_melt) = dh(ice_s_melt) - dh_melt
      h1 = max(c0,abs(dh(ice_s_melt) + dh(ice_s_subl)))

    endif

    ! basal boundary
    if (F_bot > c0) then ! basal growth

      ice%flux(q_totalfreeze) = ice%flux(q_totalfreeze) + F_bot*ice%af
      call sia2_basal_growth(ice,t_ocn,s_ocn,F_bot/dtt_s,dtt_s,dh_grow, &
        s_new,dh(ice_b_dhdt))
      ! protect against zero error/precision problems here - above can return zero growth
      if (dh_grow > c0) then
        h2 = dh(ice_b_grow_ml) + dh(ice_b_grow_con)
        dh(ice_b_sal) = (dh(ice_b_sal)*h2 + s_new*dh_grow)/(h2 + dh_grow) ! salinity of new ice
        dh(ice_b_grow_con) = dh(ice_b_grow_con) + dh_grow   ! add to new ice
        h2 = max(c0,abs(dh(ice_b_melt_ml) + dh(ice_b_melt_con) + dh(ice_b_grow_ml) &
          + dh(ice_b_grow_con)))
      endif
      F_bot = c0   ! is this right?

    elseif (th_ice .gt. th_min) then ! basal melt

      F_temp = abs(F_bot)
      call sia2_melt_subl_ice(ice,F_temp,h1,h2,t_ocn,.false., &
        .true.,.false.,dh_melt,f_salt)
      ice%flux(q_totalmelt) = ice%flux(q_totalmelt) + (-F_bot-F_temp)*ice%af
      dh(ice_b_melt_con) = dh(ice_b_melt_con) - dh_melt   ! add to new ice
      h2 = max(c0,abs(dh(ice_b_melt_ml) + dh(ice_b_melt_con) + dh(ice_b_grow_ml) &
        + dh(ice_b_grow_con)))

      F_bot = F_bot - min(c0,Fc_bot*dtt_s) + F_temp ! result is negative if there is mixed layer potential left over
            ! (-)   - (-)                  + (+)
      F_bot = min(c0,F_bot) ! if there is unused conductive melt + mixed layer melt potential

    endif

    ! Add frazil flux (Fm) to ice bottom (if ni=0, then all frazil added in dAic)
!    if (Fm_bot .gt. c0) then
!     call sia2_basal_growth(ice,t_ocn,s_ocn,Fm_bot,dtt_s,dh_grow,s_new,dh(ice_b_dhdt))
!     if (dh_grow .gt. c0) then
!       h2 = dh(ice_b_grow_ml) + dh(ice_b_grow_con)
!       dh(ice_b_sal) = (dh(ice_b_sal)*h2 + s_new*dh_grow)/(h2 + dh_grow) ! salinity of new ice
!       dh(ice_b_grow_ml) = dh(ice_b_grow_ml) + dh_grow   ! add to new ice
!       h2 = max(c0,abs(dh(ice_b_melt_ml) + dh(ice_b_melt_con) + dh(ice_b_grow_ml) &
!         + dh(ice_b_grow_con)))
!     endif
!   endif


    ! ------------------------------------------------------------------
    ! Update Energy Fluxes - need to set: Fin,Fbcon,Fbmelt
    ! ------------------------------------------------------------------
    ! update potential fluxes back to ocean, if not accounted for
    ! by boundary changes
    Fe_tot = min(Fe_tot,c0) ! latent heating accounted for in surface heat balance - deposition not accounted for yet
    if (F_surf .gt. c0) then
      F_surf = nc1*F_surf ! ice melt not accounted for
    else
      F_surf = c0 ! ice getting colder, should have been taken care of in heat transfer
    endif
    ice%flux(q_ssmelt) = ice%flux(q_ssmelt) + (Fe_tot + F_surf)*ice%af

    ice%flux(q_bcon) = ice%flux(q_bcon) + F_bot*ice%af ! Fbot < 0
    !ice%flux(q_mlmelt) = ice%flux(q_mlmelt) + Fbmelt*ice%af
    !ice%flux(q_mlmelt) = ice%flux(q_mlmelt) + Fm*dtt_s*ice%af  ! no longer adding frazil ice to bottom - passing back in ocean
    !if (fml .lt. 0.) then
    ! print *,'Fraction Fml (W/m^2) used to melt ice : ', F_bot/dtt_s, '/',Fml
    !endif

    !ce%z = int_z+z_sk  ! restore ice%z - KEI does not use skeletal layers yet

  END SUBROUTINE boundary_calcs

!***********************************************************************


!***********************************************************************

  SUBROUTINE flood(ice,d_ocn,flooded)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_parameters, only: vb_crit,z_sk
    use sia2_constants, only: c0
    use sia2_desalination
    use sia2_state
    use sia2_types

    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    type(ice_type), intent(inout) :: &
      ice
    real, intent(in) :: &
      d_ocn
    real, intent(out) :: &
      flooded

  ! Internal Variables
  ! --------------------------------------------------------------------
    logical :: &
      porous_ice = .true.   ! are brine channels open? switch
    integer :: &
      ii,                 & ! iterator
      vb_open               ! highest layer away from ocean where brien channels are open
    real :: &
      snow_mass,          & ! (g)
      ice_mass,           & ! (g)
      fbh_new                ! new freeboard height of ice (m)

    flooded = c0

    ! find mean density
    ice_mass = c0
    do ii = 1,ice%z
        ice_mass = ice_mass + ice%d(ii)*ice%th(ii)
    enddo
    ! find snow mass
    snow_mass = c0
    do ii = 1,ice%snow%z
        snow_mass = snow_mass + ice%snow%d(ii)*ice%snow%th(ii)
    enddo

    ! Calculate freeboard height of mean snow depth
    fbh_new = ice%id(ice%z) - (snow_mass + ice_mass)/d_ocn ! m

    print *, 'isostatic ice freeboard = ', fbh_new

    !vb_open = sia2_bv_open(ice%bv,ice%z-z_sk,vb_crit)
    vb_open = sia2_bv_open(ice%bv,ice%z,vb_crit)
    if (vb_open .ne. 0 .and. fbh_new .lt. -0.2) then
      print *,'ice not-porous, surface flooding not allowed: ',ti(1:ni)
      porous_ice = .false.
    endif

    if (fbh_new .le. fs_crit .and. porous_ice) then
        flooded = abs(fbh_new)
        fbh_new = fbh_new + flooded
    endif

    ice%fbh = fbh_new

  END SUBROUTINE flood

!***********************************************************************


!***********************************************************************
  subroutine boundary_adjust(ice,ff,dh)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_constants, only : &
      c0
    use sia2_parameters, only : &
      z_sk,         &
      z_th_min,     &
      z_max_ice
    use sia2_flux_heat
    use sia2_state
    use sia2_grid

    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    type(ice_type), intent(inout) :: &
      ice
    type(forcing_type), intent(in) :: &
      ff
    real, intent(inout) :: &
      dh(n_dh)          !

  ! Internal Variables
  ! --------------------------------------------------------------------
    type(meta_type) :: &
      m_dummy           ! dummy meta type for feeding sia2 functions
    type(snow_type) :: &
      new_snow          ! snow type to hold new "snow" created from drained ice
    logical  :: &
      did_flood,      & !
      flood_1x,       & !
      did_drain         !
    integer :: &
      i,              & ! iterator
      j,              & ! iterator
      z_new             !
    real :: &
      th_diff,          & ! change in total ice thickness
      th_new(z_max_ice),          & !
      id_new(z_max_ice),          & !
      f_temp,         & !
      fh20,           & !
      fsalt,          & !
      dh_melt,        & !
      temp,           & ! dummy flux used in sia2_melt_subl_snow
      snow_gl,        & !
      h1,             & !
      h2,             & !
      s_gl,           & !
      c_gl,           & !
      c_gl_sal,       & !
      flooded,        & !
      rain_fraction,  & !
      sn_depth,       & !
      dhdt_cm_s,      & !
      r_depth,        & !
      f_depth,        & !
      t_flood,        & !
      sal_flood,      & !
      melt_flood,     & !
      melt_drain,     & !
      dh20,           & !
      d_snow,         & !
      dsalt,          & !
      q_sum

  ! Function Code
  ! ------------------------ --------------------------------------------

    ! initializations
    if ((ff%at-kelvin0) .lt. den_s_switch) then
       d_snow = den_s_dry*c1e6
    else
       d_snow = den_s_wet*c1e6
    endif
    snow_gl = dh(sn_precip) + dh(sn_melt) + dh(sn_subl)
    s_gl = dh(ice_s_melt) + dh(ice_s_subl)
    c_gl = dh(ice_b_melt_ml) + dh(ice_b_melt_con) + dh(ice_b_grow_ml) &
      + dh(ice_b_grow_con)
    c_gl_sal = dh(ice_b_sal)
    flooded = dh(sn_flood) + dh(sn_rain)
    if (flooded > 0.) then
      rain_fraction = dh(sn_rain)/flooded
    else
      rain_fraction = 0.
    endif
    dhdt_cm_s = dh(ice_b_dhdt)

    call sia2_pre_grid(ice,ff,m_dummy,flooded,rain_fraction,c_gl,s_gl, &
      snow_gl,r_depth,f_depth,did_flood,flood_1x,t_flood,sal_flood, &
      melt_flood,melt_drain,did_drain,sn_depth,new_snow)

    ! record eco ice->ocn fluxes
    th_diff = r_depth - ice%id(ice%z)
    if (th_diff .lt. c0) then
      th_diff = abs(th_diff)  ! sources, not sinks!
      ice_to_ocean_eflux(Fe_ind) = ice_to_ocean_eflux(Fe_ind) &
        + ice_Fe_local*th_diff*ice%af
      ice_to_ocean_eflux(diatChl_ind) = ice_to_ocean_eflux(diatChl_ind) &
        + ice_diatChl*th_diff*ice%af
    endif

    ! wipe out ice if less than 1 layer thick, return melt fluxes to mixed layer
    if (r_depth .lt. z_th_min) then

      if (r_depth .gt. c0) then
        ! below setting ignore_f=.true. updates energy fluxes as well
        F_temp = c0
        h1 = abs(min(c0,dh(sn_melt) + dh(sn_subl)))
        call sia2_melt_subl_snow(ice,temp,h1,ff%t,.true., &
            .false.,dh_melt)
        dh(sn_melt) = dh(sn_melt) - dh_melt
        h1 = abs(min(c0,s_gl+melt_drain))  ! adding back in melt_drain here b/c has not been accounted for in flux yet
        h2 = abs(min(c0,c_gl))
        call sia2_melt_subl_ice(ice,temp,h1,h2,ff%t,.true., &
            .true.,.false.,dh_melt,dh(ice_s_sal))
        dh(ice_b_melt_ml) = dh(ice_b_melt_ml) - dh_melt
      endif

      ! kill ice structure - fluxes accounted for in boundary_calcs
      call sia2_null_ice(ice)

    else
    ! regrid ice
      call sia2_new_grid(r_depth,z_th_min,10.,2,1,z_max_ice, &
        .false.,th_new,id_new,z_new)
      ! copy skeletal layers
!     j=z_new + 1
!     do i=ice%z-z_sk+1,ice%z
!       th_new(j) = ice%th(i)
!       id_new(j) = id_new(j-1) + th_new(j)
!       j = j + 1
!     enddo
      call sia2_ice_remap(ice,ff,m_dummy,z_new,th_new,id_new,t_flood, &
        sal_flood,flooded,melt_flood,s_gl,c_gl,c_gl_sal,dhdt_cm_s, &
        f_depth)

      ! regrid snow
      if ((ice%snow%depth .ne. sn_depth) .or. did_flood .or. did_drain) then

        ! new snow thickness is below minimum
        if (sn_depth .gt. c0) then

          if (sn_depth .lt. z_th_min) then

            ! convert discrete ice layer(s) to thin layer
            if (ice%snow%z .gt. 0) then
              ice%snow%d_small = ice%snow%d(1)
              ice%snow%depth = sn_depth
              ice%snow%t = c0
              ice%snow%d = c0
              ice%snow%heat = c0
              ice%snow%melt = c0
              ice%snow%z = c0
            else

            ! add/subtract from thin layer
              temp = sn_depth - ice%snow%depth
              if (sn_depth .gt. ice%snow%depth) then
                ice%snow%d_small = (ice%snow%d_small*ice%snow%depth + &
                  temp*d_snow)/sn_depth
              endif
              ice%snow%depth = sn_depth
            endif

          else

            ! convert thin layer -> dicrete layer snow
            if (ice%snow%z .eq. 0) then
              ice%snow%t(1) = ice%snow%ts
              ice%snow%d(1) = ice%snow%d_small
              ice%snow%heat(1) = sia2_snow_heat(ice%snow%d(1),ice%snow%t(1))
              ice%snow%melt(1) = c0
              ice%snow%z = 1
              ice%snow%th(1) = ice%snow%depth
            endif

            ! regrid snow
            call sia2_new_grid(sn_depth,z_th_min,5.,2,1,z_max_snow, &
              .false.,th_new,id_new,z_new)
            call sia2_snow_remap(ice,ff,new_snow,z_new,th_new,flooded, &
              melt_drain)

          endif ! end of < z_th_min? test

        else

          ! zero snow
          if (ice%snow%z .gt. 0) then
            ice%snow%t = c0
            ice%snow%d = c0
            ice%snow%heat = c0
            ice%snow%melt = c0
            ice%snow%z = c0
            !ice%snow%ts = ice%t(1) ! move surface temp to prevent numerical heat issues - needed?
          endif
          ice%snow%depth = c0
          ice%snow%d_small = c0

        endif

      endif

      ! reset dh vars
      dh = c0  ! vector operation

    endif

  end subroutine boundary_adjust

!***********************************************************************


!***************************************************************

  SUBROUTINE ice_flux(ice,sflux,NSFLXS_local,hc_old,hc_new,dti,flx,X)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_types

    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    type(ice_type), intent(inout) :: &
      ice               ! ice
    real, dimension (NSFLXS_local,5), intent(inout) :: &
      sflux             ! external flux matrix
    integer, intent(in) :: &
      NSFLXS_local      ! sflux parameter
    real, intent(inout) :: &
      hc_old,         & !
      hc_new            !
    real, intent(in) :: &
      dti               !
    real :: &
      flx(11)   ! diagnostic flux data structure
    real, intent(in) :: &
      X(NZP1,NSCLR)

  ! Internal Variables
  ! --------------------------------------------------------------------
    real :: &
      Fin,            & !
      Fbmelt,         & !
      Fbcon,          & !
      Flmf,           & !
      fsalt,          & !
      fh2o,           & !
      dice = -1.8,    & ! solar extinction depth for ice (m-1)
      dsnw = -15.,    & ! solar extinction depth for snow (m-1)
      QSNW              !
    double precision :: &
      ratio,            & !
      dHCice            !

  ! Function Code
  ! ------------------------ --------------------------------------------


    ! debugging print statements
      print *,' '
      print *,'w_smelt: ',ice%flux(w_smelt)/dti
      print *,'w_bmelt: ',ice%flux(w_bmelt)/dti
      print *,'w_bfreeze: ',ice%flux(w_bfreeze)/dti
      print *,'w_flood: ',ice%flux(w_flood)/dti
      print *,'w_frfreeze: ',ice%flux(w_frfreeze)/dti
      print *,'w_latmelt: ',ice%flux(w_latmelt)/dti
      print *,'w_latfreeze: ',ice%flux(w_latfreeze)/dti
      print *,'w_snowmelt: ',ice%flux(w_snowmelt)/dti
      print *,'w_desal: ',ice%flux(w_desal)/dti
      print *,'s_smelt: ',ice%flux(s_smelt)/dti
      print *,'s_bmelt: ',ice%flux(s_bmelt)/dti
      print *,'s_bfreeze: ',ice%flux(s_bfreeze)/dti
      print *,'s_flood: ',ice%flux(s_flood)/dti
      print *,'s_frfreeze: ',ice%flux(s_frfreeze)/dti
      print *,'s_latmelt: ',ice%flux(s_latmelt)/dti
      print *,'s_latfreeze: ',ice%flux(s_latfreeze)/dti
      print *,'s_desal: ',ice%flux(s_desal)/dti
      print *,'q_ssmelt: ',ice%flux(q_ssmelt)/dti
      print *,'q_mlmelt: ',ice%flux(q_mlmelt)/dti
      print *,'q_bcon: ',ice%flux(q_bcon)/dti
      print *,'q_latmelt: ',ice%flux(q_latmelt)/dti
      print *,'q_totalmelt: ',ice%flux(q_totalmelt)/dti
      print *,'q_totalfreeze: ',ice%flux(q_totalfreeze)/dti


    ! assign bookeeping vars
    total_ice_melt = ice%flux(q_totalmelt)
    total_ice_freeze = ice%flux(q_totalfreeze)

    ! ice to ocean fluxes: sflux(n,4),n=1,6

    ! momentum
    sflux(1,4) = sflux(1,3)
    sflux(2,4) = sflux(2,3)

    ! shortwave irradiance
    QSNW    = (1.-albice)* sflux(3,3) * exp( hsn(iold)/DSNW)
    sflux(3,4) = QSNW            * exp(hice(iold)/DICE)
     sflux(3,4) = 0.0

    ! let some light through 'fake' ice
    if (ic_conform .eq. 2 .and. fice .eq. 0.) then
        QSNW    = (1.-albice)* sflux(3,3) * exp( 0.05/DSNW)
        sflux(3,4) = QSNW  * exp(0.3/DICE)
    endif

    flx(6) = flx(6)/dti
    flx(7) = flx(7)/dti
    flx(8) = flx(8)/dti
    flx(9) = flx(9)/dti

    ! IO Heat Flux
    ! Fin: amt of net Famt not used in surface melting (i.e., if all snow/ice melted)
    ! Fbmelt: amt of sflux(8,4) used to melt ice/snow from bottom, <=0
    ! Fbcon:  amt of surface cooling not used to grow basal ice, <=0
    ! Flmf:   amt of lateral HF used in melt, <=0
    Fin = ice%flux(q_ssmelt)/dti
    Fbmelt = ice%flux(q_mlmelt)/dti
    Fbcon = ice%flux(q_bcon)/dti
    Flmf = ice%flux(q_latmelt)/dti

    sflux(4,4) = Fin + Fbmelt + Fbcon + Flmf          ! [W/m2]
    flx(10) = sflux(4,4)

    print *, ' '
    print *, 'Fin    = ', Fin
    print *, 'Fbmelt = ', Fbmelt
    print *, 'Fbcon  = ', Fbcon
    print *, 'Flmf   = ', Flmf
    print *, 'Total (4,4) = ', sflux(4,4)

!      ( i.e., remaining net heat flux to ocean after snow/ice melt )

!     Another angle:  budgeting the IO Heat Flux

    dHCice = hc_new - hc_old            ! change in heat content [J/m2]

    Qsfc = (1.-albice)*sflux(3,3)       & !  net solar, [W/m2]
    + sflux(4,3) + sflux(9,3)      & !  net longwave, [W/m2]
    + sflux(5,3)                   & !  sensible, [W/m2]
    + sflux(6,3) * SL              ! sublimation, [W/m2]

    print *, ' '
    print *, 'dHCice/dtflx = ', dHCice/dti
!     print *, 'dhm          = ', dhm
!     print *, 'hice(iold)   = ', hice(iold)
!     print *, 'hice(inew)   = ', hice(inew)
!     print *, '(dhm*FL*CPice)/dtflx = ', (dhm*FL*CPice)/dtflx
!     print *, 'budgeted sflux(4,4)  = ', budg44


    if (fice .gt. 0.) then

      ! MASS FLUX -> OCN -------------------------------------------------
      fh2o =    ice%flux(w_snowmelt) +    & !
                ice%flux(w_smelt) +       & !
                ice%flux(w_bmelt) +       & !
                ice%flux(w_bfreeze) +     & !
                ice%flux(w_flood) +       & !
                ice%flux(w_frfreeze) +    & !
                ice%flux(w_latmelt) +     & !
                ice%flux(w_latfreeze) +   & !
                ice%flux(w_desal)           !

      fsalt =   ice%flux(s_smelt) +       & !
                ice%flux(s_bmelt) +       & !
                ice%flux(s_bfreeze) +     & !
                ice%flux(s_flood) +       & !
                ice%flux(s_frfreeze) +    & !
                ice%flux(s_latmelt) +     & !
                ice%flux(s_latfreeze) +   & !
                ice%flux(s_desal)           !
      fsalt = fsalt/dti ! kg m-2 s-1
      fh2o = fh2o/dti   ! kg m-2 s-1

      ! if melting, use scale factor
      if (fh2o > 0.) then
        fh2o = fh2o*import_melt_f
        fsalt = fsalt*import_melt_f
      endif

      sflux(5,4) = ( fsalt - sflux(9,4) )
      print *,'fsalt, sflux(9,4), sflux(5,4): ',fsalt, sflux(9,4), sflux(5,4)
      sflux(6,4) = ( fh2o - sflux(7,4) )
      print *,'fh2o, sflux(7,4), sflux(6,4): ',fh2o, sflux(7,4), sflux(6,4)

    else ! no ice case

      ! hmm, what about the fice partitioning in o2ice? I guess it's zero anyway...
      sflux(5,4) = - sflux(9,4)
      sflux(6,4) = - sflux(7,4)

    endif

    ! convert eco fluxes
    ice_to_ocean_eflux(Fe_ind) = ice_to_ocean_eflux(Fe_ind)*1.0e-3 / dti  ! µMoles/m^3 * 1 mMole / 1000 µMoles

    ice_to_ocean_eflux(diatChl_ind) = ice_to_ocean_eflux(diatChl_ind) / dti

    ratio =  X(1,diatC_ind+2) / X(1,diatChl_ind+2)
    ice_to_ocean_eflux(diatC_ind) = &
      ice_to_ocean_eflux(diatChl_ind) * ratio  ! ? / s

    ratio =  X(1,diatSi_ind+2) / X(1,diatChl_ind+2)
    ice_to_ocean_eflux(diatSi_ind) = &
      ice_to_ocean_eflux(diatChl_ind) * ratio ! ? / s

    ratio =  X(1,diatFe_ind+2) / X(1,diatChl_ind+2)
    ice_to_ocean_eflux(diatFe_ind) = &
      ice_to_ocean_eflux(diatChl_ind) * ratio ! ? / s



    ! done with icestep and flux reporting to sflux - reset all fluxes
    ice%flux = 0.

  end SUBROUTINE ice_flux

!***************************************************************


! **********************************************************************
! SUBROUTINE: ice_hfs
! finds total heat content, freshwater content, and salt content of
! ice/snow pack
! ----------------------------------------------------------------------
  subroutine ice_hfs(ice,hc_ice,fc_ice,sc_ice)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_state
    use sia2_types

    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    type (ice_type), intent(in) :: &
      ice             ! ice structure
    real, intent(out) :: &
      hc_ice,          & ! heat content of ice and snow (J relative to 0 deg C)
      fc_ice,          & ! fresh water content of ice and snow(kg/m^2)
      sc_ice             ! salt content of ice (kg/m^2)

  ! Internal Variables
  ! --------------------------------------------------------------------
    real :: &
      temp              ! dunny var for sia2_hc subroutine call

    call sia2_hc(ice,temp,hc_ice) ! temp returned is ice hc w/out snow
    call sia2_ice_mass(ice,temp,sc_ice) ! temp returned here is ice fresh water
    call sia2_snow_mass(ice,fc_ice)
    fc_ice = fc_ice + temp  ! sum ice and snow freshwater (kg/m^2)

  end subroutine ice_hfs

! **********************************************************************

! **********************************************************************
! SUBROUTINE: ice_hfs
! finds total heat content, freshwater content, and salt content of
! ice/snow pack
! ----------------------------------------------------------------------
  subroutine ice_hfs_pack(ice_pack,hc_ice,fc_ice,sc_ice)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_state
    use sia2_types

    implicit none

  ! Function Arguments
  ! --------------------------------------------------------------------
    type (ice_pack_type), intent(in) :: &
      ice_pack           ! ice structure
    real, intent(out) :: &
      hc_ice,          & ! heat content of ice and snow (J relative to 0 deg C)
      fc_ice,          & ! fresh water content of ice and snow(kg/m^2)
      sc_ice             ! salt content of ice (kg/m^2)

    hc_ice = ice_pack%hc_ice_mean + ice_pack%hc_snow_mean
    fc_ice = ice_pack%ice_mass_mean + ice_pack%snow_mass_mean
    sc_ice = ice_pack%s_mass_mean

  end subroutine ice_hfs_pack

! **********************************************************************


! **********************************************************************
! SUBROUTINE: summarize_ice_pack
! ----------------------------------------------------------------------
  SUBROUTINE summarize_ice_pack(ice_pack)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
     USE sia2_parameters
     USE sia2_types
     USE sia2_state

  ! Function Arguments
  ! --------------------------------------------------------------------
     type (ice_pack_type), intent(inout) :: &
      ice_pack             ! ice structure

  ! Internal Variables
  ! --------------------------------------------------------------------
     integer :: &
      i,j
     real :: &
      hc_ice, &             ! ice heat
      hc_total, &           ! total snow+ice heat
      fc_ice, &             ! fresh ice mass total
      fc_snow, &            ! fresh snow mass total
      sc_ice, &             ! salt mass total
      af, &                 ! temporary ice volume
      th_af_sum, &          ! volume divisor for calculating means
      th_snow_af_sum, &      ! volume divisor for calculating means
      s_mean, &             ! mean bulk salinity (psu)
      d_mean, &             ! mean ice density
      th_mean, &            ! mean ice thickness
      q_mean, &             ! mean ice enthalpy (J/g)
      ts_mean, &            ! mean ice or snow surface temperature
      t1_mean, &            ! mean snow-ice interface temperature
      smalg_sum, &          ! total algal biomass
      af_total, &           ! area fraction - records fraction of grid cell that this ice type represents
      ff_mean, &            ! area fraction of ice with surface flooding
      d_snow_mean, &        ! mean snow density
      th_snow_mean, &       ! mean snow thickness
      q_snow_mean,  &       ! mean ice enthalpy (J/g)
      hc_snow_mean,  &      ! mean heat of snow melt
      hc_ice_mean,  &       ! mean heat of ice melt (J/g)
      s_mass_mean, &        ! mean salt (kg m-2)
      ice_mass_mean,  &     ! mean ice mass (J/g)
      snow_mass_mean        ! mean snow mass (J/g)
    real(kind=dbl_kind), dimension(n_flx) :: &
      flux

    th_af_sum = c0
    s_mean  = c0
    d_mean  = c0
    q_mean = c0
    smalg_sum = c0
    d_snow_mean = c0
    q_snow_mean = c0
    af_total = c0
    ts_mean = c0
    t1_mean = c0
    ff_mean = c0
    hc_ice_mean = c0
    hc_snow_mean = c0
    s_mass_mean = c0
    ice_mass_mean = c0
    snow_mass_mean = c0
    th_mean = c0
    th_snow_mean = c0

    do i=1,n_max_floes
      if (ice_pack%ice(i)%af .gt. c0) then

         ! record ice state

        af = ice_pack%ice(i)%af
        af_total = af_total + af
        ts_mean = ts_mean + ice_pack%ice(i)%ts * af
        t1_mean = t1_mean + ice_pack%ice(i)%t(1) * af
        if (ice_pack%ice(i)%bv(1) .gt. bv_conv) then
            ff_mean = ff_mean + af
        endif

        call sia2_hc(ice_pack%ice(i),hc_ice,hc_total)
        call sia2_ice_mass(ice_pack%ice(i),fc_ice,sc_ice)
        call sia2_snow_mass(ice_pack%ice(i),fc_snow)

        hc_ice_mean = hc_ice_mean + hc_ice * af
        hc_snow_mean = hc_snow_mean + (hc_total - hc_ice) * af
        s_mass_mean = s_mass_mean + sc_ice * af        ! mean salt (kg m-2)
        ice_mass_mean = ice_mass_mean + fc_ice * af     ! mean ice mass (J/g)
        snow_mass_mean = ice_mass_mean + fc_snow * af       ! mean snow mass (J/g)

        th_mean = th_mean * ice_pack%ice(i)%id(ice_pack%ice(i)%z) * af
        th_snow_mean = th_snow_mean * ice_pack%ice(i)%snow%depth * af

        do j=1,ice_pack%ice(i)%z

          s_mean = s_mean + ice_pack%ice(i)%s(j) *  ice_pack%ice(i)%th(j) * af
          d_mean = d_mean + ice_pack%ice(i)%d(j) *  ice_pack%ice(i)%th(j) * af
          q_mean = q_mean + ice_pack%ice(i)%heat(j) *  ice_pack%ice(i)%th(j) * af
          smalg_sum = smalg_sum + ice_pack%ice(i)%smalg(j) * &
            ice_pack%ice(i)%bv(j) * c_001 *ice_pack%ice(i)%th(j) * af

          th_af_sum = th_af_sum + ice_pack%ice(i)%th(j) * af

        enddo

        if (ice_pack%ice(i)%snow%z .gt. 0) then
        do j=1,ice_pack%ice(i)%z

          d_snow_mean = d_snow_mean + ice_pack%ice(i)%snow%d(j) *  ice_pack%ice(i)%snow%th(j) * af
          q_snow_mean = q_snow_mean + ice_pack%ice(i)%snow%heat(j) *  ice_pack%ice(i)%snow%th(j) * af

          th_snow_af_sum = th_snow_af_sum + ice_pack%ice(i)%snow%th(j) * af

        enddo
        endif

        ! record fluxes


      endif

    enddo

    ice_pack%s_mean = s_mean / th_af_sum
    ice_pack%d_mean = d_mean / th_af_sum
    ice_pack%q_mean = q_mean / th_af_sum
    ice_pack%t_mean = sia2_ice_temp(ice_pack%s_mean,ice_pack%d_mean,ice_pack%q_mean)

    ice_pack%smalg_mean = smalg_sum / th_af_sum
    ice_pack%smalg_sum = smalg_sum

    ice_pack%d_snow_mean = d_snow_mean / th_snow_af_sum
    ice_pack%q_snow_mean = q_snow_mean / th_snow_af_sum
    ice_pack%t_snow_mean = sia2_snow_temp(ice_pack%d_mean,ice_pack%q_mean)

    ice_pack%af_total = af_total
    ice_pack%ts_mean = ts_mean / af_total
    ice_pack%t1_mean = t1_mean / af_total
    ice_pack%ff_mean = ff_mean

    ice_pack%hc_ice_mean = hc_ice_mean / af_total
    ice_pack%hc_snow_mean = hc_snow_mean / af_total
    ice_pack%s_mass_mean = s_mass_mean / af_total
    ice_pack%ice_mass_mean = ice_mass_mean / af_total
    ice_pack%snow_mass_mean = snow_mass_mean / af_total

    ice_pack%th_mean = th_mean  / af_total
    ice_pack%th_snow_mean = th_snow_mean  / af_total

  end subroutine summarize_ice_pack

! **********************************************************************
! SUBROUTINE: ice_report
! translate sia2 ice back to kpp_eco_ice common varibles
! ----------------------------------------------------------------------
  subroutine ice_report(ice,nt)

  ! Use Statements and Variables/Globals/Parameters
  ! --------------------------------------------------------------------
    use sia2_constants, only: &
      c0,           & ! zero
      kelvin0         !
    use sia2_parameters, only: &
      z_sk           !
    use sia2_types

    interface
     real pure function qsat(mode, TK)
     implicit none
     integer, intent(in) :: mode
     real, intent(in) :: TK
     end function qsat
    end interface

  ! Function Arguments
  ! --------------------------------------------------------------------
    type (ice_type), intent(in) :: &
      ice             ! ice structure
    integer :: &
      nt              ! timestep
  ! Internal Variables
  ! --------------------------------------------------------------------
  ! sw albedos : ice-solid, ice-melt, snow-solid, snow-melt
    integer :: &
      ialb,i,j
    real :: &
      ! albsw(4) = (/ 0.75     , 0.66    , 0.85      , 0.75     /), & -- old values
      albsw(4) = (/ 0.41     , 0.41    , 0.81      , 0.75     /), &  ! -- this and below from Brant et al. 2005
      albsw_cloudy(4) = (/ 0.45     , 0.45    , 0.87      , 0.82     /), &
      t_switch = -1.0

  ! Function Code
  ! --------------------------------------------------------------------

    ! output variables
    !-------------------------------------------------------------------

    ! init and clear out old stuff
    hi = c0
    hs = c0
    dzi = c0
    ti = c0
    si = c0
    dzs = c0
    ts = c0
    ns = 0

    !ni = max(0,ice%z-z_sk)
    ni = max(0,ice%z)
    if (ni .gt. 0) then

      ! ice
      hi = ice%id(ni)
      dzi = c0
      dzi(1:ni) = ice%th(1:ni)
      ti(1:ni) = ice%t(1:ni)
      si(1:ni) = ice%s(1:ni)
      ti = ti + kelvin0

      ! snow
      ns = ice%snow%z
      hs = ice%snow%depth
      if (ns .gt. 0) then
        do i = ns,1,-1
          j = ns-i+1
          dzs(i) = ice%snow%th(j)
          ts(i) = ice%snow%t(j)
        enddo
        ts = ts + kelvin0
      endif

      ! update external ice/snow layer numbers
      ni_cur = ni
      ns_cur = ns

    endif

    ! surface temp
    Tas = ice%snow%ts + kelvin0
    TI0 = ice%snow%ts

    ! flooded snow not used
    hfs = c0

    ! variables by parts besides ice model
    !-------------------------------------------------------------------

    ! update main common varibles
    qI0           = QSAT(1,C2k+TI0)    ! qsat in fluxes.f; mode 1 is for ice
    hsn(inew)     = hs
    hice(inew)    = hi
    hfsnow(inew)  = hfs
    fice          = ice%af
    dhi           = hice(inew) - hice(iold)
    dhs           = hsn(inew) - hsn(iold)
    dhfs          = hfsnow(inew) - hfsnow(iold)

    ! set albedo
    if ( hsn(inew) > 0.0 ) then
!        if ( TI0 > t_switch .or. nt .gt. 6000) then
        if ( TI0 > t_switch) then
            ialb = 4
        else
            ialb = 3
        endif
    else
!        if ( TI0 >= t_switch .or. nt .gt. 6000) then
        if ( TI0 > t_switch) then
            ialb = 2
        else
            ialb = 1
        endif
    endif
    albice = albsw_cloudy(ialb)

  end subroutine ice_report

! **********************************************************************
!
!! **********************************************************************
!! SUBROUTINE: ice_report
!! translate sia2 ice back to kpp_eco_ice common varibles
!! ----------------------------------------------------------------------
!  subroutine ice_report_pack(ice_pack,nt)
!
!  ! Use Statements and Variables/Globals/Parameters
!  ! --------------------------------------------------------------------
!    use sia2_constants, only: &
!      c0,           & ! zero
!      kelvin0         !
!    use sia2_parameters, only: &
!      z_sk           !
!    use sia2_types
!
!    interface
!     real pure function qsat(mode, TK)
!     implicit none
!     integer, intent(in) :: mode
!     real, intent(in) :: TK
!     end function qsat
!    end interface
!
!  ! Function Arguments
!  ! --------------------------------------------------------------------
!    type (ice_pack_type), intent(in) :: &
!      ice             ! ice structure
!    integer :: &
!      nt              ! timestep
!  ! Internal Variables
!  ! --------------------------------------------------------------------
!  ! sw albedos : ice-solid, ice-melt, snow-solid, snow-melt
!    integer :: &
!      ialb,i,j
!    real :: &
!      albsw(4) = (/ 0.75     , 0.66    , 0.85      , 0.75     /), &
!      t_switch = -2.0
!
!  ! Function Code
!  ! --------------------------------------------------------------------
!
!    ! output variables
!    !-------------------------------------------------------------------
!
!    ! init and clear out old stuff
!    hi = c0
!    hs = c0
!    dzi = c0
!    ti = c0
!    si = c0
!    dzs = c0
!    ts = c0
!    ns = 0
!
!    !ni = max(0,ice%z-z_sk)
!    ni = max(0,ice%z)
!    if (ni .gt. 0) then
!
!      ! ice
!      hi = ice%id(ni)
!      dzi = c0
!      dzi(1:ni) = ice%th(1:ni)
!      ti(1:ni) = ice%t(1:ni)
!      si(1:ni) = ice%s(1:ni)
!      ti = ti + kelvin0
!
!      ! snow
!      ns = ice%snow%z
!      hs = ice%snow%depth
!      if (ns .gt. 0) then
!        do i = ns,1,-1
!          j = ns-i+1
!          dzs(i) = ice%snow%th(j)
!          ts(i) = ice%snow%t(j)
!        enddo
!        ts = ts + kelvin0
!      endif
!
!      ! update external ice/snow layer numbers
!      ni_cur = ni
!      ns_cur = ns
!
!    endif
!
!    ! surface temp
!    Tas = ice%snow%ts + kelvin0
!    TI0 = ice%snow%ts
!
!    ! flooded snow not used
!    hfs = c0
!
!    ! variables by parts besides ice model
!    !-------------------------------------------------------------------
!
!    ! update main common varibles
!    qI0           = QSAT(1,C2k+TI0)    ! qsat in fluxes.f; mode 1 is for ice
!    hsn(inew)     = hs
!    hice(inew)    = hi
!    hfsnow(inew)  = hfs
!    fice          = ice%af
!    dhi           = hice(inew) - hice(iold)
!    dhs           = hsn(inew) - hsn(iold)
!    dhfs          = hfsnow(inew) - hfsnow(iold)
!
!    ! set albedo
!    if ( hsn(inew) > 0.0 ) then
!!        if ( TI0 > t_switch .or. nt .gt. 6000) then
!        if ( TI0 > t_switch) then
!            ialb = 4
!        else
!            ialb = 3
!        endif
!    else
!!        if ( TI0 >= t_switch .or. nt .gt. 6000) then
!        if ( TI0 > t_switch) then
!            ialb = 2
!        else
!            ialb = 1
!        endif
!    endif
!    albice = albsw(ialb)
!
!  end subroutine ice_report_pack
!
!! **********************************************************************



end module kei_ice
