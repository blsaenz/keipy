
module sia2_flux_heat

	use sia2_constants, only : log_kind,int_kind,real_kind,dbl_kind
	
	implicit none

	type, public :: heat_pre_calc
			real(kind=dbl_kind) :: &
					Tair, 					&	! 2m air temperature (kelvin)
					potT, 					&	! potential air temperature
					shum, 					&	! specific humidity (kg m-3)
					rhum, 					&	! relative humidity (percentage)
					rho_air, 				&	! air density (kg m-3)
					Fe_sub1, 				&	! 1st static part of turbulent heat flux: J/kg/K * kg/m^3 * (dimensionless) * m/s * 1/mbar = J/m^2/mbar/K
					Fe_sub2, 				&	! 2nd static part of turbulent heat flux
					Fs_sub, 				&	! static part of sensible latent heat flux (Fs)
					v10, 						&	! 10m wind speed magnitude (m s-1)
					eps0,						&	! surface emissivity
					pres,						&	! surface pressure
					f_snow,					&	! surface fraction covered by snow (fraction)
					Fr,							&	! surface downward shortwave irradiance (W m-2)
					Fl, 						& ! surface downward longwave irradiance (W m-2)		
					F0_constants,		&
					F0_4,						&
					F0_3,						&
					F0_2,						&
					F0_1
	end type heat_pre_calc


	public 
	
	contains

! **********************************************************************
! subroutine sia2_conductivity 
! Purpose: find thermal conductivity (ki)and inter-layer midpoint distances
! (dth) 
! NOTE - SNOW PARAMETERS ARE ASSUMED uPSIDE DOWN, WITH THE BOTTOM
! SNOW LAYER AS SNOW_TH(1), SNOW_D(1), etc...
! ----------------------------------------------------------------------

	subroutine sia2_conductivity(ice_t,ice_s,ice_th,ice_d, &
		snow_t,snow_d,snow_th,ice_z,snow_z,ki,dth)

 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_parameters, only : &
			z_max_ice,						&	! maximum ice layers
			z_max_snow,						&	! maximum snow layers
			z_max_pack,						&	! maximum ice+snow layers
			ki_min									! minimum ice thermal conductivity J g-1 K-1			
		use sia2_constants, only : &
			kelvin0, 							&	!	0degC in K			
			IceD,									&	! density of pure ice (g m-3)
			c_5,									&	! constant = 0.5
			c1e6,									&	! constant = 1e6
			c_1e6 									! constant = 1e-6

		implicit none

	! Subroutine Arguments
	! --------------------------------------------------------------------
		real(kind=dbl_kind), intent(in) :: &
			ice_t(z_max_ice),					& ! ice temperature (deg C)
			ice_s(z_max_ice),					& ! ice bulk salinity (psu)
			ice_d(z_max_ice),					& ! ice density (g m-3)
			ice_th(z_max_ice),				&	! ice layer thickness (m)
			snow_t(z_max_snow),				&	! snow density (g m-3)  NOTE - topmost snow layer is lowest in vector (upside down!!)
			snow_d(z_max_snow),				&	! snow density (g m-3)  NOTE - topmost snow layer is lowest in vector (upside down!!)
			snow_th(z_max_snow)					! snow layer thickness
		integer(kind=int_kind), intent(in) :: &
			ice_z,										& ! current number of ice layers
			snow_z								      ! current number of snow layers
		real(kind=dbl_kind), intent(out) :: &
			ki(z_max_pack+1),					& ! thermal conductivity profile (W m-1 K-1) 
			dth(z_max_pack+1)					  ! midpoint-midpoint layer thicknesses

	! Internal Variables
	! --------------------------------------------------------------------
		integer(kind=int_kind) :: &
			z,									&	! ice_z + snow_z
			ice_z1,							&	! 1st ice layer
			mo,									&	! switch for solution of surface temperature
			info,								&	! return value from LAPACK solver function
			ii,									&	! iterator
			jj										! iterator

		real(kind=dbl_kind) :: &
			rs,      						&	! rvec surface value
			s_mean,							&
			d_mean,							&
			ki1,								&
			ki_up,							&
			ki_down
				
	! Subroutine Code
	! --------------------------------------------------------------------
			! find internal indices
			ice_z1 = snow_z + 1
			z = ice_z + snow_z

			! top ice layer conductivity
			s_mean = ice_s(1) !(ice(sc,ic,mi)%s(1)+s_new(1))*c_5
			d_mean = ice_d(1) !(ice(sc,ic,mi)%d(1)+d_new(1))*c_5
			ki1 = sia2_ice_ki(ice_t(1),s_mean,d_mean)
			ki1 = max(ki1,ki_min)
	
		  if (snow_z .eq. 0) then

			  ! find intra-layer thickness
			  dth(1) = ice_th(1)*c_5	

        ! top layer conductivity
 			  ki(1) = ki1
			  ki_down = ki1

	    else			  			  

			  ! find intra-layer thickness
			  dth(1) = snow_th(snow_z)*c_5 
			  dth(ice_z1) = (snow_th(1) + ice_th(1))*c_5
              
              ! top layer conductivity
			  d_mean = c_1e6*snow_d(snow_z)
			  ki_down = sia2_snow_ki(d_mean)  ! sturm 1997 + 0.1495
              !ki_down = ksnow
			  ki(1) = ki_down
	
	          ! remaning snow intra-layer conductivities and thicknesses
	      if (snow_z .gt. 1) then
				  jj = 2			  
				  do ii=snow_z-1,1,-1
				      ! snow conductivity
				      ki_up = ki_down
                      d_mean = c_1e6*snow_d(ii)
                      ki_down = sia2_snow_ki(d_mean)  ! sturm 1997 + 0.1495				      
                      !ki_down = ksnow
                      ! equate fluxes to find intra-layer ki
                      ki(jj) = ki_up*ki_down/ &
                        (snow_th(ii+1)*ki_down + &
                        snow_th(ii)*ki_up)* &
                        (snow_th(ii+1)+snow_th(ii))
                      ! intra-layer snow thickness
					  dth(jj) = (snow_th(ii) + snow_th(ii+1))*c_5
					  jj = jj+1
				  enddo
			  endif	

			  ! equate fluxes to find intra-ice/snow ki
			  ki(ice_z1) =  ki1*ki_down/(snow_th(1)*ki1 + &
				ice_th(1)*ki_down)*(snow_th(1)+ice_th(1))

		  endif
	
		  ! find inter-layer ice thicknesses, conductivities
		  do ii=2,ice_z

              ! mid-point to mid-point layer thickness
			  dth(ii+snow_z) = (ice_th(ii-1) + ice_th(ii))*c_5

              ! intra-layer conductivity
			  ki_up = ki_down
			  s_mean = ice_s(ii) !(ice(sc,ic,mi)%s(ii)+s_new(ii))*c_5
			  d_mean = ice_d(ii) !(ice(sc,ic,mi)%d(ii)+d_new(ii))*c_5
			  ki_down = sia2_ice_ki(ice_t(ii),s_mean,d_mean)
	      ki_down = max(ki_down,ki_min)

			  ! equate fluxes to find inter-layer ki
			  ki(ii+snow_z) = ki_up*ki_down/(ice_th(ii-1)*ki_down + &
				ice_th(ii)*ki_up)*(ice_th(ii-1)+ice_th(ii))

		  enddo
	
		  ! last thickness, conductivity to water
		  dth(z+1) = ice_th(ice_z)*c_5
		  ki(z+1) = ki_down


end subroutine sia2_conductivity

! **********************************************************************


! **********************************************************************
! function sia2_snow_ki 
! Purpose: find snow thermal conductivity (W m-1 K-1)
! ----------------------------------------------------------------------

	real(kind=dbl_kind) function sia2_snow_ki(d_mean)

		implicit none

		! Function Arguments
		! --------------------------------------------------------------------
 		real(kind=dbl_kind), intent(in) :: d_mean   ! snow density

!  	sia2_snow_ki = 0.138_dbl_kind - 1.01_dbl_kind*d_mean + &
!    	3.233_dbl_kind*d_mean**2 + 0.1495_dbl_kind  ! sturm 1997 + 0.1495
    if (d_mean .ge. 0.4_dbl_kind) then
    	sia2_snow_ki = 0.5_dbl_kind
    else
    	sia2_snow_ki = 0.33_dbl_kind
		endif

	end function sia2_snow_ki

! **********************************************************************

! **********************************************************************
! function sia2_ice_ki 
! Purpose: find ice thermal conductivity (W m-1 K-1)
! ----------------------------------------------------------------------

	real(kind=dbl_kind) function sia2_ice_ki(t_mean,s_mean,d_mean)

		use sia2_constants, only : &
			c1e6,									&	! constant = 1e6
			IceD										! density of pure ice (g m-3)

		implicit none

		! Function Arguments
		! --------------------------------------------------------------------
 		real(kind=dbl_kind), intent(in) :: &
 			t_mean,									&   ! ice mean temperature (degC)
 			s_mean,									&   ! ice bulk salinity (psu)
 			d_mean										  ! ice density (g m-3)

  	sia2_ice_ki = (d_mean/IceD)*(2.11_dbl_kind-0.011_dbl_kind*t_mean + &
				0.09_dbl_kind*s_mean/t_mean - (d_mean-IceD)/c1e6)  ! pringle et al. 2006?				      

	end function sia2_ice_ki

! **********************************************************************


! **********************************************************************
! subroutine sia2_heat_solver 
! Purpose: calculates new temperature forcings given interface heat fluxes
! and ice column parameters
! ----------------------------------------------------------------------

	subroutine sia2_heat_solver(ts_next,t_next,ts_last,t_last,ki,dth, &
		ice_s,ice_d,ice_bd,ice_bv,ice_th,snow_d,snow_th, &
		z_ice,z_snow,ed_w_ice,ed_w_snow,dcheat,F0,dF0,T_bot,T_air, &
		T_diff,T_step,Fc_top,Fc_bot,dtt_s)

	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only : &
			ci0,									&	! heat capacity of pure ice (J g-1 K-1)
			cw, 									&	! heat capacity of cold seawater (J g-1 K-1)
			Lf,										&	! latent heat of fusion of water (J g-1)
			c0,										&	! zero
			c_5,									&	! constant = 0.5
			c_001,   							&	! constant = 0.001
			IceD,									&	! density of pure ice (g m-3)
			mu,										&	! linear liquidus constant for seawater
			kelvin0 								!	0degC in K
		use sia2_parameters, only : &
			z_max_ice,						&	! maximum ice layers
			z_max_snow,						&	! maximum snow layers
			z_max_pack,						&	! maximum ice+snow layers
			snow_model,						&	! switch for snow treatment
			ts_is_at,							&	! if non-zero, surface temperature is fixed at air temperature
			bb_f										! constant bubble fraction in sea ice (fraction)
			
		implicit none
		
	! Function Arguments
	! --------------------------------------------------------------------
			real(kind=dbl_kind), intent(inout) :: &
				ts_last,						&	! previous surface temperature (degC)
				ts_next,						&	! calculated surface temp (degC) - must be estimated before calling function
				t_next(z_max_pack)					! calculated snow/ice temperature profile - (degC) must be estimated before calling function
			real(kind=dbl_kind), intent(in) :: &
				t_last(z_max_pack),			&	! previous snow/ice temperature profile
				ki(z_max_pack+1),				&	! thermal conductivity profile (W m-1 K-1) 
				dth(z_max_pack+1),				&	! midpoint-midpoint layer thicknesses
				ed_w_ice(z_max_ice),				&	!	ice shortwave irradiance absorption profile (W m-2)
				ed_w_snow(z_max_snow),			&	!	snow shortwave irradiance absorption profile (W m-2) 
				ice_s(z_max_ice),				&	! ice bulk salinity (psu)
				ice_d(z_max_ice),				&	! ice density (g m-3)
				ice_bd(z_max_ice),			&	! ice brine density (g m-3)
				ice_bv(z_max_ice),			&	! ice density (g m-3)
				ice_th(z_max_ice),			&	! ice layer thickness (m)
				dcheat(z_max_ice),			&	! additional heat source, typically due to desalination/latent heat (W m-2)
				snow_d(z_max_snow),				&	! snow density (g m-3)  NOTE - topmost snow layer is lowest in vector (upside down!!)
				snow_th(z_max_snow),				&	! snow layer thickness
				F0,									&	! surface heat balance, less conductive flux (W m-2)
				dF0,								&	! derivative of F0 with respect to surface temperature (W m-2 K-1)
				T_bot,							&	! bottom ice temperature (fixed, degC)
				T_air,							&	! surface air temperature (K)
				dtt_s									! time step length (s)
			integer(kind=int_kind), intent(in) :: &
				z_ice,							&	! current number of ice layers
				z_snow								! current number of snow layers
			real(kind=dbl_kind), intent(out) :: &
				T_diff,							&	! max difference in layer temperature of solution
				T_step,							&	! max difference in layer temperature of convergence
				Fc_top,							&	! surface conductive heat flux (W m-2)
				Fc_bot								! basal conductive heat flux (W m-2)

	! Internal Variables
	! --------------------------------------------------------------------
			integer(kind=int_kind) :: &
				z,									&	! z_ice + z_snow
				z1,									&	! z+mo
				mo,									&	! switch for solution of surface temperature
				info,								&	! return value from LAPACK solver function
				ii,									&	! iterator
				jj										! iterator
	
			real(kind=dbl_kind) :: &
				rs,      						&	! rvec surface value
				r1,									&	! rvec 1st layer value
				xl,									&	! tridiagonal calc intermediate
				Fm0,								&	! surface flux balance (W m-2)
				bv_mean,						&	! mean brine volume (ppt)
				s_mean,							&	! mean bulk salinity (psu)
				d_mean,							&	! mean density (g m-3)
				t_melt,							&	! melt temperature (K)
				mat_a,							&	! tridiagonal term temporary
				mat_b,							&	! tridiagonal term temporary
				mat_c,							&	! tridiagonal term temporary
				mat_sc,							&	! tridiagonal term temporary
				mat_sb,							&	! tridiagonal term temporary
				tmp1,								&	! temporary var
				tmp2,								&	! temporary var
				tmp3,								&	! temporary var
				ci(z_max_pack+1),				&	! layer heat capacity * density (J K-1 m-3)
				r_vec(z_max_pack+2),			&	! right-hand-side (rhs) vector
				DU(z_max_pack+1),				&	! upper tridiagonal vector
				DC(z_max_pack+2),				&	! central tridiagonal vector
				DL(z_max_pack+1),				&	! lower tridiagonal vector				
				DU_calc(z_max_pack+1),		&	! upper tridiagonal vector
				DC_calc(z_max_pack+2),		&	! central tridiagonal vector
				DL_calc(z_max_pack+1)			! lower tridiagonal vector				
		
			! total snow/ice column layers
			z = z_ice + z_snow

      ! surface matrix constant - same for all cases
      rs = -1.0_dbl_kind*F0 + dF0*Ts_last		

      ! top layer is SNOW
      ! ----------------------------------------------------------------
      if (z_snow .gt. 0) then              

          ! mean temperature
          tmp1 = (T_last(1)+T_next(1))*c_5
          tmp1 = min(tmp1,c0)
          
          ! snow heat capacity averaged between T_last and T_next, * snow density
          tmp2 = 0.185_dbl_kind + 0.689e-2_dbl_kind*(tmp1+kelvin0)
          tmp2 = max(0.02_dbl_kind,tmp2)
          ci(1) = snow_d(z_snow)*tmp2
          xl = 2.*dtt_s/((dth(1)+dth(2))*ci(1))

          ! assign constant rhs vector 
          tmp3 = dtt_s/(ci(1)*snow_th(z_snow))
          r1 = T_last(1) + Ed_W_snow(1)*tmp3              

          ! determine of snow is at the melting point.  If so, toggle mo variable
          ! for different surface heat balance equation
          Fc_top = ki(1)*(Ts_last-T_last(1))/dth(1)
          Fm0 = F0-Fc_top
          if (snow_model .eq. 1 .or. snow_model .eq. 3) then
              if ((Ts_last .ge. c0) .and. (Fm0 .gt. c0)) then
                  ! melting - no need to calc new surface temp
                  mo = 0
                  Ts_next = c0
              else
                  mo = 1
              endif
          else
              mo = 1
          endif

      ! top layer is ICE
      ! ----------------------------------------------------------------
      else

!          d_mean = (d_new(1) + ice(sc,ic,mi)%d(1))*c_5
!          s_mean = (ice(sc,ic,mi)%s(1) + s_new(1))*c_5
					d_mean = ice_d(1)
          s_mean = ice_s(1)
          bv_mean = ice_bv(1)*c_001          
          t_melt = s_mean*mu

          tmp1 = min(T_last(1),T_melt)
          tmp2 = min(T_next(1),T_melt)

          ci(1) = (1.0_dbl_kind-bb_f) * ( &         ! remove air fraction
              ice_bd(1)*bv_mean*cw +   &    ! heat capacity of pure ice
              IceD*(1.0_dbl_kind-bv_mean)*ci0 - &    ! heat capacity of brine 
              IceD*Lf*T_melt/(tmp1*tmp2))   ! latent heat of freezing based on temp change

          xl = 2.*dtt_s/((dth(1)+dth(2))*ci(1))
          ! assign constant vector - include irradiance absorbed in ice
          tmp3 = dtt_s/(ci(1)*ice_th(1))
          r1 = T_last(1) + tmp3*Ed_W_ice(1) + xl*dcheat(1)      

          ! look at Fc beforehand to see whether melting is occuring
          Fc_top = ki(1)*(Ts_last-T_last(1))/dth(1)
          Fm0 = F0-Fc_top
					if ((Ts_last .ge. c0) .and. (Fm0 .gt. c0)) then
							! melting - no need to calc new surface temp
							mo = 0
							Ts_next = c0
					else
							mo = 1
					endif
 
      endif             

      if (ts_is_at .gt. 0) then
          Ts_last = T_air - kelvin0
          Ts_next = Ts_last
          mo = 0
      endif
      
      if (mo .eq. 1) then

          ! surface matrix coeffs. 
          mat_sc = ki(1)/dth(1)                   
          mat_sb = dF0 - ki(1)/dth(1) 

          mat_a = -1.*xl*ki(1)/dth(1)   
          mat_c = -1.*xl*ki(2)/dth(2)
          mat_b = 1. - mat_a - mat_c
        
          ! assign constant vector
          r_vec(1) = rs
          r_vec(2) = r1
         
          ! assign top row of matrix
          DU(1) = mat_sc
          DC(1) = mat_sb

          ! assign second row of matrix
          DU(2) = mat_c
          DC(2) = mat_b
          DL(1) = mat_a

     else             
          ! assign constant vector
          r_vec(1) = r1 + Ts_last*xl*ki(1)/dth(1) ! <-- not needed when Ts_last = 0
         
          ! assign top row of matrix
          mat_a = -1.*xl*ki(1)/dth(1)   
          mat_c = -1.*xl*ki(2)/dth(2)
          mat_b = 1. - mat_a - mat_c

          DU(1) = mat_c
          DC(1) = mat_b

     endif

     do ii=2,z

         ! COMPUTE INTERNAL LAYER MATRIX COEFFICIENTS
         ! ----------------------------------------------------             
         if (ii .lt. z_snow+1) then  
             tmp1 = (T_last(ii)+T_next(ii))*c_5
             tmp1 = min(tmp1,c0)

             ! snow heat capacity averaged between T_last and T_next, *ice_d
             tmp2 = 0.185_dbl_kind + 0.689e-2_dbl_kind*(tmp1+kelvin0)
             tmp2 = max(0.02_dbl_kind,tmp2)
             ci(ii) = snow_d(z_snow-ii+1)*tmp2
             xl = 2.*dtt_s/((dth(ii)+dth(ii+1))*ci(ii))

             ! assign constant vector w/ absorbed irradiance
             !r_vec(ii+mo) = T_last(ii) + xl*Ed_W_snow(ii)
             tmp3 = dtt_s/(ci(ii)*snow_th(z_snow-ii+1))
             r_vec(ii+mo) = T_last(ii) + tmp3*Ed_W_snow(ii)

         else

              ! heat flux dependent desalination
              ! ------------------------------------------------------------
              jj = ii-z_snow
    !          d_mean = (d_new(jj) + ice_d(jj))*c_5
    !          s_mean = (ice(sc,ic,mi)%s(jj) + s_new(jj))*c_5
    					d_mean = ice_d(jj)
              s_mean = ice_s(jj)
              t_melt = s_mean*mu
              bv_mean = ice_bv(jj)*c_001
    
              tmp1 = min(T_last(ii),T_melt)
              tmp2 = min(T_next(ii),T_melt)

              ci(ii) = (1.0_dbl_kind-bb_f) * ( &            ! remove air fraction
                  ice_bd(jj)*bv_mean*cw +  &    ! heat capacity of pure ice
                  IceD*(1.0_dbl_kind-bv_mean)*ci0 - &    ! heat capacity of brine 
                  IceD*Lf*T_melt/(tmp1*tmp2))   ! latent heat of freezing based on temp change

              xl = 2.0_dbl_kind*dtt_s/((dth(ii)+dth(ii+1))*ci(ii))
              tmp3 = dtt_s/(ci(ii)*ice_th(jj))

              ! assign constant vector - include irradiance absorbed in ice
               r_vec(ii+mo) = T_last(ii) &
!                + tmp3*(dcheat(jj) + Ed_W_temp(jj))
                + tmp3*Ed_W_ice(jj) + xl*dcheat(jj)

         endif

         ! assign row of matrix
         DU(ii+mo) = -1.*xl*ki(ii+1)/dth(ii+1) ! replaced above to reduce computation
         DL(ii-1+mo) = -1.*xl*ki(ii)/dth(ii) 
         DC(ii+mo) = 1. - DU(ii+mo) - DL(ii-1+mo)   

         ! ammend matrix values for final ice layer, when bottom temp is fixed
         if (ii .eq. z) then

             r_vec(z+mo) = r_vec(z+mo) - DU(ii+mo)*T_bot

             ! add heat absorbed in skeletal layers to bottom layer 
             ! to maintain conservation of heat
!             do jj=sk_1,sk_z
!                 r_vec(z+mo) = r_vec(z+mo) + tmp3*(Ed_W_ice(jj))
!             enddo
 
          endif
                       
     enddo                 


     ! SOLVE MATRIX PROBLEM FOR IMPLICIT HEAT CONDUCTION
     ! ---------------------------------------------------------
     z1 = z+mo
     DC_calc = DC
     DU_calc = DU
     DL_calc = DL      

!             print *,'Iteration: ',jjj,z_snow,sk_z,z
!             print *,'heat_stop: ',heat_stop,int_z
!             print *,'T_last: ',t_last
!             print *,'T_next: ',t_next
!             print *,'DU: ',du
!             print *,'DC: ',dc
!             print *,'DL: ',dl
!             print *,'r_vec: ',r_vec

		if (DBL_KIND .eq. 8) then
	    call DGTSV(z1,1,DL_calc,DC_calc,DU_calc,r_vec,z1,info)
    else
	    call SGTSV(z1,1,DL_calc,DC_calc,DU_calc,r_vec,z1,info)
		endif
 !            print *,'info: ',info
 !            print *,'r_vec: ',r_vec
     
     ! triage new temps
     if (mo .eq. 1) then
         T_diff = abs(r_vec(1) - Ts_next)
         T_step = abs(r_vec(1) - Ts_last)
         Ts_next = r_vec(1)
     else
         T_diff = c0
         T_step = c0
         !Ts_next = 0.
     endif
     do ii=1,z
         T_diff=max(T_diff,abs(r_vec(ii+mo)-T_next(ii)))
         T_step=max(T_step,abs(r_vec(ii+mo)-T_last(ii)))
         T_next(ii) = r_vec(ii+mo)
     enddo

    ! conductive heat flux at boundaries
    Fc_bot = ki(z+1)*(T_bot-T_next(z))/dth(z+1) !  flux over sub-dt timestep
    Fc_top = ki(1)*(Ts_next-T_next(1))/dth(1)  !  flux over sub-dt timestep

	end subroutine sia2_heat_solver

! **********************************************************************

! **********************************************************************
! subroutine sia2_F0 
! Purpose: find balance (excepting conductive heat) of surface heat 
! flux (F0) and derivative (dF0) 
! ----------------------------------------------------------------------

	!pure subroutine sia2_F0(Ts,pc,F0,dF0,Fe)
	subroutine sia2_F0(Ts,pc,F0,dF0,Fe,F_s,F_little_L)

	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only : &
			c1e5,										&	! 100,000
			steph_boltz, 						&	! stephan boltzmann constant
			qq1,										&	! latent heat constant 1
			qq2,										&	! latent heat constant 2
			c4,											&	! 4
			c3,											&	! 3
			c2,											&	! 2
			c1,											&	! 1
			nc1,										&	! -1
			fe_A, 									&	! latent turbulent heat eq. constant
			fe_B, 									&	! latent turbulent heat eq. constant
			fe_C, 									&	! latent turbulent heat eq. constant
			fe_D, 									&	! latent turbulent heat eq. constant
			fe_E, 									&	! latent turbulent heat eq. constant
			kelvin0 									!	0degC in K
		use sia2_parameters, only : &
			ch,											&	! 
			Cp,											&	! 
			atmo											! atmospheric/ice heat flux switch (1 = CICE, 0=kottmeier et al., 2003)
			
		implicit none
		
	! Function Arguments
	! --------------------------------------------------------------------
		real(kind=dbl_kind), intent(in) :: &
			Ts											! surface temperature (Kelvin)
		type (heat_pre_calc), intent(in) :: &
			pc											! pre-calculated surface heat flux parameters
		real(kind=dbl_kind), intent(out) :: &
			F0,										&	! surface heat balance (W m-2)
			dF0, 									&	! surface heat balance derivative (W m-2 K-1)
			F_little_L,						&	! outgoing longwave irradiance (W m-2)
			F_s,									&	! sensible turbulent flux (W m-2)
			Fe 											! latent turbulent flux (W m-2)

	! Internal Variables
	! --------------------------------------------------------------------
		real(kind=dbl_kind) :: &
			shcoef,								&	! sensible turbulent heat flux coefficient
			lhcoef,								&	! latent turbulent heat flux coefficient
			F0_constants,					&	! F0 constant pasts, Kottmeier et al. 2003
			F0_4,									&	! F0 polynomial coefficient
			F0_3,									&	! F0 polynomial coefficient
			F0_2,									&	! F0 polynomial coefficient
			F0_1										! F0 polynomial coefficient

		
		if (atmo .eq. 1) then
	
			call sia2_env_atmo(Ts-kelvin0,pc%potT,pc%v10,10.,pc%shum, &
				pc%rho_air,lhcoef,shcoef)
			!call sia2_env_atmo_const(Ts-kelvin0,pc%v10,pc%rho_air,lhcoef,shcoef)

			F_little_L = nc1*pc%eps0*steph_boltz*Ts**4 ! outgoing longwave - W/m^2/K^4 * K^4 = W/m^2		
			F_s = shcoef*(pc%potT - Ts)		! sensible heat flux
			Fe = lhcoef*(pc%shum - (qq1/pc%rho_air)*exp(-qq2/Ts))	! latent heat flux/sublimation	
	
			F0 = pc%Fr + pc%FL + F_little_l + F_s + Fe  
			dF0 = nc1*( shcoef +  &  ! term from sensible turbulent flux (F_s)
						c4*pc%eps0*steph_boltz*Ts**3 +  & ! term from outgoing longwave (F_little_l)
						(lhcoef*qq1*qq2)/(pc%rho_air*exp(qq2/Ts)*Ts**2) &  ! term from turbulent latent flux (Fe)
						)

			print *,'Fr FL Fl F_s Fe: '
			print *,pc%Fr,pc%Fl,F_little_l,F_s,Fe
			print *,'Ts F0: ',Ts,F0
		else
		
			! total surface heat flux minus conductive part
			F0 = pc%F0_constants + pc%F0_4*Ts**4 + pc%F0_3*Ts**3 + &
				pc%F0_2*Ts**2 + pc%F0_1*Ts

			! slope of above equation
			dF0 = c4*pc%F0_4*(Ts**3) + c4*pc%F0_3*(Ts**2) + &
				c4*pc%F0_2*Ts + pc%F0_1            
    
			! debug vars ---------------
!			Fm0 = Fr + FL + Fe_sub2 - (eps0*steph_boltz+Fe_sub1*fe_A)*ts_kelvin**4 - &
!							Fe_sub1*fe_B*ts_kelvin**3 - Fe_sub1*fe_C*ts_kelvin**2 - (Fe_sub1*fe_D + &
!							rho_air*cp*Ch*v10)*ts_kelvin + Fs_sub 
!
			F_little_L = -1*pc%eps0*steph_boltz*Ts**4 ! W/m^2/K^4 * K^4 = W/m^2


			F_s = pc%rho_air*cp*Ch*pc%v10*(pc%Tair-Ts)  ! kg/m^3 * J/kg/K * (dimensionless) * m/s * K = J/s/m^2 = W/m^2
			! debug vars ---------------                 

			Fe = pc%Fe_sub1*(fe_A*(pc%rhum* &
					pc%Tair**4-Ts**4)+fe_B*(pc%rhum* &
					pc%Tair**3-Ts**3)+fe_C*(pc%rhum* &
					pc%Tair**2-Ts**2)+fe_D*(pc%rhum* &
					pc%Tair-Ts)+fe_E*(pc%rhum-c1)) ! W/m^2

		endif		

	end subroutine sia2_F0

! **********************************************************************



! **********************************************************************
! subroutine sia2_F0 
! Purpose: find balance (excepting conductive heat) of surface heat 
! flux (F0) and derivative (dF0) 
! ----------------------------------------------------------------------

	pure subroutine sia2_heat_pre_calc(Tair,v10,pres,shum,rhum,dewpoint, &
		sn_depth,fc,Fr,Fl,pc)

	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only : &
			c1e5,										&	! 100,000
			steph_boltz, 						&	! stephan boltzmann constant
			R_air,									&	! latent heat constant 1
			R_h2o,									&	! latent heat constant 2
			Lv,											&	! latent heat of vaporization of water (J g)
			nc1,										&	! -1
			c0,											&	! 0
			c_01,										&	! 0.01
			c1,											&	! 1
			c100,										&	! 100
			c1e3, 									&	! 1000
			fe_A, 									&	! latent turbulent heat eq. constant
			fe_B, 									&	! latent turbulent heat eq. constant
			fe_C, 									&	! latent turbulent heat eq. constant
			fe_D, 									&	! latent turbulent heat eq. constant
			fe_E 											! latent turbulent heat eq. constant
		use sia2_parameters, only : &
			Ce,											&	! 
			ch,											&	! 
			Cp,											&	! 
			h_snowpatch,						& ! 50% coverage snow depth (m) 
			eps_ice,								&	! emissivity of ice
			eps_snow									! emissivity of snow
			
		implicit none
		
	! Function Arguments
	! --------------------------------------------------------------------
		real(kind=dbl_kind), intent(in) :: &
			Tair,									&	! air temperature (Kelvin)
			v10,									&	! 10m wind speed (m/s)
			pres,									&	! surface air pressure (mbar)
			shum,									&	! specific humidity (kg kg-1)
			rhum,									&	! relative humidity (fraction)
			dewpoint,							&	! (Kelvin)
			sn_depth,							&	! snow depth (m)
			fc,										&	! cloud fraction
			Fr,										&	! downward shortwave irradiance
			Fl											! downward longwave irradiance
		type (heat_pre_calc), intent(out) :: &
			pc											! pre-calculated surface heat flux parameters

	! Internal Variables
	! --------------------------------------------------------------------
		real(kind=dbl_kind) :: &
			p_h2o_sat,						&	! saturation water pressure of air (mb)
			p_h2o,								&	! partial pressure of water vaptor (mb)
			shum_sat,								&	! saturation specific humidity (kg kg-1)
			epsa										! emissivity of air

		! these are required for calcs below
		pc%Tair = Tair
		pc%v10 = v10
		pc%pres = pres

		! water vapor saturation pressure
		p_h2o_sat = p_sat(Tair) ! mb (1 mb = 100 Pa = 1 hPa)

		! find humidities and rho_air
		shum_sat = qsat(1,Tair,pres) ! kg/kg
		if (shum .lt. c0) then
			! use dewpoint to find humidities
			if (rhum .lt. c0) then				
				!p_h2o = p_sat(dewpoint)
				!pc%rhum = max(c1,p_h2o/p_h2o_sat)
				pc%shum = qsat(1,dewpoint,pres) !kg/kg
				pc%rhum = pc%shum/shum_sat
			else
				pc%shum = shum_sat*rhum				
				pc%rhum = rhum
			endif		
			!p_h2o = p_h2o_sat*pc%rhum
			!pc%shum = 0.622_dbl_kind*p_h2o/(pres*c_01 - 1.622_dbl_kind*p_h2o)
			pc%rho_air = rhoair(Tair,pres*c100,pc%shum)	 ! covert pres to Pa
		else	
			pc%shum = shum
			pc%rho_air = rhoair(Tair,pres*c100,pc%shum)	
      !pc%rhum = pc%rho_air*pc%shum*R_h2o*Tair / (p_h2o_sat*c100)  
      pc%rhum = pc%shum/shum_sat
    endif
    
    ! potential temperature
    pc%potT = Tair*(c1e5/(pres*c100))**(R_air/cp) ! T(kelvin)*(P0/P)^(R/cp)

		! estimate fractional snow coverage (from CICE)
		if (sn_depth .gt. c0) then
			if (h_snowpatch .eq. c0) then
				pc%f_snow = c1
			else
				pc%f_snow = sn_depth/(sn_depth + h_snowpatch)			  
			endif
		else
			pc%f_snow = c0
		endif

		! pass along shortwave irradiance, if any calculated here
		if (Fr .gt. c0) then
			pc%Fr = Fr
		else
			pc%Fr = c0
		endif

		! incoming longwave irradiance
		if (Fr .gt. c0) then
			! find longwave emissivity basced on atmospheric cloud cover;  (König-Langlo & Augstein 1994)
			epsa = 0.765_dbl_kind + 0.22_dbl_kind*(fc*c_01)**3 ! longwave emissivity of air
			pc%Fl = epsa*steph_boltz*Tair**4
    else
      pc%Fl = Fl
		endif
		
		! estimate ice pack emissivity
		pc%eps0 = eps_snow*pc%f_snow + eps_ice*(c1 - pc%f_snow) 

		! find first static part of turbulent heat flux equation
		pc%Fe_sub1 = 0.622_dbl_kind*(Lv*c1e3)*pc%rho_air*Ce*v10/(pres)  ! J/kg/K * kg/m^3 * (dimensionless) * m/s * 1/mbar = J/m^2/mbar/K

		! find second static part of turbulent heat flux equation
		pc%Fe_sub2 = pc%Fe_sub1*(pc%rhum*(fe_A*Tair**4 &
			+ fe_B*Tair**3 &
			+ fe_C*Tair**2 &
			+ fe_D*Tair + fe_E) - fe_E)

		! static part of sensible latent heat flux (Fs)
		pc%Fs_sub = pc%rho_air*cp*Ch*v10*Tair	

		! update pre-calc kottmeier et al. 2003 polynomial multipliers
		pc%F0_constants = pc%Fr + pc%FL + pc%Fe_sub2 + pc%Fs_sub
		pc%F0_4 = nc1*(pc%eps0*steph_boltz+pc%Fe_sub1*fe_A)
		pc%F0_3 = nc1*pc%Fe_sub1*fe_B
		pc%F0_2 = nc1*pc%Fe_sub1*fe_C
		pc%F0_1 = nc1*(pc%Fe_sub1*fe_D + pc%rho_air*cp*Ch*v10)



	end subroutine sia2_heat_pre_calc

! **********************************************************************


! **********************************************************************
! function p_sat 
! Purpose: saturation vapor pressure at temp  
! ----------------------------------------------------------------------
	real(kind=dbl_kind) pure function p_sat(T) 
		
		implicit none
		
		real(kind=dbl_kind), intent(in) :: &
			T							! Air temperature (K)
			
		p_sat = 6.1078_dbl_kind*10**( &
			(7.5_dbl_kind*T-2048.625_dbl_kind) / (T-35.85_dbl_kind))  ! (mb)

	end function p_sat


! **********************************************************************

! **********************************************************************
! function qsat 
! Purpose: calculate saturation humidity in (kg/kg) for given
!     temperature in T (K), and pressure P(mb) over 
!			either water (mode = 0), or ice (mode = 1)
! ----------------------------------------------------------------------
	real(kind=dbl_kind) pure function qsat(mode, T, P)

		implicit none

		integer, intent(in) :: &
			mode							! 1=ice, 0 = water
		real(kind=dbl_kind), intent(in) :: &
			T,							& ! temperature (K)
			P								! mean sea level pressure (mb)
		real(kind=dbl_kind) :: &
			P_h2o							! partial pressure of h2o in air

    !if(mode == 1) THEN
    !    qsat  =  11637800./exp(5897.8/T) ! kg/m^3
    !else 
    !    qsat  =  640380./exp(5107.4/T)   ! kg/m^3
    !endif
    
    ! alternate qsat - Saenz 4/2012
    P_h2o = 6.1078*10**((7.5*T-2048.625)/(T-35.85))  ! vapor pressure (mb)
    qsat = 0.622*P_h2o/(P - 1.622*P_h2o)  ! kg/kg

  end FUNCTION qsat

! **********************************************************************



! **********************************************************************
! function rhoair 
! Purpose: air density (kg m-3) at temperature T, pressure P, 
! specific humidity H 
! ----------------------------------------------------------------------
	real(kind=dbl_kind) pure function rhoair(T,P,H) 
		
		use sia2_constants, only : &
			R_air,				& ! universal gas constant for air
			R_h2o					! universal gas constant for water vapor
		
		implicit none

		real(kind=dbl_kind), intent(in) :: &
			T,						&	! Air temperature (K)
			P,						&	! Air pressure (Pa)
			H								! Specific humidity (kg kg-1)
			
			rhoair = P / (T*(R_air - H*R_air + H*R_h2o))

	end function rhoair
! **********************************************************************


! **********************************************************************
! function ice_albedo 
! Purpose: find ice albedo according to CICEv3 (originally tuned using SHEBA data) 
! ----------------------------------------------------------------------
	real(kind=dbl_kind) pure function ice_albedo(Ts,ice_depth,f_snow) 
		
		use sia2_constants, only : &
			nc1,				& ! -1
			c_5,				& ! 0.5
			c1						! 1
		use sia2_parameters, only : &
			alb_i_wet,				& ! wet ice albedo
			alb_i_dry,				& ! dry ice albedo
			alb_s_wet,				& ! wet snow albedo
			alb_s_dry,				&	! dry snow albedo
			atan_c_i						! pre-calculated inverse tangent
		
		implicit none

	! Function Arguments
	! --------------------------------------------------------------------
		real(kind=dbl_kind), intent(in) :: &
			Ts,						&	! Ice/snow surface temperature (mb)
			ice_depth,		&	! sea ice thickness
			f_snow					! sea ice thickness
			
	! Internal Variables
	! --------------------------------------------------------------------
		real(kind=dbl_kind) :: &
			albedo_i,			&	! intermediate ice albedo
			albedo_s				! intermediate snow albedo

	! Function Code
	! --------------------------------------------------------------------

		! find albedo(s)
		if (Ts .gt. nc1) then
			! scale albedo by difference between wet and dry ice between 0 and -1 degC
			albedo_i = alb_i_wet + (alb_i_dry - alb_i_wet)*min(c1,abs(Ts))
			albedo_s = alb_s_wet + (alb_s_dry - alb_s_wet)*min(c1,abs(Ts))
		else
			albedo_i = alb_i_dry
			albedo_s = alb_s_dry
		endif  
		if (ice_depth .lt. c_5) then
			! scale using atan function toward open water albedo as thickness approaches 0
			albedo_i = albedo_i - (albedo_i-0.06_dbl_kind)* &        ! 0.06 = open water albedo
				(c1-atan(atan_c_i*ice_depth))
		endif	
		ice_albedo = f_snow*albedo_s + (c1-f_snow)*albedo_i

	end function ice_albedo
! **********************************************************************



! **********************************************************************
!=======================================================================
!BOP
!
! %DESCRIPTION:
!
! Atmospheric boundary interface (stability based flux calculations)
!
! !REVISION HISTORY:
!  SVN:$Id: ice_atmo.F90 140 2008-07-25 20:15:53Z eclare $
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!=======================================================================
!BOP
!
! !IROUTINE: atmo_boundary_layer - compute coefficients for atm-ice fluxes, 
!                                  stress and Tref/Qref
!
! !INTERFACE:
!
	pure subroutine sia2_env_atmo(Tsf,potT,wind,zlvl,Qa,rhoa,lhcoef,shcoef)
 
		implicit none

! %DESCRIPTION:
!
! Compute coefficients for atm/ice fluxes, stress, and reference
! temperature and humidity. NOTE: \\
! (1) all fluxes are positive downward,  \\
! (2) here, tstar = (WT)/U*, and qstar = (WQ)/U*,  \\
! (3) wind speeds should all be above a minimum speed (eg. 1.0 m/s). \\
!
! ASSUME:
!  The saturation humidity of air at T(K): qsat(T)  (kg/m**3)
!
! Code originally based on CSM1
!
! %REVISION HISTORY: same as module
!
! %USES:
!
! %INPUT/OUTPUT PARAMETERS:
!

		 real(kind=dbl_kind), intent(in) :: &
				Tsf      , & ! surface temperature of ice or ocean (degC)
				potT     , & ! air potential temperature  (K)          
				wind     , & ! wind speed (m/s)
				zlvl     , & ! atm level height (m)
				Qa       , & ! specific humidity (kg/kg)
				rhoa         ! air density (kg/m^3)

		 real(kind=dbl_kind), intent(out) :: &
				shcoef   , & ! transfer coefficient for sensible heat
				lhcoef       ! transfer coefficient for latent heat
! !
! !EOP
! !
! 
		 integer :: &
				k

		 real(kind=dbl_kind) :: &
				TsfK  , & ! surface temperature in Kelvin (K)
				xqq   , & ! temporary variable
				psimh , & ! stability function at zlvl   (momentum)
				psimhs, & ! stable profile
				ssq   , & ! sat surface humidity     (kg/kg)
				qqq   , & ! for qsat, dqsfcdt
				TTT   , & ! for qsat, dqsfcdt
				qsat  , & ! the saturation humidity of air (kg/m^3)
				delt  , & ! potential T difference   (K)
				delq  , & ! humidity difference      (kg/kg)
				Lheat     ! Lvap or Lsub, depending on surface type

		 real(kind=dbl_kind) :: &
				ustar , & ! ustar (m/s)
				tstar , & ! tstar
				qstar , & ! qstar
				rdn   , & ! sqrt of neutral exchange coefficient (momentum)
				rhn   , & ! sqrt of neutral exchange coefficient (heat)
				ren   , & ! sqrt of neutral exchange coefficient (water)
				rd    , & ! sqrt of exchange coefficient (momentum)
				re    , & ! sqrt of exchange coefficient (water)
				rh    , & ! sqrt of exchange coefficient (heat)
				vmag  , & ! surface wind magnitude   (m/s)
				alz   , & ! ln(zlvl  /z10)
				thva  , & ! virtual temperature      (K)
				cp    , & ! specific heat of moist air
				hol   , & ! H (at zlvl  ) over L
				stable, & ! stability factor
				psixh     ! stability function at zlvl   (heat and water)

			 ! from CICE 4 ice_constants.F90
		real(kind=dbl_kind), parameter :: &
			 c0 = 0.e0, &
			 c1 = 1.e0, &
			 c2 = 2.e0, &
			 c8 = 8.e0, &
			 c10 = 10.e0, &
			 c16 = 16.e0, &
			 p5 = 0.5e0, &
			 iceruf   = 0.0005e0   ,&! ice surface roughness (m)
			 zref   = 10.e0   ,&! reference height for stability (m)
			 qqqice  = 11637800.e0   ,&! for qsat over ice
			 TTTice  = 5897.80e0      ,&! for qsat over ice
			 Lsub      = 2.835e6 ,&! latent heat, sublimation freshwater (J/kg)
			 Lvap      = 2.501e6 ,&! latent heat, vaporization freshwater (J/kg)
			 zvir      = 0.606e0   ,&! rh2o/rair - 1.0
			 R_air = 285.05e0, & ! J/kg/K
			 R_h2o = 461.495e0, &
			 vonkar    = 0.4e0     ,&! von Karman constant
			 cp_wv     = 1.81e3  ,&! specific heat of water vapor (J/kg/K)
			 cp_air    = 1005.0e0  ,&! specific heat of air (J/kg/K)
			 ! (Briegleb JGR 97 11475-11485  July 1992)
			 Tffresh   = 273.15e0  ,&! freezing temp of fresh ice (K)
			 gravit    = 9.80616e0 , &     
			 cpvir = cp_wv/cp_air-c1 , & ! defined as cp_wv/cp_air - 1.
			 umin  = c1, &  !             , &    ! minimum wind speed (m/s)
			 pih = 0.5*3.14159265358979323846e0 ! p5*pi
								
		 ! local functions
		 real(kind=dbl_kind) :: &
				xd    , & ! dummy argument
				psimhu, & ! unstable part of psimh
				psixhu    ! unstable part of psimx

		!------------------------------------------------------------
		! Define functions
		!------------------------------------------------------------

		 psimhu(xd)  = log((c1+xd*(c2+xd))*(c1+xd*xd)/c8) &
								 - c2*atan(xd) + pih
!ech                  - c2*atan(xd) + 1.571_dbl_kind

		 psixhu(xd)  =  c2 * log((c1 + xd*xd)/c2)

		!------------------------------------------------------------
		! Initialize
		!------------------------------------------------------------

!          print *,'------------------------------------------------------------'
!          print *,'Tsf: ',Tsf      ! surface temperature of ice or ocean (degC)
!          print *,'potT: ', potT    ! air potential temperature  (K)          
!          print *,'wind: ', wind    ! wind speed (m/s)
!          print *,'zlvl: ', zlvl   ! atm level height (m)
!          print *,'Qa: ', Qa      ! specific humidity (kg/kg)
!          print *,'rhoa: ', rhoa         ! air density (kg/m^3)



			 shcoef = c0
			 lhcoef = c0

		!------------------------------------------------------------
		! define some more needed variables
		!------------------------------------------------------------

			 qqq  = qqqice          ! for qsat
			 TTT  = TTTice          ! for qsat
			 if (Tsf .ge. c0) then
					 Lheat = Lvap           ! water to vapor
			 else
					 Lheat = Lsub           ! ice to vapor
			 endif
			 vmag = max(umin, wind)
			 rdn  = vonkar/log(zref/iceruf) ! neutral coefficient


			 TsfK       = Tsf + Tffresh     ! surface temp (K)
			 qsat       = qqq * exp(-TTT/TsfK)   ! saturation humidity (kg/m^3)
			 ssq        = qsat / rhoa       ! sat surf hum (kg/kg)

			 thva   = potT * (c1 + zvir * Qa) ! virtual pot temp (K)
			 delt  = potT - TsfK       ! pot temp diff (K)
			 delq  = Qa - ssq          ! spec hum dif (kg/kg)
			 alz    = log(zlvl/zref)
			 cp     = cp_air*(c1 + cpvir*ssq)

		!------------------------------------------------------------
		! first estimate of Z/L and ustar, tstar and qstar
		!------------------------------------------------------------

			 ! neutral coefficients, z/L = 0.0
			 rhn = rdn
			 ren = rdn

			 ! ustar,tstar,qstar
			 ustar = rdn * vmag
			 tstar = rhn * delt
			 qstar = ren * delq

		!------------------------------------------------------------
		! iterate to converge on Z/L, ustar, tstar and qstar
		!------------------------------------------------------------

		do k=1,5

			! compute stability & evaluate all stability functions
					hol = vonkar * gravit * zlvl &
								 * (tstar/thva &
									+ qstar/(c1/zvir+Qa)) &
								 / ustar**2

					hol    = sign( min(abs(hol),c10), hol )
					stable = p5 + sign(p5 , hol)
					xqq    = max(sqrt(abs(c1 - c16*hol)) , c1)
					xqq    = sqrt(xqq)

					! Jordan et al 1999
					psimhs = -(0.7d0*hol &
								 + 0.75d0*(hol-14.3d0) &
								 * exp(-0.35d0*hol) + 10.7d0)
					psimh  = psimhs*stable &
									+ (c1 - stable)*psimhu(xqq)
					psixh  = psimhs*stable &
									+ (c1 - stable)*psixhu(xqq)

			! shift all coeffs to measurement height and stability
					rd = rdn / (c1+rdn/vonkar*(alz-psimh))
					rh = rhn / (c1+rhn/vonkar*(alz-psixh))
					re = ren / (c1+ren/vonkar*(alz-psixh))

			! update ustar, tstar, qstar using updated, shifted coeffs
					ustar = rd * vmag
					tstar = rh * delt
					qstar = re * delq

		enddo                     ! end iteration

		!------------------------------------------------------------
		! coefficients for turbulent flux calculation
		!------------------------------------------------------------
		! add windless coefficient for sensible heat flux
		! as in Jordan et al (JGR, 1999)
		!------------------------------------------------------------

			 shcoef = rhoa * ustar * cp * rh + c1
			 lhcoef = rhoa * ustar * Lheat  * re


	end subroutine sia2_env_atmo

! **********************************************************************



! **********************************************************************
!=======================================================================
!BOP
!
! !IROUTINE: atmo_boundary_const - compute coeeficients for atm-ice fluxes
!
!
! !INTERFACE:
!
	pure subroutine sia2_env_atmo_const(Tsf,wind,rhoa,lhcoef,shcoef)

! !DESCRIPTION:
!
! Compute coefficients for atm/ice fluxes, stress
! NOTE: \\
! (1) all fluxes are positive downward,  \\
! (2) reference temperature and humidity are NOT computed
!
! !REVISION HISTORY: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:

	real(kind=dbl_kind), intent(in) :: &
		 Tsf      , & ! surface temperature in Kelvin (K)
		 wind     , & ! wind speed (m/s)
		 rhoa         ! air density (kg/m^3)

	real(kind=dbl_kind), intent(out):: &
		 shcoef   , & ! transfer coefficient for sensible heat
		 lhcoef       ! transfer coefficient for latent heat

! internal vars

	real(kind=dbl_kind) :: &
		 Lheat   ! Lvap or Lsub, depending on surface type

	! from CICE 4 ice_constants.F90
	real(kind=dbl_kind), parameter :: &
		 cp_air    = 1005.0  ,& ! specific heat of air (J/kg/K)
		 Lsub      = 2.835e6 ,&   ! latent heat, sublimation freshwater (J/kg)
		 Lvap      = 2.501e6      ! latent heat, vaporization freshwater (J/kg)

	!------------------------------------------------------------
	! Initialize
	!------------------------------------------------------------

		 shcoef = 0.e0
		 lhcoef = 0.e0

	!------------------------------------------------------------
	! define variables that depend on surface type
	!------------------------------------------------------------
		 if (Tsf .ge. 0.e0) then
				 Lheat = Lvap           ! water to vapor
		 else
				 Lheat = Lsub           ! ice to vapor
		 endif

	!------------------------------------------------------------
	! coefficients for turbulent flux calculation
	!------------------------------------------------------------

		 shcoef = (1.20e-3)*cp_air*rhoa*wind
		 lhcoef = (1.50e-3)*Lheat*rhoa*wind

	end subroutine sia2_env_atmo_const


! **********************************************************************
! function sia2_ohf 
! Purpose: find ocean heat flux to ice, > 0 is melting flux
! ----------------------------------------------------------------------
	pure subroutine sia2_ohf(t_ocn,s_ocn,d_ocn,um_ocn,ud_ocn,um_ice, &
	ud_ice,fixed_mu,F_io)

	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only : &
			nc1,									&	! -1
			mu,										&	! linear liquidus constant for seawater
			cw 											! heat capacity of cold seawater (J g-1 K-1)
		use sia2_parameters, only : &
			chw, 									&	! ice-ocean heat transfer coefficient
			mu_w_min								! 
			
		implicit none
		
	! Function Arguments
	! --------------------------------------------------------------------
		real(kind=dbl_kind), intent(in) :: &
			t_ocn,									&	! ocean temperature (degC)
			s_ocn,									&	! ocean salinity (psu)
			d_ocn,									&	! ocean density (g m-3)
			um_ocn,									&	! ocean velocity magnitude (m/s)
			ud_ocn,									&	! ocean velocity direction (radians)
			um_ice,									&	! ice velocity magnitude (m/s)
			ud_ice										! ice velocity direction (radians)
		logical(kind=log_kind), intent(in) :: &
			fixed_mu								! downward longwave irradiance
		real(kind=dbl_kind), intent(out) :: &
			F_io										! ocean temperature (degC)

	! Internal Variables
	! --------------------------------------------------------------------
		real(kind=dbl_kind) :: &
			t_frz,									&	! freezing temp of ocean water
			tau_w,									&	! ice-ocean stress
			mu_w											! friction velocity

	! Function Code
	! --------------------------------------------------------------------
	
		! find mu
		mu_w = mu_w_min
		if (.not. fixed_mu) then
			tau_w = 2585.0_dbl_kind    ! placeholder - equation not yet implemented 
			mu_w = sqrt(abs(tau_w)/d_ocn)
			mu_w = max(mu_w,mu_w_min)
		endif

		! find freezing temperature of seawater
		t_frz = s_ocn*mu
		
		!t_frz = (-0.0575_dbl_kind +1.710523e-3_dbl_kind *sqrt(s_ocn) &
		!	-2.154996e-4_dbl_kind *s_ocn) * s_ocn  ! - 7.53e-4 *Db  ! <-- not including pressure part since this is at "surface"
		
		! find ocean heat flux, > 0 is melting....
		F_io = d_ocn*cw*chw*mu_w*(T_ocn - T_frz)
		!print *,d_ocn,cw,chw,mu_w,T_ocn, T_frz
		!stop

	end subroutine sia2_ohf


end module sia2_flux_heat

	
	
	
	
	
	
	