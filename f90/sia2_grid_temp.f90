module sia2_grid

	use sia2_constants
	use sia2_parameters
	use sia2_types
	use sia2_desalination
	
	implicit none

	! required parameters, globals
	! grid_model

	public 
	
	contains

! **********************************************************************
! sia2_create_ice: initialize ice type upon ice creation, using minumum
! values for kpp_eco_ice usage
!
! ----------------------------------------------------------------------
	subroutine sia2_create_ice(ice,iis,hi_new,hs_new,af_new,Tair, &
		t_ocn,s_ocn)
		
 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only: &
			c0,						& !
			c_5,					& !
			c1,						& !
			c1e6						!
		use sia2_parameters, only: 
			z_th_min,			& !
			z_max_snow,		& !
			z_max_ice,		& !
			s_const,			& !
			den_s_dry,		& !
			den_s_wet				!
		use sia2_types

		implicit none

	! Subroutine Arguments
	! --------------------------------------------------------------------
		type(ice_type), intent(inout) :: &
			ice 							  !
		integer(kind=int_kind), intent(in) :: &
			iis 								!
		real(kind=dbl_kind), intent(in) :: &
			hi_new,	 					& !
			hs_new,	 					& !
			Tair,	 						& !
			t_ocn,	 					& !
			s_ocn	 							!

	! Internal Variables
	! --------------------------------------------------------------------
		integer(kind=int_kind) :: &
			z_new,							&	!
			ii,									&	!
			int_z,							&	!
			sk_1,								&	!
			sk_z,								&	!
			z_new									!
		real(kind=dbl_kind) :: &
			tmp1,	 							& !
			t_interface,	 			& !
			snowfall_d,	 				& !
			snowfall_melt,	 		& !
			hs1	 							 !
		real(kind=dbl_kind), dimension(z_max_ice) :: &
			th_new, 						& !
			id_new								!

	! Subroutine Code
	! --------------------------------------------------------------------
		call sia2_null_ice(ice)
		ice%af = af_new

		!ice%pur = c0
		!ice%ed_w = c0
		!ice%snow%ed_w = c0

		! initially, assume no infiltration, but zero freeboard
		ice%fbh = c0

		! initially, snow thickness
		ice%sh_prev = hs_new
		
		! find mean ice surface temp from snow thickness
		tmp1 = 0.33*2.2/(2.2*hs_new + 0.33*hi_new) ! gamma in Arrigo et al. 1993 appendix
		t_interface = (0.2456008*Tair + tmp1*t_ocn)/(0.2456008+tmp1)
		t_interface = min(f(mi)%t,t_interface)
      
		! new snow density and 'melt' status
		if (Tair .lt. des_s_switch) then
			 snowfall_d = den_s_dry*c1e6
			 snowfall_melt = c0
		else
			 snowfall_d = den_s_wet*c1e6
			 snowfall_melt = c1
		endif
		
		! find new snow grid
		call sia2_new_grid(hs_new,z_th_min,5.e0,2,1,z_max_snow,.false., &
			th_new,id_new,z_new)
		ice%snow%depth = hs_new
		ice%snow%th = th_new
		ice%snow%z = z_new
				
	 	! fill new snow parameters
		hs1 = c0	 	
		do ii=1,z_new
			hs1 = hs1 + th_new(ii)*c_5  							   
	 
			! find snow state, interpolated between ice interface temp and airtemp	
			ice%snow%t(ii) = t_interface + hs1*(min(Tair,c0) - &
				t_interface)/ice%snow%depth
			ice%snow%d(ii) = snowfall_d
			ice%snow%heat = sia2_snow_heat(ice%snow%d(ii),ice%snow%t(ii))
			ice%snow%melt = snowfall_melt
										 
			! add other half to current layer depth
			hs1 = hs1 + th_new(ii)*c_5
		enddo  							   
		
		! find new ice grid			
		call sia2_new_grid(hi_new,z_th_min,10.e0,2,1,z_max_ice-sk_z,.false., &
			th_new,id_new,z_new)
		sk_1 = z_new + 1
		sk_z = z_new + z_sk
		int_z = z_new
		ice%z = sk_z

	 	! fill new ice parameters
	 	
		do ii=1,sk_z				
		
			if (ii .le. int_z) then

				! update thickness & ice depth
				ice(sc,ic,mi)%th(ii) = th_new(ii)
				ice(sc,ic,mi)%id(ii) = id_new(ii)
	
				! new ice temp & salinity
				if (iis .eq. 2) then
					! interpolate default new-ice salinity curve to get initial salinity
					call sia2_interpolate_salinity2( & 
					(ice%id(ii)-ice%th(ii)/2)/hi_new,ice%s(ii))
					! interpolate linearly between t_interface and water temp
					ice%t(ii) = t_interface - (t_interface-t_ocn) * (ice%id(ii) - &
					ice%th(ii)/c2)/hi_new
				elseif (iis .eq. 1) then
					! interpolate default new-ice salinity curve to get initial salinity
					call sia2_interpolate_salinity1( & 
					(ice%id(ii)-ice%th(ii)/2)/hi_new,ice%s(ii))
					! interpolate linearly between t_interface and water temp
					ice%t(ii) = t_interface - (t_interface-t_ocn) * (ice%id(ii) - &
					ice%th(ii)/c2)/hi_new
				else
					! constant ice salinity
					ice%s(ii) = s_const
					ice%t(ii) = t_ocn
				endif
			
			else 

				! constant skeletal layer  - this will be ignored in heat/mass
				! balance in KPP_ECO_ICE, but is required by regridding routines
				ice(sc,ic,mi)%th(ii) = z_th_min
				ice(sc,ic,mi)%id(ii) = id_new(ii-1) + z_th_min
			
				! constant skeletal layer 
				ice%s(ii) = s_ocn/c2
				ice%t(ii) = t_ocn
				
			endif
									
			! set new ice state
			call sia2_ice_state(ice%t(ii),ice%s(ii),ice%bs(ii),ice%bd(ii), &
				ice%d(ii),ice%bv(ii),ice%heat(ii))
				
		enddo
								
		! ignoring nutrients/tracers platelet layer, age, snow distribution for now!

	end subroutine sia2_create_ice

! **********************************************************************


! **********************************************************************
! ----------------------------------------------------------------------
! sia2_ice_grid finds new layer heights for a new total ice depth
!
! Requires:
! r_depth - the new depth to be assigned
!
! Results:
! th_new - contains new heights, up to th_z (all layers except skeletal layers)
! id_new - contains new ice depths for th_new
! int_z_new - new count of internal layers
! ----------------------------------------------------------------------
	subroutine sia2_new_grid(r_depth,h_min,h_max,grid_type,z_min,z_max, &
		ignore_h_max,th_new,id_new,z_new)
		
 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------

		implicit none

	! Subroutine Arguments
	! --------------------------------------------------------------------
		real(kind=dbl_kind), intent(inout) :: &
			r_depth 							  !
		real(kind=dbl_kind), intent(in) :: &
			h_min,	 							& !
			h_max										!
		integer(kind=int_kind), intent(in) :: &
			grid_type, 						& !
			z_min, 								& !
			z_max										!
		logical(kind=log_kind), intent(in) :: &
			ignore_h_max						!
		real(kind=dbl_kind), dimension(z_max), intent(out) :: &
			th_new, 							& !
			id_new									!
		integer(kind=int_kind), intent(out) :: &
			z_new								!
	
	! Internal Variables
	! --------------------------------------------------------------------
		integer(kind=int_kind) :: &
			ii,										&	!
			jj											!
		logical(kind=log_kind) :: &
			keep_growing						! 
		integer(kind=dbl_kind) :: &
			layer_divisor						!

		keep_growing = .true.         
!          z_last = ice(sc,ic,mi)%z
!          ice(sc,ic,mi)%snow%depth = r_depth   

!		 if (r_depth .gt. snow_min .and. r_depth .gt. z_th_min) then

		! test for max_depth
!		if ((r_depth .gt. h_max) .and. (.not. ignore_h_max)) then
		!    max_d_exceeded = max_d_exceeded + 1
		!    maxed_depth = .true.
	!			r_depth = h_max
!				print *,'Max depth exceeded - sc ic mi:',sc,ic,mi
		!else
		!    maxed_depth = .false.
!		endif        

		if (grid_type .eq. 0) then

!				if (r_depth .lt. h_min) then    
!						if ((r_depth+0.000001) .lt. h_min) then
!								print *,'Warning ***** minimum depth exceeded - unknown results may occur'
!								print *,'mi: ',mi,'in sia2_ice_grid line 26'
!						else
!								r_depth = h_min ! tollerate some floating point slop here, reset to h_min to avoid boundary creep
!						endif
!				endif

				z_new = z_max
											
				! divide up depth equally between layers 
				th_new = r_depth/z_new

		elseif (grid_type .ge. 1) then

				if ((grid_type .eq. 2) .and. (r_depth .le. h_crit)) then
						! determine ice layers required if less than z_int_max
						r_depth = r_depth - z_th_min*c2  ! top and bottom layers don't change thickness                   
						ii = z_min-1
						do while ((ii .le. z_max-2) .and. keep_growing)
								if (r_depth .le. z_th_min*dble(ii)) then
										ii = ii-1 ! take it back one layer so that we don't end up with layers less than z_th_min
										z_new = ii+2 ! include top and bottom layers that are at z_th_min always...

										! update middle layers
											 do jj=2,ii+1
												th_new(jj) = r_depth/dble(ii)
										enddo
										! update top/bottom layers
										th_new(1) = z_th_min
										th_new(z_new) = z_th_min
										keep_growing = .false.
								endif
								ii = ii+1
						enddo
						if (z_new .lt. z_max) then
								do ii=z_new+1,z_max
										th_new(ii) = c0
								enddo
						endif
				else

!						if (r_depth .lt. h_min .and. (grid_type .eq. 1)) then    
!								if ((r_depth+c_1e6) .lt. h_min) then
!										print *,'Warning ***** minimum depth exceeded - unknown results may occur'
!										print *,'mi: ',mi,'in sia2_ice_grid line 26'
!								else
!										r_depth = h_min ! tollerate some floating point slop here, reset to h_min to avoid boundary creep
!								endif
!						endif
				
						! number of layers is z_max
						z_new = z_max    

						! remove minimum depth & find unit to be added to 
						r_depth = r_depth - dble(z_new)*z_th_min
	
						jj = z_int_max/2
						! find divisor by which we will divide depth bins
						! geometric series sum (minus 2 b/c top and bottom are constant)
						layer_divisor = c2*(c1 - 1.3_dbl_kind**jj)/(c1 - 1.3_dbl_kind) - c2 

						! find new thicknesses
						th_new(1) = z_th_min
						th_new(z_new) = z_th_min
						do ii=2,z_new/2
								jj = ii-1
								th_new(ii) = z_th_min + (1.3_dbl_kind**jj)/layer_divisor*r_depth
								th_new(z_new-ii+1) = th_new(ii)
						enddo

				endif

		endif

!		 else ! end of r_depth > snow_min test
!		    ice(sc,ic,mi)%snow%z = 0 ! no snow layers
!		    ice(sc,ic,mi)%snow%depth = 0.d0  
!
!            ! vector assignments below
!            ice(sc,ic,mi)%snow%t = 0.d0 
!            ice(sc,ic,mi)%snow%d = 0.d0 
!            ice(sc,ic,mi)%snow%heat = 0.d0 
!            ice(sc,ic,mi)%snow%th = 0.d0
!            ice(sc,ic,mi)%snow%melt = 0.d0 
!		    
!		    ! save r_depth into sh_prev so mass is not wiped out
!		    ice(sc,ic,mi)%sh_prev = r_depth
!		    
!		 endif


           ! deal with skeletal thicknesses
!           if (z_last .eq. 0) then
!               ! new ice - create skeletal layer
!               do ii=int_z_new+1,int_z_new+z_sk
!                   th_new(ii) = sk_th_min
!               enddo
!           else
!               ! re-grid existing skeletal layer to new position
!               jj=z_last-z_sk+1
!               do ii=int_z_new+1,int_z_new+z_sk
!                   th_new(ii) = ice(sc,ic,mi)%th(jj)
!                   jj=jj+1
!               enddo
!           endif

           ! update id_new
           ! -----------------------------------------------------------
!           id_new(1) = th_new(1)
!           do ii=2,int_z_new+z_sk
!               id_new(ii) = id_new(ii-1) + th_new(ii)
!           enddo
	
	end subroutine sia2_new_grid

! **********************************************************************


! **********************************************************************
! sia2_pre_grid: calculate grid surface options, new thickness, and new
! snow to be added if porous ice was drained above freeboard 
! ----------------------------------------------------------------------

	subroutine sia2_pre_grid(ice,f,m,flooded,c_gl,s_gl, &
		r_depth,f_depth,did_flood,flood_1x,t_flood,s_flood, &
		melt_flood,melt_drain,did_drain,sn_depth,new_snow)

 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_state

		implicit none
		
	! Function Arguments
	! --------------------------------------------------------------------
		type(ice_type), intent(in) :: &
			ice											!
		type(forcing_type), intent(in) :: &
			f												!
		type(meta_type), intent(inout) :: &
			m												!
		real(kind=dbl_kind), intent(in) :: &
			flooded									!
		real(kind=dbl_kind), intent(inout) :: &
			c_gl,									& !
			s_gl										!
		real(kind=dbl_kind), intent(out) :: &
			r_depth,							& !
			f_depth,							& !
			t_flood,							& !
			s_flood, 							& !
			melt_flood,						& !
			melt_drain							!
		logical(kind=log_kind), intent(out) :: &
			did_flood,						& !
			flood_1x,							& !
			did_drain								!
		type(snow_type), intent(out) :: &
			new_snow								!

	! Internal Variables
	! --------------------------------------------------------------------
		integer(kind=int_kind) :: &
			ii,										& !
			jj,										& !
			z_snow,								& !
			int_z										!
		real(kind=dbl_kind) :: &
			top_f,								& !
			bot_f,								& !
			tmp1,									& !
			d_mean,								& !
			heat_mean								!

	! Subroutine Code
	! --------------------------------------------------------------------

	 	! find new thickness/icedepth structure with th_new and id_new. flooded will be added later
	 	int_z = ice%z - z_sk
	 	r_depth = ice%id(int_z) + c_gl + s_gl 

		! check to make sure we don't grow below minimum
!		if (r_depth .ge. h_min) then

				   ! record bm_lost and freshwater input for melting
!				   if (c_gl .lt. 0.d0) then
!					   dz_mean = abs(c_gl)
!					   ii = sk_z
!					   do while(dz_mean .gt. 0.d0)
!						   smalg_lost = ice%smalg(ii) - min_alg
!						   smalg_lost = max(0.,smalg_lost)
!						   if (dz_mean .gt. ice%th(ii)) then
!							  ! converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
!							  m(mi)%bm_lost =  m(mi)%bm_lost + smalg_lost* &
!							  ice%th(ii)*cell_area*1.e3*ice%af   ! g/pixel
!							  dz_mean = dz_mean - ice%th(ii)
!							  ii = ii-1
!						   else
!							  ! converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
!							  m(mi)%bm_lost =  m(mi)%bm_lost + smalg_lost* &
!							  dz_mean*cell_area*1.e3*ice%af   ! g/pixel
!							  dz_mean = 0.d0
!						   endif
!					   enddo
!				   endif


		! Determine surface flooding - only one type of surface flooding
		! may occur at once (seawater flooding or meltwater flooding)
		! and seawater flooding takes precedence
	               
		melt_flood = c0
		melt_drain = c0
		z_snow = ice%snow%z
				   
!		if ((flooded .ge. fl_max .and. ice%bv(1) .ge. vb_crit) &
!			.or. (flooded .ge. fl_max*2.)) then				   
		if (flooded .ge. fl_max) then
					  
			did_flood = .true.
			flood_1x = .false.

			! compress snow to minimum of 1/2 pure ice density
			tmp1 = max(IceD/c2,ice%snow%d(1))
			f_depth = flooded*ice%snow%d(1)/tmp1
			
			! change new ice depth
			r_depth = r_depth + f_depth
			
			! record ice volume change due to snow-ice gain
			m%snow_growth = m%snow_growth + f_depth* &
				ice%af*cell_area*c1e6

			! Conservation of heat: mix snow and seawater
			! assuming mass/volume changes are negligable inside a layer initially...
			!top_f = snow_mass/snowh1m/IceD  ! volume fraction of fresh ice in snow
			
			top_f = tmp1/IceD ! snow/water volume ratio					  
			bot_f = c1-top_f ! water/snow volume ratio

			s_flood = max(min_sal,f%s*bot_f)
			d_mean = f%d*bot_f + IceD*top_f

		 ! snow heat is assumed to be just latent heat of melting + 0.02*t*d
			heat_mean = cw*f%t*f%d*bot_f - &
				(ice%snow%t(z_snow)*0.02_dbl_kind + Lf)*IceD*top_f
			!dz_mean = 1.
			
			t_flood = sia2_ice_temp(s_flood,d_mean,heat_mean)

			!call sia2_ice_state(t_flood,s_flood,bs_flood,bd_flood,d_flood, &
			!	bv_flood,heat_flood		


		elseif (use_drain .eq. 1) then
				
			 jj=0
			 new_snow%z = 0
			 new_snow%th = c0
			 new_snow%d = c0
			 new_snow%heat = c0
			 do ii=1,ice%z
					if (ice%bv(ii) .gt. bv_conv .and. &
							ice%id(ii) .le. ice%fbh) then
							jj = ii
					endif
			 enddo

			 ! drain ice layers above and eq. to jj of brine, turn it into snow
			 if (jj .gt. 0) then
					 did_drain = .true.
					 new_snow%z = jj
					 do ii=1,jj
							 melt_drain = melt_drain + ice%th(ii)
							 new_snow%th(ii) = ice%th(ii)
							 new_snow%d(ii) = (c1-ice%bv(jj)*c_001)*IceD
							 new_snow%heat(ii) = sia2_snow_heat(new_snow%d(ii),ice%t(ii))
					 enddo
					 new_snow%depth = melt_drain												 

					 ! update s_gl to drop layer from ice (and move to snow)
					 s_gl = s_gl - melt_drain
					 
					 r_depth = r_depth - melt_drain

			 endif ! end of is there draining?

    endif
    
    ! deal with snow 
    sn_depth = max(c0,ice%snow%depth + snow_gl)

! snow melt in KPP_ECO_ICE goes directly to mixed layer
		! add fresh water to melt water pond
!		if ((use_ponds .eq. 1) .and. (snow_melt .gt. 0.) .and. &
!		(z_snow .gt. 0)) then
!			snow_melt_now = snow_melt
!			new_layer_top = c0
!#include "sia2_env_pond_update.inc.f90"               
!			 snow_melt = c0			
!		endif
		
		if (did_flood) then
!      if (validation .eq. 1) then
!        ice(sc,ic,mi)%sh_offset = ice(sc,ic,mi)%sh_offset + flooded
!      endif
			sn_depth = max(sn_depth - flooded, c0)
		elseif (did_drain) then
      ! draining surface ice layer and turning to snow
			if(melt_drain .gt. c0) then
				sn_depth = sn_depth + melt_drain
!				if (validation .eq. 1) then
!						ice(sc,ic,mi)%sh_offset = ice(sc,ic,mi)%sh_offset + th_snow_new(ii)
!				endif
		endif

    

	end subroutine sia2_pre_grid

! **********************************************************************


! **********************************************************************
! ----------------------------------------------------------------------
! subroutine: sia2_ice_remap
! Purpose: maps ice to new boundaries
!
! INPUT VARIALBES:
!
! ice - current physical ice vars
! bio - current biology vars
! int_z - current internal ice depth
! th_new - new thicknesses
! id_new - new ice depths of layers
! int_z_new - the new number of internal layers
! new_layer_top - where among the old layers that the new ice column should 
!     begin (used during sublimation)
! c_gl - congelation ice growth layer (m), to be added to bottom internal layer
!    If c_gl is greater than zero, the salinity of the newly grown ice is also required
!    as the variable:
! c_gl_sal - salinity of congelation growth layer 
! flooded - if part of the ice column has been flooded, this is the depth of 
!    flooding in meters. If flooded is greater than zero, then these 
!    additional variables must be defined for assignment to flooded 
!    portions of layers:
! T_flood - temp (degC)
! S_flood - salinity (psu)
! d_flood - density of bulk ice (g/m^3)
! bs_flood - brine salinity (psu)
! bd_flood - density of brine (g/m^3)
! bv_flood - brine volume (ppt)
! heat_flood - heat required to bring ice from T_ref to T_flood
!
! OUTPUT VARIALBES:
!
! ice(mi)th
! ice(mi)id
! ice(mi)t
! ice(mi)s
! ice(mi)d
! ice(mi)bs
! ice(mi)bd
! ice(mi)bv
! ice(mi)heat
! ice(mi)smalg
! ice(mi)NO3
! ice(mi)NH4
! ice(mi)PO4
! ice(mi)SiOH4
!
! INTERNAL VARS REQUIRED:
! 
! integer :: ii,jj
! real(kind=dbl_kind) :: interim_top,interim_bot,old_layer_top,new_layer_bot
! real(kind=dbl_kind) :: old_layer_bot,z1,z2,dz,heat_total,dz_total
! real(kind=dbl_kind) :: t_mean,s_mean,d_mean,bs_mean,bd_mean,bv_mean,heat_mean
! real(kind=dbl_kind), dimension(z_max) :: t_new,s_new,d_new,bs_new,bd_new
! real(kind=dbl_kind), dimension(z_max) :: bv_new,smalg_new
! real(kind=dbl_kind), dimension(z_max+1) :: NO3_new,NH4_new,PO4_new,SiOH4_new
!
! ----------------------------------------------------------------------

	subroutine sia2_ice_remap(ice,f,m,int_z_new,th_new,id_new,t_flood, &
		s_flood,flooded,melt_flood,s_gl,c_gl,c_gl_sal,dhdt_cm_s,f_depth)
	
 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_state

		implicit none
		
		! Function Arguments
		! --------------------------------------------------------------------
		type(ice_type), intent(inout) :: &
			ice											!
		type(forcing_type), intent(in) :: &
			f												!
		type(meta_type), intent(inout) :: &
			m												!
		integer(kind=int_kind), intent (in) :: &
			int_z_new								!
		real(kind=dbl_kind), dimension(z_max_ice), intent(in) :: &
			th_new,								&	!
			id_new									!
		real(kind=dbl_kind), intent(in) :: &
			t_flood,							&	!
			s_flood,							&	!
			flooded,							&	!
			melt_flood,						&	!
			dhdt_cm_s								!
		real(kind=dbl_kind), intent(inout) :: &
			s_gl,									&	!
			c_gl,									&	!
			c_gl_sal,							&	!
			f_depth									!
	
		! Internal Variables
		! --------------------------------------------------------------------
		logical(kind=log_kind) :: &
			remapping,						&	!
			alg_fixed  							! 
		integer(kind=int_kind) :: &
			ii,										&	!
			jj,										&	!
			int_z,								&	!
			sk_1,									&	!
			sk_z										!		
		real(kind=dbl_kind) :: &
			mf_depth,  						&	! 
			tmp1,									&	!
			tmp2,									&	!
			tmp3,									&	!
			dz,										&	!
			z1,										&	!
			z2,										&	!
			new_layer_top, 				&	!
			old_layer_top, 				&	!
			new_layer_bot, 				&	!
			old_layer_bot, 				&	!
			interim_bot, 					&	!
			interim_top,	 				&	!
			dz_total, 						&	!
			dz_mean, 							&	!
			heat_total, 					&	! 
			d_flood,  						&	!
			bs_flood, 						&	!
			bd_flood, 						&	!
			bv_flood, 						&	!
			t_mean, 							&	! 
			s_mean,  							&	!
			d_mean,  							&	!
			bs_mean, 							&	!
			bd_mean, 							&	!
			bv_mean, 							&	!
			heat_mean, 							&	!
			heat_flood  						!
		
		real(kind=dbl_kind), dimension(z_max_ice) :: &
			d_new, 							&	!
			t_new, 							&	!
			s_new, 							&	!
			debug_z, 						&	!
			NO3_new, 						&	!
			NH4_new, 						&	!
			PO4_new, 						&	!
			SiOH4_new, 					&	!
			smalg_pre_map, 			&	!
			smalg_new, 					&	!
			poc_new 							!
		

		! find snow flooding/draining modification parameters
		!f_depth = flooded*ice%snow%d(1)/max(IceD/c2,ice%snow%d(1))
		mf_depth = melt_flood

		if (f_depth .gt. c0) then
			! record ice volume change due to snow-ice gain
			m%snow_growth = m%snow_growth + f_depth* &
				ice%af*cell_area*c1e6		
			call sia2_ice_state(t_flood,s_flood,bs_flood,bd_flood, &
				d_flood,bv_flood,heat_flood)
		endif
			
		! setup vars to track new and old layer heights								
		old_layer_top = c0
		new_layer_bot = c0
		old_layer_bot = c0
	
		! record flooded into var that we can change during boundary adjustment
		f_depth = flooded*ice%snow%d(1)/max(IceD/2.d0,ice%snow%d(1))
		mf_depth = melt_flood
	
		! determine if we are to move algae with growing ice
		if ((c_gl .lt. c0) .or. (alg_mig .eq. 0) .or. (flooded .gt. c0)) then
			alg_fixed = .true.
		else
			alg_fixed = .false.
		endif						
	
		! don't grow if depth will be maxed out by growth
		if (.NOT. maxed_depth) then
				 
			! need to know when re-mapping of old ice column begins - set var to tell when this happens
			remapping = .false.
	
			! calculate algal migration, and store in an array for latter use in
			! re-mapping
			tmp1 = dhdt_cm_s*86400.
			if ((alg_mig .eq. 1) .and. (abs(tmp1) .lt. 1.5) .and. (c_gl .ne. c0) .and. dhdt_cm_s .gt. c0) then
				! remap algae into temp array that reflects migration after
				! bulk ice remapping below
				tmp1 = c_gl*(1.-abs(tmp1)/1.5) ! height of "retainment"
				tmp2 = (ice%id(sk_z)+tmp1)/ice%id(sk_z) ! percentage height change of virtual ice
				dz = c0
				do ii=1,sk_z
					if (tmp1 .lt. c0) then
						dz = dz + ice%th(ii)*abs(1.-tmp2)	 ! increment inter-layer size
						tmp3 = dz/ice%th(ii) ! precentage of layer that reaches into adjacent layer
						! case where melting is occuring, but slowly enough that algae can retreat into ice
						if (ii .lt. sk_z) then
							smalg_pre_map(ii) = (ice%smalg(ii)*(1.-tmp3) + &
							ice%smalg(ii+1)*tmp3)/tmp2
						else
							smalg_pre_map(ii) = ice%smalg(ii)/tmp2 ! just use concentrated algae for last layer
						endif
					else
						! case where we are growing, but slowly enough the algae can remain at ice bottom
						if (ii .eq. 1) then
							smalg_pre_map(ii) = ice%smalg(ii)/tmp2 ! just use diluted algae for 1st layer
						else
							dz = dz + ice%th(ii-1)*abs(1.-tmp2)	 ! increment inter-layer size
							tmp3 = dz/ice%th(ii) ! precentage of layer that reaches into adjacent layer
							smalg_pre_map(ii) = (ice%smalg(ii)*(1.-tmp3) + &
							ice%smalg(ii-1)*tmp3)/tmp2
						endif
					endif
				enddo
			else
					! default for when algae is not migrating
					smalg_pre_map = ice%smalg
			endif
				 
			do ii=1,int_z_new+z_sk
	
				! These vars keep track of new ice physics while layers are mixing
				dz_total = c0
				d_new(ii) = c0				 
				s_new(ii) = c0			
				NO3_new(ii) = c0
				NH4_new(ii) = c0
				PO4_new(ii) = c0
				SiOH4_new(ii) = c0
				smalg_new(ii) = c0
				poc_new(ii) = c0
				heat_total = c0
	
				! flood top layers first, then integrate old ice into new
				if (f_depth .gt. c0) then
					f_depth = f_depth - th_new(ii)
					if (f_depth .lt. c0) then
						dz = f_depth + th_new(ii)
					else
						dz = th_new(ii)
					endif
					 
					! record tracers
					tmp1 = S_flood*dz
					tmp2 = tmp1/f(mi)%s
					s_new(ii) = s_new(ii) + tmp1 
					d_new(ii) = d_new(ii) + d_flood*dz
					heat_total = heat_total + heat_flood*dz
					NO3_new(ii) = NO3_new(ii) + f(mi)%no3*tmp2
					NH4_new(ii) = NH4_new(ii) + f(mi)%nh4*tmp2
					PO4_new(ii) = PO4_new(ii) + f(mi)%po4*tmp2
					SiOH4_new(ii) = SiOH4_new(ii) + f(mi)%sioh4*tmp2
					poc_new(ii) = poc_new(ii) + f(mi)%poc*tmp2
					smalg_new(ii) = smalg_new(ii) + alg_wc*dz
	
					! keeping track of emerging new layer thickness 
					dz_total = dz_total + dz
					 
					! assume layers that flood have a kick-ass brine channel network for flooding/desal
					ice%bced(ii) = 1
					 
	!			elseif (mf_depth .gt. 0.) then
	!				mf_depth = mf_depth - th_new(ii)
	!				if (mf_depth .lt. 0.) then
	!					dz = mf_depth + th_new(ii)
	!				else
	!					dz = th_new(ii)
	!				endif
	!				 
	!				! record tracers
	!				z_snow = ice%snow%z
	!				tmp1 = 1. - ice%snow%d(z_snow)/IceD
	!				 
	!				s_new(ii) = s_new(ii) + S_flood*dz 
	!				d_new(ii) = d_new(ii) + d_flood*dz
	!				heat_total = heat_total + heat_flood*dz
	!				NO3_new(ii) = NO3_new(ii) + no3_pond_new*dz*tmp1
	!				NH4_new(ii) = NH4_new(ii) + nh4_pond_new*dz*tmp1
	!				PO4_new(ii) = PO4_new(ii) + po4_pond_new*dz*tmp1
	!				SiOH4_new(ii) = SiOH4_new(ii) + sioh4_pond_new*dz*tmp1
	!				poc_new(ii) = poc_new(ii) + poc_pond_new*dz*tmp1
	!				smalg_new(ii) = smalg_new(ii) + smalg_pond_new*dz*tmp1
	!
	!				! keeping track of emerging new layer thickness 
	!				dz_total = dz_total + dz
	
				endif					 
	
				if (dz_total .lt. th_new(ii)) then
	
					! find heights of new layer with regard to old thicknesses
					! new_layer_top should be zero even if flooding has occured,
					! because it is really a reference to the old ice pack.
					if (.NOT. remapping) then
						remapping = .TRUE.
						if (new_layer_top .lt. c0) then
							new_layer_top = c0
						endif
							new_layer_bot = new_layer_top + th_new(ii) - dz_total
					else
						new_layer_top = new_layer_bot
						new_layer_bot = new_layer_bot + th_new(ii)											 
					endif
	
					! find 1st old layer that contains part of new layer, going down 
					jj = 0
					old_layer_bot = c0
					do while(new_layer_top .ge. old_layer_bot .and. jj .le. sk_z)
						jj=jj+1
						if (jj .eq. 43) then
	!											print *,'jj exceeded max layers in ice remapping'
	!												 print *,'mi: ',mi,'layer: ',ii,'in sia2_env_ice_remap line 228'
						endif
						! old thickness...
						old_layer_bot = old_layer_bot + ice%th(jj)
					enddo
					old_layer_top = old_layer_bot - ice%th(jj)
					! now jj is OLD layer where NEW layer ii starts...
	
					! find total heat/salt from multiple layers/partial layers that make up
					! the new layer
					do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. sk_z))
					
						
						! ----- NEW LAYER GEOMETRIES ------------------------------		
						interim_top = max(old_layer_top,new_layer_top)
						interim_bot = min(old_layer_bot,new_layer_bot)
	
						z1 = interim_top - old_layer_top
						z2 = interim_bot - old_layer_top
						dz = z2-z1
							
						if (ii .le. int_z_new) then
							if (jj .le. int_z .and. ii .le. int_z_new) then
							! redistribution of congelation ice to congelation ice (basal growth)
								tmp1 = dz
								s_mean = ice%s(jj)
								d_mean = ice%d(jj)
								heat_mean = ice%heat(jj)
							else
							! redistribution of skeletal ice into congelation ice (basal growth)
								tmp1 = min(1.,c_gl_sal/ice%s(jj))*dz ! remove brine volume
								t_mean = f(mi)%t
								s_mean = c_gl_sal
								! find heat_mean
								call sia2_ice_state(t_mean,s_mean,bs_mean,bd_mean,d_mean, &
									bv_mean,heat_mean)       
							endif
						else
							! redistribution of congelation ice to skeletal ice (basal melt)
							tmp1 = dz
							t_mean = f(mi)%t
							s_mean = f(mi)%s/c2
							! find heat_mean
							call sia2_ice_state(t_mean,s_mean,bs_mean,bd_mean,d_mean, &
								bv_mean,heat_mean)       
						endif
																			
						! record tracers into new layer
						s_new(ii) = s_new(ii) + s_mean*dz 
						d_new(ii) = d_new(ii) + d_mean*dz
						heat_total = heat_total + heat_mean*dz
						no3_new(ii) = no3_new(ii) + ice%no3(jj)*tmp1
						nh4_new(ii) = nh4_new(ii) + ice%nh4(jj)*tmp1
						po4_new(ii) = po4_new(ii) + ice%po4(jj)*tmp1
						sioh4_new(ii) = sioh4_new(ii) + ice%sioh4(jj)*tmp1
						poc_new(ii) = poc_new(ii) + ice%poc(jj)*tmp1
						! assuming algae doesn't dilute(is not lost) during skeletal->congelation freeze
						smalg_new(ii) = smalg_new(ii) + smalg_pre_map(jj)*dz
	
						! setup variable for next partial layer addition
						dz_total = dz_total + dz 
	
						jj=jj+1
						if (jj .le. sk_z) then
							! find boundaries of next old layer
							old_layer_top = old_layer_bot
							old_layer_bot = old_layer_bot + ice%th(jj)																										
						endif
	
					enddo
	
					if (dz_total .lt. th_new(ii) .and. (jj .gt. sk_z)) then
					! making new ice out of seawater
	
						if (ii .le. int_z_new .and. c_gl .gt. c0) then ! these tests should alway be true together
						! new congelation ice out of seawater - fast growth!!
							dz = max(c_gl,th_new(ii)-dz_total)
							t_mean = f(mi)%t
							s_mean = c_gl_sal						 
						else
							! new skeletal ice out of seawater
							dz = th_new(ii) - dz_total
							t_mean = f(mi)%t
							s_mean = f(mi)%s/c2
						endif			 
						! fraction seawater concentrations by new brine volume
						call sia2_ice_state(t_mean,s_mean,bs_mean,bd_mean,d_mean, &
							 bv_mean,heat_mean)       					
						tmp1 = dz*bv_mean*c_001 
	
						! record tracers into new layer
						s_new(ii) = s_new(ii) + s_mean*dz 
						d_new(ii) = d_new(ii) + d_mean*dz
						heat_total = heat_total + heat_mean*dz
						no3_new(ii) = no3_new(ii) + f%no3*tmp1 ! extending into mixed layer 
						nh4_new(ii) = nh4_new(ii) + f%nh4*tmp1 ! extending into mixed layer
						po4_new(ii) = po4_new(ii) + f%po4*tmp1 ! extending into mixed layer
						sioh4_new(ii) = sioh4_new(ii) + f%sioh4*tmp1 ! extending into mixed layer
						smalg_new(ii) = smalg_new(ii) + alg_wc*tmp1
						poc_new(ii) = poc_new(ii) + f(mi)%poc*tmp1
						
						! update boundaries for next layer					
						c_gl = c_gl - dz
						dz_total = th_new(ii) ! round out layer - this is the last way to add mass
						
					endif
	
					! round out dz_total whether we are shrinking or growing,
					! otherwise division below can create/take away tracers
					dz_total = th_new(ii)
	
				endif	 ! end of "if dz_total .lt. th_new(ii)								 
	
				if (ii .le. int_z_new) then
					! solve for new layer temp. from heat/salt/density	
					dz_mean = dz_total
	!				if (dz_total .eq. 0.) then
	!					print *,mi,ic,i_ease,j_ease,'**********************************************'
	!					print *,'t---',ice%t
	!					print *,'t_prev---',t_prev
	!					print *,'t_last---',t_last
	!					print *,'t_next---',t_next
	!					print *,'ki---',ki
	!				endif												 
					d_new(ii) = d_new(ii)/dz_total
					s_new(ii) = s_new(ii)/dz_total
					d_mean = d_new(ii)
					s_mean = s_new(ii)
					heat_mean = heat_total ! division by dz is included in calc below
					if (s_mean .le. 0.) then
						print *,'Salinity (smean) is less than zero.'
						print *,'mi: ',mi,'layer: ',iii,'in sia2_env_ice_remap line 363'
					endif								 
			
					t_mean = sia2_ice_temp(s_mean,d_mean,heat_mean)
		
					! check to make sure we don't get warmer than the ocean
					!if (ice%id(ii) .le. ice%fbh) then
					!		t_mean = min(t_mean,S_mean*mu)
					!else
					!		t_mean = min(t_mean,f(mi)%t)
					!endif
		
					t_new(ii) = t_mean
				
				else
					! we know skeletal temp, save computations
					t_new(ii) = f(mi)%t
					s_new(ii) = s_new(ii)/dz_total
					d_new(ii) = d_new(ii)/dz_total
				endif
	
			enddo ! end of new grid portioning
	
			! update layer indices
			int_z = int_z_new
			sk_1 = int_z + 1
			sk_z = int_z + z_sk
			ice%z = sk_z
	
			! update ice structure
			do ii=1,sk_z
			
				! assign new layer boundaries
				ice%th(ii) = th_new(ii)
				ice%id(ii) = id_new(ii)
	
				! assign new ice physics
				ice%t(ii) = t_new(ii)
				ice%s(ii) = s_new(ii)
				call sia2_ice_state(t_new(ii),s_new(ii),ice%bs(ii),ice%bd(ii), &
					ice%d(ii),ice%bv(ii),ice%heat(ii))
				
				! update tracers
				ice%no3(ii) = NO3_new(ii)/th_new(ii)
				ice%nh4(ii) = NH4_new(ii)/th_new(ii)
				ice%po4(ii) = PO4_new(ii)/th_new(ii)
				ice%sioh4(ii) = SiOH4_new(ii)/th_new(ii)
				ice%poc(ii) = poc_new(ii)/th_new(ii)
				ice%smalg(ii) = smalg_new(ii)/th_new(ii)
	
			enddo
	
	!		dz_sk = c0
	!#include "sia2_env_update_skeletal.inc.f90"
	
			! adjust surface temp if there was no snow, but flooding occured (to prevent numerical instability)
			! uh, this can't happen...
			!if (flooded .gt. 0 .and. m(mi)%z_snow .eq. 0) then
			!		 ice%ts = ice%t(1)	! make equal to 1st layer temp
			!endif
							
		endif	 ! end of "don't grow if depth will be maxed out by growth"
	
	end subroutine sia2_ice_remap

! **********************************************************************


! **********************************************************************

	! Adjust boundaries of snow layers and calculate new layer temperatures
	! using snow heat capacity to conserve total heat.  
	! ------------------------------------------------------------------

	subroutine sia2_snow_remap(ice,f,new_snow,z_snow_new,th_new, &
		flooded,melt_drain)
	
 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_state

		implicit none
		
		! Function Arguments
		! --------------------------------------------------------------------
		type(ice_type), intent(inout) :: &
			ice											!
		type(forcing_type), intent(in) :: &
			f												!
		type(snow_type), intent(in) :: &
			new_snow								!
		integer(kind=int_kind), intent (in) :: &
			z_snow_new							!
		real(kind=dbl_kind), dimension(z_max_snow), intent(in) :: &
			th_new								 	!
		real(kind=dbl_kind), intent(in) :: &
			flooded 							 	!
		real(kind=dbl_kind), intent(inout) :: &
			melt_drain							!
	
		! Internal Variables
		! --------------------------------------------------------------------
		logical(kind=log_kind) :: &
			keep_growing  					! 
		integer(kind=int_kind) :: &
			ii,										&	!
			jj,										&	!
			z_last									!		
		real(kind=dbl_kind) :: &
			tmp1,									&	!
			dz,										&	!
			z1,										&	!
			z2,										&	!
			new_layer_top, 				&	!
			old_layer_top, 				&	!
			new_layer_bot, 				&	!
			old_layer_bot, 				&	!
			interim_bot, 					&	!
			interim_top,	 				&	!
			dz_total, 						&	!
			heat_total, 					&	! 
			t_mean, 							&	! 
			d_mean,  							&	!
			airtemp1c,  					&	!
			snowfall_d,  					&	!
			snowfall_melt,  			&	!
			heat_mean 							!
		
		real(kind=dbl_kind), dimension(z_max_snow) :: &
			d_new, 							&	!
			t_new, 							&	!
			debug_z, 						&	!
			melt_new 							!
		
		! Subroutine Code
		! --------------------------------------------------------------------

		z_last = ice%snow%z        ! 
		keep_growing = .TRUE.  ! set this variable to trigger rest of layer-tracking vars
		
		if (z_snow_new .gt. 0) then
						
			! setup vars to track new and old layer heights               
			old_layer_top = d0_
			new_layer_bot = d0_
			old_layer_bot = d0_
	
			! new snow density and 'melt' status
			airtemp1c = f%at - kelvin0
			if (airtemp1c .lt. des_s_switch) then
				 snowfall_d = den_s_dry*c1e6
				 snowfall_melt = c0
			else
				 snowfall_d = den_s_wet*c1e6
				 snowfall_melt = c1
			endif
			
			! allocate new layer tracers
			do ii=1,z_snow_new			   
			
				! These vars keep track of new ice physics while layers are mixing
				dz_total = c0
				d_new(ii) = c0
				heat_total = c0
				melt_new(ii) = c0
				
				if (melt_drain .gt. c0) then
			
					 ! re-adjust new layer bottom for next snow layer temp adjustment
					 new_layer_top = new_layer_bot
					 new_layer_bot = new_layer_bot + th_new(ii) ! th_new comes from snow_grid calc
	
					 ! find 1st old layer that contains part of new layer, going down 
					 jj = new_snow%z+1
					 old_layer_bot = c0
					 do while((new_layer_top .ge. old_layer_bot) .and. (jj .gt. 0))
							 jj=jj-1
							 ! old thickness...
							 old_layer_bot = old_layer_bot + new_snow%th(jj)
					 enddo
	
					 ! add heat to new layer from old layer 
					 if ((jj .le. new_snow%z) .and. (old_layer_bot .gt. new_layer_top)) then
							 old_layer_top = old_layer_bot - new_snow%th(jj)
							 ! jj is OLD layer where NEW layer ii starts...
				
							 ! add heat to new layer from old layer
							 do while ((old_layer_top .lt. new_layer_bot) .and. (jj .gt. 0))
									 ! make sure that drained layer still has some ice left in it
									 if (new_snow%d(jj) .gt. c0) then
											 ! ----- NEW LAYER GEOMETRIES ------------------------------ 
																 
											 interim_top = max(old_layer_top,new_layer_top)
											 interim_bot = min(old_layer_bot,new_layer_bot)
		
											 z1 = interim_top - old_layer_top
											 z2 = interim_bot - old_layer_top
											 dz = z2-z1
		
											 ! record density*thickness, heat*thickness
											 d_new(ii) = d_new(ii) + new_snow%d(jj)*dz
											 heat_total = heat_total + new_snow%heat(jj)*dz
											 melt_new(ii) = melt_new(ii) + dz
		
											 ! ----- SETUP VARIABLES FOR NEXT PARTIAL LAYER ------------
											 ! keeping track of emerging new layer thickness, for 
											 dz_total = dz_total + dz 
		
									 endif
									 ! find boundaries of next old layer
									 jj=jj-1
									 old_layer_top = old_layer_bot
									 old_layer_bot = old_layer_bot + ice%snow%th(jj)
							 enddo
	
						endif
			 
				endif
							 
			 
				if (dz_total .lt. th_new(ii)) then				       
			
				 if (keep_growing) then
						 ! setup vars to track new and old layer heights               
						 old_layer_top = c0
						 if (flooded .gt. c0) then
								 new_layer_bot = flooded ! start new layers higher up if some snow flooded
						 else
								 new_layer_bot = c0
						 endif
						 old_layer_bot = c0
						 keep_growing = .FALSE.
						 melt_drain = c0
				 endif
			
				 ! add newly-fallen snow to layer
				 if (z_last .eq. 0) then
					
					 ! fill layer with new snow
					 dz = th_new(ii) - dz_total						   
				 
					 ! new snow density and 'melt' status
					 !if (airtemp1c .lt. des_s_switch) then
					 !	 d_mean = den_s_dry*c1e6
					 !else
					!	 d_mean = den_s_wet*c1e6
					!	 melt_new(ii) = melt_new(ii) + dz
					 !endif
					 d_mean = snowfall_d
					 melt_new(ii) = snowfall_melt
					 
					 ! find median depth of new snow sub-layer, store in tmp1
					 tmp1 = c0  							   
					 do jj=1,ii
						 tmp1 = tmp1 + th_new(jj)
					 enddo
					 tmp1 = tmp1 - dz_total - th_new(ii)*c_5
					 
					 ! find temp of snow, interpolated between surface temp and airtemp
					 t_mean = ice%snow%ts + tmp1*(min(airtemp1c,c0) - &
						 ice%snow%ts)/ice%snow%depth
											 
					 ! find heat mean							      
					 !t_mean = t_mean + kelvin0 
					 !heat_mean = d_mean*(heat_snow0 - 0.2309*T_mean - &
					 !	0.0034*T_mean**2)
					 heat_mean = sia2_snow_heat(d_mean,t_mean)
				 
					 ! record tracers (melt already done above)
					 d_new(ii) = d_new(ii) + d_mean*dz
					 heat_total = heat_total + heat_mean*dz
					 dz_total = th_new(ii)						       						   
				 
				 ! map old snow to new snow layer
				 else	
			
					 ! re-adjust new layer bottom for next snow layer temp adjustment
					 new_layer_top = new_layer_bot
					 new_layer_bot = new_layer_bot + th_new(ii) ! th_new comes from snow_grid calc
			
					 ! find 1st old layer that contains part of new layer, going down 
					 jj = 0
					 old_layer_bot = c0
					 do while((new_layer_top .ge. old_layer_bot) .and. (jj .lt. z_last))
						 jj=jj+1
						 ! old thickness...
						 old_layer_bot = old_layer_bot + ice%snow%th(jj)
					 enddo
			
					 ! add heat to new layer from old layer 
					 if ((jj .le. z_last) .and. (old_layer_bot .gt. new_layer_top)) then
						 old_layer_top = old_layer_bot - ice%snow%th(jj)
						 ! jj is OLD layer where NEW layer ii starts...
				
						 ! add heat to new layer from old layer
						 do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. z_last))
							 ! ----- NEW LAYER GEOMETRIES ------------------------------ 
			
							 interim_top = max(old_layer_top,new_layer_top)
							 interim_bot = min(old_layer_bot,new_layer_bot)
			
							 z1 = interim_top - old_layer_top
							 z2 = interim_bot - old_layer_top
							 dz = z2-z1
			
							 ! record density*thickness, heat*thickness
							 d_new(ii) = d_new(ii) + ice%snow%d(jj)*dz
							 heat_total = heat_total + ice%snow%heat(jj)*dz
							 melt_new(ii) = melt_new(ii) + ice%snow%melt(jj)*dz
			
							 ! ----- SETUP VARIABLES FOR NEXT PARTIAL LAYER ------------
							 ! keeping track of emerging new layer thickness, for 
							 dz_total = dz_total + dz 
			
							 ! find boundaries of next old layer
							 jj=jj+1
							 old_layer_top = old_layer_bot
							 old_layer_bot = old_layer_bot + ice%snow%th(jj)
						 enddo
			
						 ! add heat from 'new snow' to layer (using surface temp)
						 if (dz_total .lt. th_new(ii)) then
							 ! record heat*thickness
							 d_new(ii) = d_new(ii) + snowfall_d*(th_new(ii) - dz_total)
							 !T_mean = ice%snow%ts+kelvin0
							 !heat_mean = snowfall_d*(heat_snow0 - 0.2309*T_mean - &
							 !0.0034*T_mean**2)
							 heat_mean = sia2_snow_heat(snowfall_d,ice%snow%ts)
							 heat_total = heat_total + heat_mean*(th_new(ii) - dz_total)
							 melt_new(ii) = melt_new(ii) + &
								snowfall_melt*(th_new(ii) - dz_total)
							 !if (airtemp1c .gt. -1.0) then
							!	 melt_new(ii) = melt_new(ii) + (th_new(ii) - dz_total)
							 !endif
							 dz_total = th_new(ii)
							 
						 endif
					 else		           
						 ! add heat from 'new snow' to layer (from surface temp)
						 ! record heat*thickness
						 d_new(ii) = snowfall_d*th_new(ii)
						 !T_mean = ice%snow%ts+kelvin0
						 !heat_mean = snowfall_d*(heat_snow0 - 0.2309*T_mean - &
						 !0.0034*T_mean**2)
						 heat_mean = sia2_snow_heat(snowfall_d,ice%snow%ts)
						 heat_total = heat_total + heat_mean*th_new(ii)
						 melt_new(ii) = melt_new(ii) + snowfall_melt*th_new(ii)
						 !if (airtemp1c .gt. -1.0) then
						!	 melt_new(ii) = melt_new(ii) + th_new(ii)
						 !endif
						 dz_total = th_new(ii)
					 endif
			
				 endif  ! end of z_last .eq. 0 test
			
			 endif ! end of dz_total .lt. th_new test
			
			 ! Solve quadradic for new snow temp
			 d_new(ii) = d_new(ii)/dz_total
			 heat_total = heat_total/dz_total		 
			 !heat_total = heat_total/dz_total/d_new(ii)		 
			 !t_new(ii) = (0.2309-sqrt(0.05331481+0.0136* &
			 !	 (heat_snow0-heat_total)))/(-0.0068)
			 !t_new(ii) = min(t_new(ii) - kelvin0,c0)
			 t_new(ii) = sia2_snow_temp(d_new(ii),heat_total)
			 melt_new(ii) = melt_new(ii)/dz_total
			enddo
			
			! update snow thickness & temp
			ice%snow%z = z_snow_new
			do ii=1,z_snow_new
				ice%snow%th(ii) = th_new(ii)
				ice%snow%t(ii) = t_new(ii)
				ice%snow%d(ii) = d_new(ii)
				ice%snow%heat(ii) = sia2_snow_heat(d_new(ii),t_new(ii))
	!			 !T_mean = t_new(ii) + kelvin0
	!			 !d_mean = d_new(ii)					 
	!		#include "sia2_env_snow_calc.inc.f90"                  
			
	!			ice%snow%heat(ii) = heat_mean
				ice%snow%melt(ii) = melt_new(ii)
			enddo
			! zero out unused array members
			if (z_snow_new .lt. z_max_snow) then
				do ii=z_snow_new+1,z_max_snow
					ice%snow%th(ii) = c0
					ice%snow%t(ii) = c0
					ice%snow%d(ii) = c0					   
					ice%snow%heat(ii) = c0
					ice%snow%melt(ii) = c0
				enddo
			endif			   
		
		endif

	end subroutine sia2_snow_remap

! **********************************************************************


! **********************************************************************
! SUBROUTINE: sia2_melt_subl_ice
! melt/sublimate ice down from the surface, with heat flux F_surf
! ----------------------------------------------------------------------
	pure subroutine sia2_melt_subl_ice(ice,F_surf,h_top,h_bot,t_ocn, &
		ignore_f,basal,subl,dh_total,s_gl_sal)
		
	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only: &
			c0								! zero 
		use sia2_parameters, only: &
			w_bmelt,				& ! basal ice h2o flux -> ocn (kg/m^2)
			s_bmelt,				& ! basal ice salt flux -> ocn (kg/m^2)
			w_smelt,				& ! surface ice h2o flux -> ocn (kg/m^2)
			s_smelt						! surface ice salt flux -> ocn (kg/m^2)
		use sia2_types

		implicit none

	! Function Arguments
	! --------------------------------------------------------------------
		type (ice_type), intent(inout) :: &
			ice             ! ice structure
		real(kind=dbl_kind), intent(inout) :: &
			F_melt          ! heat forcing (W/m^2) > 0
		real(kind=dbl_kind), intent(in) :: &
			h_top,        &	! start of considered ice, in this direction specified by 'basal' (m)
			h_bot,				&	! end of considered ice, in this direction specified by 'basal' (m) (m)
			t_ocn						! ocean mixed layer temperature (degC)
		logical(kind=log_kind), intent(in) :: &
			ignore_f,     &	! switch in indicate whether to melt according to a flux (false) or not (true)
			basal,        &	! switch in indicate sublimation (true) or melt (false)
			subl            ! switch in indicate sublimation (true) or melt (false)
		real(kind=dbl_kind), intent(out) :: &
			dh_total,			&	! total ice thickness loss (m)
			s_gl_sal				! total salt from sublimed surface layers (m)
			
	! Internal Variables
	! --------------------------------------------------------------------
		integer(kind=int_kind) :: &
			wf,						&	! h2o flux type
			sf,						&	! salt flux type
			i,						&	! iterator
			incr					  ! iteration increment
		real(kind=dbl_kind) :: &
			dh_melt,			&	! temporary ice thickness (m)
			h_start,			&	! start thickness of ice eligible to melt (m)
			h_end,				&	! end thickness of ice eligible to melt (m)
			h1,						&	! current layer start (m)
			h2,						&	! current layer end (m)
			dz,						&	! thickness of layer eligible to melt (m)
			q_melt					! heat of melting (J/m^2)
			
	! Function Code
	! --------------------------------------------------------------------
			dh_total = c0
			if (basal .and. subl) then
				print *,'Error in sia2_melt_subl_ice: cannot sublime basal ice!'
			else
				if ((F_melt .gt. c0 .or. ignore_f) .and. ice%z .gt. 0) then

					! find parameters for determining melt eligibility of ice layers
					if (basal) then
						surf = ice%z
						incr = -1
						h_start = h_bot
						h_end = ice%th - h_top
					else
						surf = 1
						incr = 1
						h_start = h_top
						h_end = ice%th - h_bot
					endif
					i = surf
					h2 = c0

					do while ((F_melt .gt. c0 .or. ignore_f) .and. &
						(i .le. ice%z .and i .ge. 1))
						
						! find current layer bounds
						h1 = h2
						h2 = h2+ice%th(i)

						if (h1 .lt. h_end .and. h2 .gt. h_start) then

							! find layer thickness available to melt
							dz = min(ice%th(i),h2-hstart)
							if (h2 .gt. h_end) then
								dz = dz - (h2 - h_end)
							endif
								
							if (dz .gt. c0) then
	
								if (ignore_f) then
									dh_melt = dz
								else
	
									! determine heat of melting
									if (subl) then
										q_melt = ice%d(i)*(Lf + Lv) ! sublimation, no temperature adjustment
									else
										q_melt = sia2_ice_heat_melt(ice%t(i),ice%s(i),ice%d(i))
									endif
			
									! determine amount of layer that melts
									dh_melt = F_melt/q_melt						
									if (dh_melt .gt. dz) then
										F_melt = (dh_melt - dz)*q_melt ! adjust remaining melt flux
										dh_melt = dz
									else
										F_melt = c0
									endif

								endif
		
								! record mass fluxes from ice->ocn during melt, salt conservation during sublimation
								if (subl) then
									s_gl_sal = s_gl_sal + &
										sia2_layer_salt(ice%d(i),ice%s(i))*dh_melt    ! record kg salt for incorporation during regridding into top layer
								else
									if (basal)
										wf = w_bmelt
										sf = s_bmelt
									else
										wf = w_smelt
										sf = s_smelt
									endif
									ice%flux(wf) = ice%flux(wf) + &
										sia2_layer_h20(ice%d(i),ice%s(i))*dh_melt*ice%af  ! report kg melt water -> ocn
									ice%flux(sf) = ice%flux(sf) + &
										sia2_layer_salt(ice%d(i),ice%s(i))*dh_melt*ice%af  ! report kg salt -> ocn
								endif
								
								! record energy flux of melting if we are not already 
								! performing flux accounting with F_melt
								if (ignore_f) then
									ice%flux(q_latmelt) = ice%flux(q_latmelt) + &
										sia2_ice_heat_melt(ice%t(i),ice%s(i), &
										ice%d(i),t_ocn)*dh_melt*ice%af
								endif
		
								! add to total melt
								dh_total = dh_total + dh_melt
	
							endif ! end of dz > 0?
							
						endif ! end of is layer eligible for melt?

						! increment for next layer
						i = i+incr
					enddo
				endif	

	end subroutine sia2_melt_subl_ice

! **********************************************************************

! **********************************************************************
! SUBROUTINE: sia2_melt_subl_snow
! melt/sublimate ice down from the surface, with heat flux F_surf
! ----------------------------------------------------------------------
	pure subroutine sia2_melt_subl_snow(ice,F_surf,already_melted,t_ocn, &
		ignore_f,subl,dh_total)
		
	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only: &
			c0,							&	! zero 
			c_001							! 0.001
		use sia2_types

		implicit none

	! Function Arguments
	! --------------------------------------------------------------------
		type (ice_type), intent(inout) :: &
			ice             	! ice structure
		real(kind=dbl_kind), intent(inout) :: &
			F_surf          	! surface heat forcing (W/m^2) > 0
		real(kind=dbl_kind), intent(in) :: &
			already_melted,	&	! ! thickness of ice already melted in this direction by other processes (m)
			t_ocn							! ocean mixed layer temperature (degC)
		logical(kind=log_kind), intent(in) :: &
			ignore_f,     	&	! switch in indicate whether to melt according to a flux (false) or not (true)
			subl            	! switch in indicate sublimation (true) or melt (false)
		real(kind=dbl_kind), intent(out) :: &
			dh_total					! total surface thickness loss (m)
			
	! Internal Variables
	! --------------------------------------------------------------------
		integer(kind=int_kind) :: &
			i
		real(kind=dbl_kind) :: &
			ya_melt,					&	! running thickness of ice already melted (m)
			dz,								&	! thickness of layer eligible to melt (m)
			dh_melt,					&	! temporary ice thickness (m)
			q_melt						 	! heat of melting (J/m^2)
			
	! Function Code
	! --------------------------------------------------------------------
			dh_total = c0
			ya_melt = already_melted
			if (F_surf .gt. c0 .and. ice%snow%depth .ge. c0) then
				if(ice%snow%z .gt. c0) then
				! melt snow from regular snow grid
					i = ice%snow%z

					! find layer thickness available to melt
					if (ya_melt .gt. c0) then
						if(ya_melt .gt. ice%snow%th(i)) then
							ya_melt = ya_melt - ice%snow%th(i)
							dz = c0
						else
							dz = min(ya_melt,ice%snow%th(i))
						endif
					else
						dz = ice%snow%th(i)
					endif
					
					if (dz. gt. c0) then
					
						do while (((F_surf .gt. c0) .or. ignore_f) .and. i .ge. 1)
	
							! find height of snow melted
							if (ignore_f) then
								dh_melt = dz
							else
								if (subl) then
									q_melt = ice%snow%d(i)*(Lf + Lv)  ! sublimation heat
								else
									q_melt = sia2_snow_heat_melt(ice%snow%t(i), &
										ice%snow%d(i),t_ocn)
								endif							
								dh_melt = F_surf/q_melt						
								if (dh_melt .gt. dz) then
									F_surf = (dh_melt - dz)*q_melt ! adjust remaining melt flux
									dh_melt = dz
								else
									F_surf = c0
								endif
							endif

							! report snow - > ocn meltwater
							if (.not. subl) then
								ice%flux(w_smelt) = ice%flux(w_smelt) + &
									ice%snow%d(i)*dh_melt*c_001*ice%af    ! report kg melt water -> ocn
							endif
							if (ignore_f) then
								ice%flux(q_latmelt) = ice%flux(q_latmelt) + &
									snow(ice%snow%t(i),ice%d(i),t_ocn)*dh_melt*ice%af
							endif						

							dh_total = dh_total + dh_melt

							i = i-1
						enddo

					endif ! end is there snow that is not already melted?

				else
				! thin snow - melt snow small accumulation layer
				
					dz = max(c0, ice%snow%depth - ya_melt)
					if (dz .gt. c0) then
						
						! find height of snow melted
						q_melt = sia2_snow_heat_melt(ice%snow%ts, &
							ice%snow%d_small,t_ocn)
						dh_total = F_surf/q_melt
						if (dh_total .gt. dz) then
							F_surf = (dh_melt - dz)*q_melt ! adjust remaining melt flux
							dh_total = dz
						else
							F_surf = c0 ! no more F_surf to melt fs or ice
						endif
						
						! report snow - > ocn meltwater
						if (.not. subl) then
							ice%flux(w_smelt) = ice%flux(w_smelt) + &
								ice%snow%d_small*dh_total*c_001*ice%af     ! report kg melt water -> ocn
						endif

					endif ! end is there snow that in not already melted?

				endif ! end is there a snow grid or not?
			
			endif	

	end subroutine sia2_melt_subl_snow

! **********************************************************************


! **********************************************************************
! SUBROUTINE: sia2_destroy_ice
! wipe out ice, returning biomass, and salt, and h20 masses contained
! ----------------------------------------------------------------------
	subroutine sia2_destroy_ice(ice,m)

	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only: &
			c0,							&	! zero 
			c_001,					& ! 0.001
			c1e3,						& ! 1000
			c1e6							! 1,000,000
		use sia2_parameters, only: &
			cell_area					! area of grid cell in km^2
		use sia2_types

		implicit none
	
	! Function Arguments
	! --------------------------------------------------------------------
		type (ice_type), intent(inout) :: &
			ice             ! ice structure
		type (meta_type), intent(inout) :: &
			m             	! meta structure
			
	! Internal Variables
	! --------------------------------------------------------------------
		real(kind=dbl_kind) :: &
			smalg_lost,			&	! layer mass (kg/m^2)
			bv_mean,				&	! layer mass (kg/m^2)
			fh20,						& ! layer mass (kg/m^2)
			fsalt	            ! layer mass (kg/m^2)
			
	! Function Code
	! --------------------------------------------------------------------
		if (ice%af .gt. c0) then

			! record ice volume as model domain loss
			int_z = ice%z - z_sk
			m%md_loss = m%md_loss + &
				ice%id(ice%z)*ice%af*cell_area*c1e6

			! record biomass mass lost from ice / model domain							      
			do ii=1,ice%z
				if (ice(sc,ic,mi)%smalg(ii) .gt. c0) then
					bv_mean = ice%bv(ii)*c_001
					smalg_lost = (ice%smalg(ii)*bv_mean - min_alg)
					smalg_lost = max(c0,smalg_lost)
					! converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
					m%bm_lost = m%bm_lost + smalg_lost* &
						ice%th(ii)*cell_area*c1e3*ice%af   ! g/pixel
				endif
			enddo

			! record freshwater and salt fluxes
			call sia2_ice_mass(ice,fh20,fsalt)
			m%h2o_flux = m%h2o_flux + fh20*ice%af*cell_area*c1e6
			m%salt_flux = m%salt_flux + fsalt*ice%af*cell_area*c1e6
			call sia2_snow_mass(ice,fh20)
			m%h2o_flux = m%h2o_flux + fh20*ice%af*cell_area*c1e6
				
		endif

		! nullify ice structure
		call sia2_null_ice(ice)

		! nullify associated PUR history					
		!pur(:,:,:,sc,ic,mi) = pur_0

	end subroutine sia2_destroy_ice

! **********************************************************************

! **********************************************************************
! SUBROUTINE: sia2_null_ice
! zero ice structure
! ----------------------------------------------------------------------
  pure subroutine sia2_null_ice(ice)

	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only: &
			c0								! zero 
		use sia2_types	

	! Function Arguments
	! --------------------------------------------------------------------
		type (ice_type), intent(inout) :: &
			ice             	! ice structure

	! Function Code
	! --------------------------------------------------------------------

		ice%af = c0

		! null 2d snow params
		ice%snow_dist = c0
		ice%snow_rd = c0
		ice%sh_prev = c0
		ice%sh_offset = c0
		ice%snow%ts = c0
		ice%snow%depth = c0
		ice%snow%d_small = c0
		ice%snow%z = 0	
		
		! null 3d snow params
		ice%snow%t = c0
		ice%snow%th = c0
		ice%snow%d = c0
		ice%snow%heat = c0
		ice%snow%melt = c0
		ice%snow%ed_w = c0
		
		! null 2d ice params
		ice%Ed0 = c0
		ice%Ed1 = c0
		ice%PAR_bot = c0
		ice%PAR_bot_pl = c0
		ice%pur = c0
		ice%ed_w = c0
		ice%z = 0
		ice%age = c0
		ice%ridged = c0

		! null 3d ice params
		ice%smalg = c0
		ice%prod = c0
		ice%poc = c0
		ice%Ik1 = c0			  
		ice%no3 = c0
		ice%nh4 = c0
		ice%po4 = c0
		ice%sioh4 = c0
		ice%no3_mean = c0
		ice%nh4_mean = c0
		ice%po4_mean = c0
		ice%sioh4_mean = c0
		ice%th = c0
		ice%bv = c0
		ice%bd = c0
		ice%bs = c0
		ice%id = c0
		ice%d = c0
		ice%s = c0
		ice%t = c0
		ice%heat = c0
		ice%llim = c0
		ice%nlim = c0
		ice%plim = c0
		ice%silim = c0
		ice%slim = c0
		ice%bced = 0
		ice%drained = 0
		ice%dsdt = c0
		ice%fbv = c0
		ice%dhdt_conv = c0			  
		ice%f0 = c0			  
		ice%dsdt3 = c0			  
		ice%tgrad = c0			  

		! zero melt pond vars                
		ice%pond%t = c0
		ice%pond%s = c0		  
		ice%pond%d = c0
		ice%pond%heat = c0
		ice%pond%th = c0
		ice%pond%perc = c0
		ice%pond%smalg = c0     ! brine-based algal concentration (mgC/m^3)
		ice%pond%poc = c0     ! particulate organic carbon (detritus) (mgC/m^3)
		ice%pond%no3 = c0        ! ice layer brine NO3 concentration (µMol)
		ice%pond%nh4 = c0      ! ice layer brine NH4 concentration (µMol)
		ice%pond%po4 = c0     ! ice layer brine PO4 concentration (µMol)
		ice%pond%sioh4 = c0    ! ice layer brine SiOH4 concentration (µMol)                  

	end subroutine sia2_null_ice

! **********************************************************************


! **********************************************************************
! SUBROUTINE: sia2_interpolate_salinity
! used to find intital conditions from a 10-value distribution
! for 1st year ice salinity
! ----------------------------------------------------------------------

      PURE SUBROUTINE sia2_interpolate_salinity1(h_rel,value)
          use sia2_constants
          implicit none

          real(kind=dbl_kind),intent(in) :: h_rel
          real(kind=dbl_kind), intent(out) :: value
          
          !h_rel is the height of the desired measuremnt normalized to 1.
          !i.e., h_rel = interp_height/total_height
 
          ! interpolate between 10 points centered on integers -
          ! this means the below 0.05 an abosse 9.5 values are
          ! simply clipped to the end values
          !------------------------------------------------
          if (h_rel .lt. 0.05) then
              value=ss0
          elseif ((h_rel .ge. 0.05) .and. (h_rel .lt. 0.15)) then
              value = ss0+(ss1-ss0)*(h_rel-0.05)/0.1
          elseif ((h_rel .ge. 0.15) .and. (h_rel .lt. 0.25)) then
              value = ss1+(ss2-ss1)*(h_rel-0.15)/0.1
          elseif ((h_rel .ge. 0.25) .and. (h_rel .lt. 0.35)) then
              value = ss2+(ss3-ss2)*(h_rel-0.25)/0.1
          elseif ((h_rel .ge. 0.35) .and. (h_rel .lt. 0.45)) then
              value = ss3+(ss4-ss3)*(h_rel-0.35)/0.1
          elseif ((h_rel .ge. 0.45) .and. (h_rel .lt. 0.55)) then
              value = ss4+(ss5-ss4)*(h_rel-0.45)/0.1
          elseif ((h_rel .ge. 0.55) .and. (h_rel .lt. 0.65)) then
              value = ss5+(ss6-ss5)*(h_rel-0.55)/0.1
          elseif ((h_rel .ge. 0.65) .and. (h_rel .lt. 0.75)) then
              value = ss6+(ss7-ss6)*(h_rel-0.65)/0.1
          elseif ((h_rel .ge. 0.75) .and. (h_rel .lt. 0.85)) then
              value = ss7+(ss8-ss7)*(h_rel-0.75)/0.1
          elseif ((h_rel .ge. 0.85) .and. (h_rel .lt. 0.95)) then
              value = ss8+(ss9-ss8)*(h_rel-0.85)/0.1
          elseif ((h_rel .ge. 0.95) .and. (h_rel .le. 1.0)) then
              value = ss9
          endif
          
      end SUBROUTINE sia2_interpolate_salinity1

! **********************************************************************


! **********************************************************************
! SUBROUTINE: sia2_interpolate_salinity2
! used to find intital conditions from a 10-value distribution
! for multi-year ice salinity
! ----------------------------------------------------------------------

      PURE SUBROUTINE sia2_interpolate_salinity2(h_rel,value)
          use sia2_constants
          implicit none

          real(kind=dbl_kind),intent(in) :: h_rel
          real(kind=dbl_kind), intent(out) :: value
          
          !h_rel is the height of the desired measuremnt normalized to 1.
          !i.e., h_rel = interp_height/total_height
 
          ! interpolate between 10 points centered on integers -
          ! this means the below 0.05 an abosse 9.5 values are
          ! simply clipped to the end values
          !------------------------------------------------
          if (h_rel .lt. 0.05) then
              value=ssm0
          elseif ((h_rel .ge. 0.05) .and. (h_rel .lt. 0.15)) then
              value = ssm0+(ssm1-ssm0)*(h_rel-0.05)/0.1
          elseif ((h_rel .ge. 0.15) .and. (h_rel .lt. 0.25)) then
              value = ssm1+(ssm2-ssm1)*(h_rel-0.15)/0.1
          elseif ((h_rel .ge. 0.25) .and. (h_rel .lt. 0.35)) then
              value = ssm2+(ssm3-ssm2)*(h_rel-0.25)/0.1
          elseif ((h_rel .ge. 0.35) .and. (h_rel .lt. 0.45)) then
              value = ssm3+(ssm4-ssm3)*(h_rel-0.35)/0.1
          elseif ((h_rel .ge. 0.45) .and. (h_rel .lt. 0.55)) then
              value = ssm4+(ssm5-ssm4)*(h_rel-0.45)/0.1
          elseif ((h_rel .ge. 0.55) .and. (h_rel .lt. 0.65)) then
              value = ssm5+(ssm6-ssm5)*(h_rel-0.55)/0.1
          elseif ((h_rel .ge. 0.65) .and. (h_rel .lt. 0.75)) then
              value = ssm6+(ssm7-ssm6)*(h_rel-0.65)/0.1
          elseif ((h_rel .ge. 0.75) .and. (h_rel .lt. 0.85)) then
              value = ssm7+(ssm8-ssm7)*(h_rel-0.75)/0.1
          elseif ((h_rel .ge. 0.85) .and. (h_rel .lt. 0.95)) then
              value = ssm8+(ssm9-ssm8)*(h_rel-0.85)/0.1
          elseif ((h_rel .ge. 0.95) .and. (h_rel .le. 1.0)) then
              value = ssm9
          endif
                    
      end SUBROUTINE sia2_interpolate_salinity2

! **********************************************************************


end module sia2_grid