module sia2_grid

	use sia2_constants
	use sia2_parameters
	use sia2_state
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
	subroutine sia2_create_ice(ice,init_is,hi_new,hs_new,af_new,Tair, &
		t_ocn,s_ocn)
		
 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only: &
			c0,						& !
			c_5,					& !
			c1,						& !
			c1e6						!
		use sia2_parameters, only: &
			z_sk,					& !
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
			init_is 						!
		real(kind=dbl_kind), intent(in) :: &
			hi_new,	 					& !
			hs_new,	 					& !
			af_new,	 					& !
			Tair,	 						& !
			t_ocn,	 					& !
			s_ocn	 							!

	! Internal Variables
	! --------------------------------------------------------------------
		integer(kind=int_kind) :: &
			ii,									&	!
			int_z,							&	!
			sk_1,								&	!
			sk_z,								&	!
			z_new									!
		real(kind=dbl_kind) :: &
			tmp1,	 							& !
			h_temp,	 						& !
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
		t_interface = min(t_ocn,t_interface)
      
		! new snow density and 'melt' status
		if (Tair .lt. den_s_switch) then
			 snowfall_d = den_s_dry*c1e6
			 snowfall_melt = c0
		else
			 snowfall_d = den_s_wet*c1e6
			 snowfall_melt = c1
		endif
		
		! new snow
		if (hs_new .gt. c0) then
			! find new snow grid
			h_temp = hs_new
			call sia2_new_grid(h_temp,z_th_min,5.e0,2,1,z_max_snow,.false., &
				th_new,id_new,z_new)
			ice%snow%depth = hs_new
			ice%snow%th = th_new(1:z_max_snow)
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
		endif
		
		! find new ice grid			
		h_temp = hi_new
		call sia2_new_grid(h_temp,z_th_min,10.e0,2,1,z_max_ice-z_sk,.false., &
			th_new,id_new,z_new)
		sk_1 = z_new + 1
		sk_z = z_new + z_sk
		int_z = z_new
		!ice%z = sk_z
		ice%z = z_new

	 	! fill new ice parameters
	 	
		!do ii=1,sk_z
		do ii=1,z_new
		
!			if (ii .le. int_z) then

				! update thickness & ice depth
				ice%th(ii) = th_new(ii)
				ice%id(ii) = id_new(ii)
	
				! new ice temp & salinity
				if (init_is .eq. 0) then
					! constant ice salinity
					ice%s(ii) = s_const
					ice%t(ii) = t_ocn
				else				
					if (init_is .eq. 2) then
						! interpolate default new-ice salinity curve to get initial salinity
						call sia2_interpolate_salinity2( & 
						(ice%id(ii)-ice%th(ii)/c2)/hi_new,ice%s(ii))
					elseif (init_is .eq. 1) then
						! interpolate default new-ice salinity curve to get initial salinity
						call sia2_interpolate_salinity1( & 
						(ice%id(ii)-ice%th(ii)/c2)/hi_new,ice%s(ii))
					endif
						! interpolate linearly between t_interface and water temp
					ice%t(ii) = t_interface - (t_interface-t_ocn) * (ice%id(ii) - &
					ice%th(ii)/c2)/hi_new
				endif

!			else
!
!				! constant skeletal layer  - this will be ignored in heat/mass
!				! balance in KPP_ECO_ICE, but is required by regridding routines
!				ice%th(ii) = sk_th_max
!				ice%id(ii) = ice%id(ii-1) + sk_th_max
!			
!				! constant skeletal layer 
!				ice%s(ii) = s_ocn/c2
!				ice%t(ii) = t_ocn
!				
!			endif
									
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
		real(kind=dbl_kind) :: &
			layer_divisor,				&	!
			r_depth_dble,					&	!
			th_temp
		real(Kind=dbl_kind) &
			h_switch								!

		keep_growing = .true.         
		z_new = 0
		id_new = c0
		th_new = c0
		h_switch = z_max*h_min

		if (r_depth .gt. c0) then

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

			if ((grid_type .eq. 2) .and. (r_depth .le. h_switch)) then
			
				if (r_depth .lt. z_th_min*c4) then

					! determine ice layers if smaller than 4
					ii = int(r_depth/z_th_min)
					th_new(1:ii) = r_depth/ii
					z_new = ii

				else
				
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
				r_depth_dble = r_depth - dble(z_new)*z_th_min

				jj = z_max/2
				! find divisor by which we will divide depth bins
				! geometric series sum (minus 2 b/c top and bottom are constant)
				layer_divisor = 2.d0*(1.d0 - 1.3d0**jj)/(1.d0 - 1.3d0) - 2.d0 

				! find new thicknesses
				th_new(1) = z_th_min
				th_new(z_new) = z_th_min
				do ii=2,z_new/2
						jj = ii-1
						th_temp = z_th_min + (1.3d0**jj)/layer_divisor*r_depth_dble
						th_new(ii) = th_temp					
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
           id_new(1) = th_new(1)
           do ii=2,z_new
               id_new(ii) = id_new(ii-1) + th_new(ii)
           enddo
		endif
	
	end subroutine sia2_new_grid

! **********************************************************************


! **********************************************************************
! sia2_pre_grid: calculate grid surface options, new thickness, and new
! snow to be added if porous ice was drained above freeboard 
! ----------------------------------------------------------------------

	subroutine sia2_pre_grid(ice,ff,m,flooded,rain_fraction,c_gl,s_gl, &
	  snow_gl,r_depth,f_depth,did_flood,flood_1x,t_flood,s_flood_out, &
		melt_flood,melt_drain,did_drain,sn_depth,new_snow)

 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_state

		implicit none
		
	! Function Arguments
	! --------------------------------------------------------------------
		type(ice_type), intent(inout) :: &
			ice											!
		type(forcing_type), intent(in) :: &
			ff											!
		type(meta_type), intent(inout) :: &
			m												!
		real(kind=dbl_kind), intent(in) :: &
			flooded,							& !
			rain_fraction						!
		real(kind=dbl_kind), intent(inout) :: &
			c_gl,									& !
			s_gl,									&	!
			snow_gl									!
	real(kind=dbl_kind), intent(out) :: &
			r_depth,							& !
			f_depth,							& !
			t_flood,							& !
			s_flood_out, 							& !
			melt_flood,						& !
			melt_drain,						&	!
			sn_depth								!
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
			f_brine,							& !
			d_mean,								& !
			s_before_melt,				& !
			t_before_melt,				& !
			d_before_melt,				& !
			heat_mean								!

	! Subroutine Code
	! --------------------------------------------------------------------

	 	! find new thickness/icedepth structure with th_new and id_new. flooded will be added later
	 	int_z = ice%z !- z_sk
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
		if (flooded .gt. c0 .and. ice%snow%z .gt. 0) then
					  
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
			
			top_f = tmp1/IceD ! snow volume ratio					  
			bot_f = c1-top_f ! brine volume ratio

      s_before_melt = ff%s*(1.-rain_fraction)
			s_flood_out = max(min_sal,(s_before_melt*bot_f))
			d_before_melt = (ff%d*(1.-rain_fraction) + 1.e6*rain_fraction)
			t_before_melt = (ff%t*ff%d*(1.-rain_fraction) + &
			  max(ff%at-273.15,0.)*1.e6*rain_fraction) &
			  / d_before_melt
			
			! density of flooded later to be reconciled with new T,S,D
			d_mean = d_before_melt*bot_f + IceD*top_f

		 ! snow heat is assumed to be just latent heat of melting + 0.02*t*d
			heat_mean = cw*t_before_melt*d_before_melt*bot_f - &
				(ice%snow%t(z_snow)*0.02_dbl_kind + Lf)*IceD*top_f
			!dz_mean = 1.
			
			t_flood = sia2_ice_temp(s_flood_out,d_mean,heat_mean)

			!call sia2_ice_state(t_flood,s_flood_out,bs_flood,bd_flood,d_flood, &
			!	bv_flood,heat_flood		

			!ice%flux(w_flood) = ice%flux(w_flood) - bot_f*f_depth * c1000 * ice%af  ! report kg flood water from ocn
			!ice%flux(s_flood) = ice%flux(s_flood) - &
			!sia2_layer_salt(d_before_melt,s_before_melt) * &
			!ice%af*c_1e3 ! report kg flooded salt from ocn      


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
					     f_brine = ice%bv(jj)*c_001*ice%th(ii) ! m
							 melt_drain = melt_drain + ice%th(ii)
							 ! all salt goes, but only brine fraction goes (which is density of 1e3kg/m^3)
								!ice%flux(w_flood) = ice%flux(w_flood) + f_brine * c1000 * ice%af  ! report kg melt water drained -> ocn
								!ice%flux(s_flood) = ice%flux(s_flood) + &
								!	sia2_layer_salt(ice%d(ii),ice%s(ii)) * &
								!	ice%af*c_1e3 ! report kg salt drained -> ocn
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
    sn_depth = max(c0,ice%snow%depth + snow_gl + melt_drain)

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

	subroutine sia2_ice_remap(ice,ff,m,int_z_new,th_new,id_new,t_flood, &
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
			ff											!
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
		
	! Function Code
	! --------------------------------------------------------------------

		! assign local bounds variables
		sk_z = ice%z
		!int_z = sk_z - z_sk
		!sk_1 = int_z + 1
		
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
		new_layer_top = abs(s_gl)
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
!		if (.NOT. maxed_depth) then
				 
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
				 
			!do ii=1,int_z_new+z_sk
			do ii = 1,int_z_new
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
					tmp2 = tmp1/ff%s
					s_new(ii) = s_new(ii) + tmp1 
					d_new(ii) = d_new(ii) + d_flood*dz
					heat_total = heat_total + heat_flood*dz
					NO3_new(ii) = NO3_new(ii) + ff%no3*tmp2
					NH4_new(ii) = NH4_new(ii) + ff%nh4*tmp2
					PO4_new(ii) = PO4_new(ii) + ff%po4*tmp2
					SiOH4_new(ii) = SiOH4_new(ii) + ff%sioh4*tmp2
					poc_new(ii) = poc_new(ii) + ff%poc*tmp2
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
							
					!	if (ii .le. int_z_new) then
					!		if (jj .le. int_z .and. ii .le. int_z_new) then
							! redistribution of congelation ice to congelation ice (basal growth)
								tmp1 = dz
								s_mean = ice%s(jj)
								d_mean = ice%d(jj)
								heat_mean = ice%heat(jj)
					!		else
					!		! redistribution of skeletal ice into congelation ice (basal growth)
					!			tmp1 = min(1.,c_gl_sal/ice%s(jj))*dz ! remove brine volume
					!			t_mean = ff%t
					!			s_mean = c_gl_sal
					!			! find heat_mean
					!			call sia2_ice_state(t_mean,s_mean,bs_mean,bd_mean,d_mean, &
					!				bv_mean,heat_mean)       
					!		endif
					!	else
							! redistribution of congelation ice to skeletal ice (basal melt)
					!		tmp1 = dz
					!		t_mean = ff%t
					!		s_mean = ff%s/c2
					!		! find heat_mean
					!		call sia2_ice_state(t_mean,s_mean,bs_mean,bd_mean,d_mean, &
					!			bv_mean,heat_mean)       
					!	endif
																			
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
	
					!	if (ii .le. int_z_new .and. c_gl .gt. c0) then ! these tests should alway be true together
						! new congelation ice out of seawater - fast growth!!
							dz = max(c_gl,th_new(ii)-dz_total)
							t_mean = ff%t
							s_mean = c_gl_sal						 
					!	else
					!		! new skeletal ice out of seawater
					!		dz = th_new(ii) - dz_total
					!		t_mean = ff%t
					!		s_mean = ff%s/c2
					!	endif			 
						! fraction seawater concentrations by new brine volume
						call sia2_ice_state(t_mean,s_mean,bs_mean,bd_mean,d_mean, &
							 bv_mean,heat_mean)       					
						tmp1 = dz*bv_mean*c_001 
	
						! record tracers into new layer
						s_new(ii) = s_new(ii) + s_mean*dz 
						d_new(ii) = d_new(ii) + d_mean*dz
						heat_total = heat_total + heat_mean*dz
						no3_new(ii) = no3_new(ii) + ff%no3*tmp1 ! extending into mixed layer 
						nh4_new(ii) = nh4_new(ii) + ff%nh4*tmp1 ! extending into mixed layer
						po4_new(ii) = po4_new(ii) + ff%po4*tmp1 ! extending into mixed layer
						sioh4_new(ii) = sioh4_new(ii) + ff%sioh4*tmp1 ! extending into mixed layer
						smalg_new(ii) = smalg_new(ii) + alg_wc*tmp1
						poc_new(ii) = poc_new(ii) + ff%poc*tmp1
						
						! update boundaries for next layer					
						c_gl = c_gl - dz
						dz_total = th_new(ii) ! round out layer - this is the last way to add mass
						
					endif
	
					! round out dz_total whether we are shrinking or growing,
					! otherwise division below can create/take away tracers
					dz_total = th_new(ii)
	
				endif	 ! end of "if dz_total .lt. th_new(ii)								 
	
!				if (ii .le. int_z_new) then
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
					heat_mean = heat_total/dz_total
					if (s_mean .le. 0.) then
						print *,'Salinity (smean) is less than zero. Layer: ',ii
					endif								 
			
					t_mean = sia2_ice_temp(s_mean,d_mean,heat_mean)
		
					! check to make sure we don't get warmer than the ocean
					!if (ice%id(ii) .le. ice%fbh) then
					!		t_mean = min(t_mean,S_mean*mu)
					!else
					!		t_mean = min(t_mean,ff%t)
					!endif
		
					t_new(ii) = t_mean
				
!				else
					! we know skeletal temp, save computations
!					t_new(ii) = ff%t
!					s_new(ii) = s_new(ii)/dz_total
!					d_new(ii) = d_new(ii)/dz_total
!				endif
	
			enddo ! end of new grid portioning
	
			! update layer indices
			!int_z = int_z_new
			!sk_1 = int_z + 1
			!sk_z = int_z + z_sk
			sk_z = int_z_new
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
			
			if (sk_z .lt. z_max_ice) then
				sk_1 = sk_z+1
				
				ice%t(sk_1:z_max_ice) = c0
				ice%s(sk_1:z_max_ice) = c0
				ice%bs(sk_1:z_max_ice) = c0
				ice%bd(sk_1:z_max_ice) = c0
				ice%d(sk_1:z_max_ice) = c0
				ice%bv(sk_1:z_max_ice) = c0
				ice%heat(sk_1:z_max_ice) = c0
				ice%no3(sk_1:z_max_ice) = c0
				ice%nh4(sk_1:z_max_ice) = c0
				ice%po4(sk_1:z_max_ice) = c0
				ice%sioh4(sk_1:z_max_ice) = c0
				ice%poc(sk_1:z_max_ice) = c0
				ice%smalg(sk_1:z_max_ice) = c0
				
			endif
	
	!		dz_sk = c0
	!#include "sia2_env_update_skeletal.inc.f90"
	
			! adjust surface temp if there was no snow, but flooding occured (to prevent numerical instability)
			! uh, this can't happen...
			!if (flooded .gt. 0 .and. m(mi)%z_snow .eq. 0) then
			!		 ice%ts = ice%t(1)	! make equal to 1st layer temp
			!endif
							
		! endif	 ! end of "don't grow if depth will be maxed out by growth"
	
	end subroutine sia2_ice_remap

! **********************************************************************


! **********************************************************************

	! Adjust boundaries of snow layers and calculate new layer temperatures
	! using snow heat capacity to conserve total heat.  
	! ------------------------------------------------------------------

	subroutine sia2_snow_remap(ice,ff,new_snow,z_snow_new,th_new, &
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
			ff											!
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
			old_layer_top = c0
			new_layer_bot = c0
			old_layer_bot = c0
			d_new = c0
			melt_new = c0
	
			! new snow density and 'melt' status
			airtemp1c = ff%at - kelvin0
			if (airtemp1c .lt. den_s_switch) then
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
					 do while((new_layer_top .ge. old_layer_bot) .and. (jj .gt. 1))
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
									 if (jj .gt. 0) then
										 old_layer_top = old_layer_bot
										 old_layer_bot = old_layer_bot + ice%snow%th(jj)
									 endif
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
							 if (jj .le. z_last) then
							 	old_layer_top = old_layer_bot
							 	old_layer_bot = old_layer_bot + ice%snow%th(jj)
							 endif
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
			dz = c0
			do ii=1,z_snow_new
				dz = dz + th_new(ii)
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
			ice%snow%depth = dz
			
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
	pure subroutine sia2_melt_subl_ice(ice,F_melt,h_top,h_bot,t_ocn, &
		ignore_f,basal,subl,dh_total,s_gl_sal)
		
	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only: &
			c0,							&	! zero 
			c_1e3							! 0.001 
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
			surf,					&	! 'top' surface layer w/ respect to melt
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
		if ( .not. (basal .and. subl)) then

			if ((F_melt .gt. c0 .or. ignore_f) .and. ice%z .gt. 0) then

				! find parameters for determining melt eligibility of ice layers
				if (basal) then
					surf = ice%z
					incr = -1
					h_start = h_bot
					h_end = ice%id(ice%z) - h_top
				else
					surf = 1
					incr = 1
					h_start = h_top
					h_end = ice%id(ice%z) - h_bot
				endif
				i = surf
				h2 = c0

				do while ((F_melt .gt. c0 .or. ignore_f) .and. &
					(i .le. ice%z .and. i .ge. 1))
					
					! find current layer bounds
					h1 = h2
					h2 = h2+ice%th(i)

					if (h1 .lt. h_end .and. h2 .gt. h_start) then

						! find layer thickness available to melt
						dz = min(ice%th(i),h2-h_start)
						if (h2 .gt. h_end) then
							dz = dz - (h2 - h_end)
						endif
							
						if (dz .gt. c0) then

							! determine heat of melting
							if (subl) then
								q_melt = ice%d(i)*(Lf + Lv) ! sublimation, no temperature adjustment
							else
								q_melt = sia2_ice_heat_melt(ice%t(i),ice%s(i), &
									ice%d(i),t_ocn)
							endif

							! determine amount of layer that melts
							if (ignore_f) then
								dh_melt = dz
							else		
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
								if (basal) then
									wf = w_bmelt
									sf = s_bmelt
								else
									wf = w_smelt
									sf = s_smelt
								endif
								ice%flux(wf) = ice%flux(wf) + &
									sia2_layer_h2o(ice%d(i),ice%s(i)) * &
									dh_melt*ice%af*c_1e3  ! report kg melt water -> ocn
								ice%flux(sf) = ice%flux(sf) + &
									sia2_layer_salt(ice%d(i),ice%s(i)) * &
									dh_melt*ice%af*c_1e3  ! report kg salt -> ocn
							endif
							
							! record energy flux of melting if we are not already 
							! performing flux accounting with F_melt
							if (ignore_f) then
								ice%flux(q_latmelt) = ice%flux(q_latmelt) + &
									q_melt*dh_melt*ice%af
							endif
	
							! add to total melt
							dh_total = dh_total + dh_melt

						endif ! end of dz > 0?
						
					endif ! end of is layer eligible for melt?

					! increment for next layer
					i = i+incr
				enddo
			endif
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
					
					if (dz .gt. c0) then
					
						do while (((F_surf .gt. c0) .or. ignore_f) .and. i .ge. 1)
	
							! find height of snow melted
							if (subl) then
								q_melt = ice%snow%d(i)*(Lf + Lv)  ! sublimation heat
							else
								q_melt = sia2_snow_heat_melt(ice%snow%t(i), &
									ice%snow%d(i),t_ocn)
							endif							
							if (ignore_f) then
								dh_melt = dz
							else
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
								ice%flux(w_snowmelt) = ice%flux(w_snowmelt) + &
									ice%snow%d(i)*dh_melt*c_001*ice%af    ! report kg melt water -> ocn
							endif
							if (ignore_f) then
								ice%flux(q_latmelt) = ice%flux(q_latmelt) + &
									q_melt*dh_melt*ice%af
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
							F_surf = (dh_total - dz)*q_melt ! adjust remaining melt flux
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
! FUNCTION: sia2_basal_growth
! returns the thickness and salinity of basal ice growth from a flux
! and a temperature
! ----------------------------------------------------------------------	
	pure subroutine sia2_basal_growth(ice,t_ocn,s_ocn,Fb,dtt_s,dhdt,dhdt_sal, &
		dhdt_cm_s)
	
 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only: &
			c_1e3,						&	! 0.001 	
			c100								! 100 	
		use sia2_desalination, only: sia2_dhdt_from_flux,sia2_keff
		use sia2_state, only: sia2_bs,sia2_bd,sia2_ice_d, &
			sia2_ice_heat_melt,sia2_layer_h2o,sia2_layer_salt
	
		implicit none
	
	! Subroutine Arguments
	! --------------------------------------------------------------------
		type (ice_type), intent(inout) :: &
			ice             	! ice structure
		real(kind=dbl_kind), intent(in) :: &
			t_ocn,            &	! underlying ocean/mixed layer temperature (degC)
			s_ocn,            &	! underlying ocean/mixed layer salinity (psu)
			Fb,            		&	! basal heat flux (W/m^2)
			dtt_s	            	! length of time step (s)
		real(kind=dbl_kind), intent(out) :: &
			dhdt,            	&	! thickness of new basal growth layer (m)
			dhdt_sal,	        &	! predicted salinity of new basal growth layer (psu)
			dhdt_cm_s						! growth rate (cm/s)

	! Internal Variables
	! --------------------------------------------------------------------
		real(kind=dbl_kind) :: &
			bs,								& !	brine salinity (psu)
			bd,								& !	brine denisty (g/m^3)
			d_new,						&	! bulk ice density (g/m^3)
			q_melt							! heat of freezing
			
	! Subroutine Code
	! --------------------------------------------------------------------
		dhdt_cm_s = sia2_dhdt_from_flux(Fb)		! predicted ice growth rate
		dhdt_sal = s_ocn*sia2_keff(dhdt_cm_s)          	! predicted new ice salinity
		bs = sia2_bs(t_ocn,dhdt_sal)                   		! new ice brine salinity
		bd = sia2_bd(bs)																! new ice brine density
		d_new = sia2_ice_d(t_ocn,dhdt_sal,bs,bd)						! new ice density
		q_melt = sia2_ice_heat_melt(t_ocn,dhdt_sal,d_new,t_ocn)	! new ice melting heat
		dhdt = Fb/q_melt															! m/s
		dhdt_cm_s = dhdt*c100													! report actual dhdt_cm_s, not predicted
		dhdt = dhdt*dtt_s															! new ice growth, using predicted salinity
		
		! update fresh water and salt fluxes (ice->ocn)		
		ice%flux(w_bfreeze) = ice%flux(w_bfreeze) - & 
			sia2_layer_h2o(d_new,dhdt_sal)*c_1e3*dhdt*ice%af ! kg m-2
		ice%flux(s_bfreeze) = ice%flux(s_bfreeze) - &
			sia2_layer_salt(d_new,dhdt_sal)*c_1e3*dhdt*ice%af ! kg m-2		
	 		 	  
	end subroutine sia2_basal_growth

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
			c1e6,						&	! 1,000,000
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
		integer(kind=int_kind) :: &
			int_z,          &	! last internal ice layer
			ii           			! iterator
		real(kind=dbl_kind) :: &
			smalg_lost,			&	! layer mass (kg/m^2)
			bv_mean,				&	! layer mass (kg/m^2)
			fh20,						& ! layer mass (kg/m^2)
			fsalt	            ! layer mass (kg/m^2)
			
	! Function Code
	! --------------------------------------------------------------------
		if (ice%af .gt. c0) then

			! record ice volume as model domain loss
			!int_z = ice%z - z_sk
			m%md_loss = m%md_loss + &
				ice%id(ice%z)*ice%af*cell_area*c1e6

			! record biomass mass lost from ice / model domain							      
			do ii=1,ice%z
				if (ice%smalg(ii) .gt. c0) then
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



! SUBROUTINE: sia2_env_merge_ice
! ======================================================================
! merges ice conservatively from two categories

    SUBROUTINE sia2_merge_ice(ice_in,ice_out,t_ocn,s_ocn)

	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
      USE sia2_constants
      USE sia2_state
      IMPLICIT NONE

    ! FUNCTION ARGUMENTS
    !----------------------------------------------------          
        
      TYPE (ice_type), INTENT(IN) :: ice_in
      TYPE (ice_type), INTENT(INOUT) :: ice_out          
      REAL(kind=dbl_kind), INTENT(IN) :: t_ocn,s_ocn

    ! INTERNAL VARIABLES
    !----------------------------------------------------                    
      integer(kind=int_kind) :: ii,jj,kk,jjj,iii,int_z,sk_1,sk_z,z_old
      real(kind=dbl_kind) :: t_mean,s_mean,d_mean,heat_mean, &
        tmp1,tmp2,heat_total,af_total,vf_total,z1,z2
      real(kind=dbl_kind), dimension(z_max_ice) :: t_new,s_new,d_new,th_new,id_new
      real(kind=dbl_kind), dimension(z_max_ice) :: smalg_new,poc_new,heat_new
      real(kind=dbl_kind), dimension(z_max_ice+1) :: NO3_new,NH4_new,PO4_new,SiOH4_new

    ! SUBROUTINE CODE
    ! ------------------------------------------------------------

      ! initializations
      !int_z = max(ice_in%z-z_sk,ice_out%z-z_sk)
      int_z = max(ice_in%z,ice_out%z)
      th_new = c0  ! vector assignment
      id_new = c0  ! vector assignment

      af_total = ice_in%af + ice_out%af ! total area fraction we are combining

      do ii=1,int_z
          
          ! calculate fractional mutlipliers for both categories                  
          !if (ii .le. ice_out%z-z_sk) then
          if (ii .le. ice_out%z) then
              z1 = ice_out%th(ii)*ice_out%af
          else
              z1 = c0
          endif
          
          z2 = c0
          if (ice_in%af .gt. c0) then
              !if (ii .le. ice_in%z-z_sk) then
              if (ii .le. ice_in%z) then
                  z2 = ice_in%th(ii)*ice_in%af
              endif
          endif                  

          ! find new bounds
          th_new(ii) = (z1+z2)/af_total
          if (ii .eq. 1) then
              id_new(ii) = th_new(ii)
          else
              id_new(ii) = id_new(ii-1) + th_new(ii)
          endif
          vf_total = af_total*th_new(ii)
          
          ! find new temp, salinity
          !dz_mean = c1
          heat_mean = (ice_out%heat(ii)*z1 + ice_in%heat(ii)*z2)/vf_total
          d_mean = (ice_out%d(ii)*z1 + ice_in%d(ii)*z2)/vf_total
          s_new(ii) = (ice_out%s(ii)*z1 + ice_in%s(ii)*z2)/vf_total
          t_new(ii) = sia2_ice_temp(s_new(ii),d_mean,heat_mean)
           
          ! check to make sure we don't get warmer than the ocean
          !if (ice_out%id(ii) .le. ice_out%fbh) then
          !   t_mean = min(t_mean,S_mean*mu)
          !else
          !   t_mean = min(t_mean,f(mi)%t)
          !endif
          
          ! record tracers
          smalg_new(ii) = (ice_out%smalg(ii)*z1 + ice_in%smalg(ii)*z2)/vf_total
          poc_new(ii) = (ice_out%poc(ii)*z1 + ice_in%poc(ii)*z2)/vf_total
          no3_new(ii) = (ice_out%no3(ii)*z1 + ice_in%no3(ii)*z2)/vf_total
          nh4_new(ii) = (ice_out%nh4(ii)*z1 + ice_in%nh4(ii)*z2)/vf_total
          po4_new(ii) = (ice_out%po4(ii)*z1 + ice_in%po4(ii)*z2)/vf_total
          sioh4_new(ii) = (ice_out%sioh4(ii)*z1 + ice_in%sioh4(ii)*z2)/vf_total

      enddo
      
      ! redistribute skeletal layers
      ! ------------------------------------------------------------
!      do ii=1,z_sk
!
!          ! find corresponding skeletal layers
!          jj = ice_out%z-z_sk+ii ! ic skeletal layer
!          if (ice_in%af .gt. 0.) then
!              jjj = ice_in%z-z_sk+ii ! ic_n skeletal layer
!          else
!              jjj=1 ! no ice in category ic
!          endif
!          kk = int_z+ii ! new skeletal layer
!
!          ! calculate fractional mutlipliers for both categories                  
!          z1 = ice_out%th(jj)*ice_out%af
!          z2 = ice_out%th(jjj)*ice_in%af
!
!          ! find new bounds
!          th_new(kk) = (z1+z2)/af_total
!          id_new(kk) = id_new(kk-1) + th_new(kk)
!          vf_total = af_total*th_new(kk)
!          
!          ! record tracers
!          smalg_new(kk) = (ice_out%smalg(jj)*z1 + ice_in%smalg(jjj)*z2)/vf_total
!          poc_new(kk) = (ice_out%poc(jj)*z1 + ice_in%poc(jjj)*z2)/vf_total
!          no3_new(kk) = (ice_out%no3(jj)*z1 + ice_in%no3(jjj)*z2)/vf_total
!          nh4_new(kk) = (ice_out%nh4(jj)*z1 + ice_in%nh4(jjj)*z2)/vf_total
!          po4_new(kk) = (ice_out%po4(jj)*z1 + ice_in%po4(jjj)*z2)/vf_total
!          sioh4_new(kk) = (ice_out%sioh4(jj)*z1 + ice_in%sioh4(jjj)*z2)/vf_total
!
!          ! keep skeletal physics constant
!          s_new(kk) = s_ocn*c_5
!          t_new(kk) = min(t_ocn,s_new(kk)*mu)
!
!      enddo

      ! update surface temperature (not heat conservative yet)
      ! ------------------------------------------------------------
      if (ice_out%snow%z .eq. 0) then
          if (ice_in%snow%z .eq. 0) then
              ! average ts if no snow - not heat conservative!!
              ice_out%snow%ts = (ice_out%snow%ts*ice_out%af &
              + ice_in%snow%ts*ice_in%af)/af_total
          else
              ! set ts equal to ts of layer w/ snow
              ice_out%snow%ts = ice_out%snow%ts
          endif
      elseif (ice_in%snow%z .gt. 0) then
          ! both have snow - take areal average (not heat conservative)!!
          ice_out%snow%ts = (ice_out%snow%ts*ice_out%af &
          + ice_in%snow%ts*ice_in%af)/af_total
      endif                      

      ! update ice physics, carry-over parameters, boundaries
      ! ------------------------------------------------------------
      z_old = ice_out%z
      !ice_out%z = int_z + z_sk
      ice_out%z = int_z 

      ! update ice structure
      do ii=1,ice_out%z

        ! assign new ice physics
        ice_out%t(ii) = t_new(ii)
        ice_out%s(ii) = s_new(ii)
        call sia2_ice_state(t_new(ii),s_new(ii),ice_out%bs(ii),ice_out%bd(ii), &
          ice_out%d(ii),ice_out%bv(ii),ice_out%heat(ii))
    
        ! update tracers
        ice_out%no3(ii) = NO3_new(ii)
        ice_out%nh4(ii) = NH4_new(ii)
        ice_out%po4(ii) = PO4_new(ii)
        ice_out%sioh4(ii) = SiOH4_new(ii)
        ice_out%poc(ii) = poc_new(ii)
        ice_out%smalg(ii) = smalg_new(ii)

          ! transfer carry-over paramaters
          if (ii .le. z_old) then
              z1 = ice_out%th(ii)*ice_out%af
          else
              z1 = 0.
          endif
          if (ii .le. ice_in%z) then
              z2 = ice_in%th(ii)*ice_in%af
          else
              z2 = 0.
          endif                  
          vf_total = af_total*th_new(ii)
          
          ! vars with units x/pixel
          ice_out%prod(ii) = ice_in%prod(ii) + ice_out%prod(ii)

          ! vars with units x/m^2
          ice_out%fbv(ii) = (ice_in%fbv(ii)*ice_in%af &
              + ice_out%fbv(ii)*ice_out%af)/af_total

          ! vars with units of x/m^3 and descriptive stats
!                      ice(sc,ic,mi)%dhdt_conv(ii) = (ice(sc,ic,mi)%dhdt_conv(ii)*z2 &
!                          + ice(sc_out,ic_out,mi)%dhdt_conv(ii)*z1)/vf_total
!                      ice(sc,ic,mi)%f0(ii) = (ice(sc,ic,mi)%f0(ii)*z2 &
!                          + ice(sc_out,ic_out,mi)%f0(ii)*z1)/vf_total
!                      do jj=1,dt_per_day
!                          do kk=1,lda_n
!                              pur(ii,jj,kk,ic,mi) = nint((pur(ii,jj,kk,ic,mi)*z2 &
!                                  + pur(ii,jj,kk,ic_out,mi)*z1)/vf_total)
!                          enddo
!                      enddo
                  
          ! instantaneous/snapshot parameters that we will ignore
          !ice(sc,ic,mi)%bced(ii) =
          !ice(sc,ic,mi)%dsdt(ii) = 0.
          !ice(sc,ic,mi)%Ik1(ii) = 0.
          !ice(sc,ic,mi)%llim(ii) = 0.
          !ice(sc,ic,mi)%nlim(ii) = 0.
          !ice(sc,ic,mi)%plim(ii) = 0.
          !ice(sc,ic,mi)%silim(ii) = 0.
          !ice(sc,ic,mi)%slim(ii) = 0.
          !ice(sc,ic,mi)%dsdt3(ii) = 0.       
          !ice(sc,ic,mi)%tgrad(ii) = 0.       
          

          ! dimensions
          ice_out%th(ii) = th_new(ii)
          ice_out%id(ii) = id_new(ii)
          
      enddo          

! remapping now done outside of merge
!          ! update these values for remapping
!          sk_z = int_z + z_sk
!          sk_1 = int_z + 1
!
!          r_depth = ice(sc,ic,mi)%id(int_z)
!
!#include "sia2_env_ice_grid.inc.f90"               
!    
!          ! start new layer top at old layer top - don't exclude any ice
!          new_layer_top = -1.
!          ! no growth/melt here
!          s_gl = 0.
!          c_gl = 0.
!          ! no flooding here
!          flooded = 0.
!          melt_flood = 0.
!          dhdt_cm_s = 0.   ! required in algal migration
!          maxed_depth = .false.
!        
!#include "sia2_env_ice_remap.inc.f90"


      ! transfer pond paramaters
      ! --------------------------------------------------------          
      if (use_ponds .eq. 1) then
          tmp1 = ice_in%af*ice_in%pond%th
          tmp2 = ice_out%af*ice_out%pond%th
          vf_total = tmp1 + tmp2
      
          if (vf_total .gt. c0) then

              ice_out%pond%s = (ice_in%pond%s*tmp1 &
                  + ice_out%pond%s*tmp2)/vf_total     
              ice_out%pond%d = (ice_in%pond%d*tmp1 &
                  + ice_out%pond%d*tmp2)/vf_total
              ice_out%pond%heat = (ice_in%pond%heat*tmp1 &
                  + ice_out%pond%heat*tmp2)/vf_total
              ice_out%pond%smalg = (ice_in%pond%smalg*tmp1 &
                  + ice_out%pond%smalg*tmp2)/vf_total     ! brine-based algal concentration (mgC/m^3)
              ice_out%pond%poc = (ice_in%pond%poc*tmp1 &
                  + ice_out%pond%poc*tmp2)/vf_total     ! particulate organic carbon (detritus) (mgC/m^3)
              ice_out%pond%no3 = (ice_in%pond%no3*tmp1 &
                  + ice_out%pond%no3*tmp2)/vf_total        ! ice layer brine NO3 concentration (µMol)
              ice_out%pond%nh4 = (ice_in%pond%nh4*tmp1 &
                  + ice_out%pond%nh4*tmp2)/vf_total      ! ice layer brine NH4 concentration (µMol)
              ice_out%pond%po4 = (ice_in%pond%po4*tmp1 &
                  + ice_out%pond%po4*tmp2)/vf_total     ! ice layer brine PO4 concentration (µMol)
              ice_out%pond%sioh4 = (ice_in%pond%sioh4*tmp1 &
                  + ice_out%pond%sioh4*tmp2)/vf_total    ! ice layer brine SiOH4 concentration (µMol)                  

              ice_out%pond%th = (tmp1+tmp2)/af_total
              ice_out%pond%perc = (ice_in%af*ice_in%pond%perc + &
                ice_out%af*ice_out%pond%perc)/af_total

          endif
      endif

      ! redistribute snow 
      ! ------------------------------------------------------------
      ! calculate fractional mutlipliers for both categories
      z1 = ice_out%af/(ice_out%af + ice_in%af)
      z2 = ice_in%af/(ice_out%af + ice_in%af)

      ice_out%sh_prev = ice_out%sh_prev*z1 + ice_in%sh_prev*z2
      ice_out%sh_offset = ice_out%sh_offset*z1 + ice_in%sh_offset*z2
      

      jjj = max(ice_out%snow%z,ice_in%snow%z)
      ice_out%snow%depth = 0.
      do ii=1,jjj
          
          ! calculate fractional mutlipliers for both categories                  
          z1 = ice_out%snow%th(ii)*ice_out%af
          z2 = ice_in%snow%th(ii)*ice_in%af

          ! find new bounds
          th_new(ii) = (z1+z2)/af_total
          if (ii .eq. 1) then
              id_new(ii) = th_new(ii)
          else
              id_new(ii) = id_new(ii-1) + th_new(ii)
          endif
          vf_total = af_total*th_new(ii)
          
          ! find new temp
          heat_total = (ice_out%snow%heat(ii)*z1 + &
              ice_in%snow%heat(ii)*z2)/vf_total
          ice_out%snow%d(ii) = (ice_out%snow%d(ii)*z1 + &
              ice_in%snow%d(ii)*z2)/vf_total
          ice_out%snow%t(ii) = sia2_snow_temp(ice_out%snow%d(ii),heat_total)
          ice_out%snow%heat(ii) = sia2_snow_heat(ice_out%snow%d(ii),ice_out%snow%t(ii))
           
          ! update all other snow parameters
          ice_out%snow%th(ii) = th_new(ii)
          ice_out%snow%melt(ii) = (ice_out%snow%melt(ii)*z1 + &
              ice_in%snow%melt(ii)*z2)/vf_total

          ! find new total snow depth
          ice_out%snow%depth = ice_out%snow%depth + th_new(ii)

      enddo

      ! record final number snow layers
      ice_out%snow%z = jjj
!
!          ! setup for snow regridding 
!          r_depth = ice(sc,ic,mi)%snow%depth
!          ice(sc,ic,mi)%sh_prev = r_depth                  
!
!#include "sia2_env_snow_grid.inc.f90"               
!
!          flooded = 0.
!          melt_drain = 0.
!
!#include "sia2_env_snow_remap.inc.f90"                                 

      tmp1 = ice_in%af + ice_out%af

      ! update age and ridged tracers
      ! --------------------------------------------------------
      ice_out%age = (ice_in%age*ice_in%af + ice_out%age*ice_out%af)/tmp1
      ice_out%ridged = (ice_in%ridged*ice_in%af + ice_out%ridged*ice_out%af)/tmp1

      ! update ml layer tracers
      ! --------------------------------------------------------
!      ice_out%no3(ml_z) = (ice_out%no3(ml_z)*ice_out%af + &
!        ice_in%nh4(ml_z)*ice_in%af)/tmp1
!      ice_out%nh4(ml_z) = (ice_out%nh4(ml_z)*ice_out%af + &
!        ice_in%nh4(ml_z)*ice_in%af)/tmp1
!      ice_out%po4(ml_z) = (ice_out%po4(ml_z)*ice_out%af + &
!        ice_in%po4(ml_z)*ice_in%af)/tmp1
!      ice_out%sioh4(ml_z) = (ice_out%sioh4(ml_z)*ice_out%af + &
!        ice_in%sioh4(ml_z)*ice_in%af)/tmp1
!         ice_in%poc(ml_z) = (ice_out%poc(ml_z)*ice_out%af + &
!           ice_in%poc(ml_z)*ice_in%af)/tmp1

      ! update af & nullify old category
      ! --------------------------------------------------------
      ice_out%af = ice_in%af + ice_out%af                  
      ! call sia2_null_ice(ice_in)

    end subroutine sia2_merge_ice


end module sia2_grid