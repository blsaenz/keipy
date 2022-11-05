module sia2_desalination

	use sia2_constants, only : log_kind,int_kind,real_kind,dbl_kind
	
	implicit none

	public 

	contains

! **********************************************************************
! subroutine sia2_desal 
! Purpose: desalinate either using the SLDM (slush layer desalination 
! method - Saenz & Arrigo 2012) or Cox and Weeks (1988) gravity drainage


	subroutine sia2_desal(dtt_s,s_ocn,t_ocn,ice,s_new,ed_w_ice,dth,ki, &
		sldm,dsdt_out,lfheat)
	
 	! Use Statements and Variables/Globals/Parameters
	! --------------------------------------------------------------------
		use sia2_constants, only: &
			c0,					& ! zero
			c_5,				& !	0.5
			c1,					& !	1.0
			c_001,			& !	0.001
			Lf						! latent heat of fusion of ice
		use sia2_parameters, only: &
			z_max_ice, 	&	! maximum ice layers
			z_max_pack, &	! maximum ice + snow layers
			vb_crit, 		&	! critical volume of brine below which there is no convection
			bb_f,				& !	bubble fraction in ice
			bv_conv				! brine volume cutoff above which the SLDM is used
		use sia2_types
	
		implicit none
		
	! Subroutine Arguments
	! --------------------------------------------------------------------
		real(kind=dbl_kind), intent(in) :: &
			dtt_s,			&	! ice physics time step length (s)
			t_ocn,			&	! under-ice ocean temperature (degC)
			s_ocn					! under-ice ocean salinity (psu)
		type(ice_type), intent(inout) :: &
			ice						! saved ice data
		real(kind=dbl_kind), dimension(z_max_ice), intent(in) :: &
			s_new,			&	! updated mid-step ice salinity for use in salinity calc
			ed_w_ice			! shortwave absorption
		real(kind=dbl_kind), dimension(z_max_pack+1), intent(in) :: &
			dth,				& ! inter-layer thickness
			ki						! conductivity
		integer(kind=int_kind), dimension(z_max_ice), intent(out) :: &
			sldm					! was the SLDM (slush layer desalination method) active in layer?
		real(kind=dbl_kind), dimension(z_max_ice), intent(out) :: &
			dsdt_out,		&	! was the SLDM (slush layer desalination method) active in layer?
			lfheat				! effective heat latent heat flux due to desalination
		
	! Internal Variables
	! --------------------------------------------------------------------
		integer(kind=int_kind) :: &		
			ii,					& ! z dimension iterator, typically
			z_snow,			& ! number of snow layers
			int_z,			& ! bottom interior layer
			sk_1,				& ! top skeletal layer
			sk_z,				&	! bottom (skeletal) layer
			vb_open				! highest layer that has access to ocean through "open" ice, for cw gravity purposes
		real(kind=dbl_kind) :: &
			bv_total,		& ! total brine volume over multiple layers
			bv_flux,		& ! total brine flux volume over multiple layers
			dsdt1,			& ! dsdt holding var, cox/weeks gravity drainage
			dsdt2,			& ! dsdt holding var, cox/weeks brine expulsion
			dsdt3,			& ! dsdt holding var, slush desalination calc
			F0,					& ! Heat flux in/out of layer, + is cooling (W/m^2)
			bv_mean,		& ! factional brine volume
			T_Grad,			& ! layer temperature gradient
			T_melt,			& ! layer melting temp (c)
			dhdt,				& ! change in ice volume due to desalination (m)
			keff,				& ! salinity segregation coefficient, cox/weeks stable salinity
			tmp1,tmp2,tmp3,tmp4
					
	! Subroutine Code
	! --------------------------------------------------------------------
		!sk_z = ice%z
		!int_z = sk_z - z_sk
		!sk_1 = int_z + 1
		z_snow = ice%snow%z 
		int_z = ice%z


		dsdt_out = c0		! vector assignment
		lfheat = c0			! vector assignment
		sldm = 0				! vector assignment
		
		! find highest potentially convecting layer
 		vb_open = sia2_bv_open(ice%bv,int_z,vb_crit)

		do ii=1,int_z

			dsdt1 = c0
			dsdt2 = c0
			dsdt3 = c0
			F0 = c0
			bv_mean = ice%bv(ii)	! ppt

			! gravity drainage: only if above critical brine volume,
			! and density (salinity) of brine is higher than of surrounding seawater
			if ((bv_mean .ge. vb_crit) .and. ((ii .gt. vb_open) &
!								 .and. (ice%bs(ii) .gt. f%s)) &
!								 .or. ((ice%bced(ii) .eq. 1) .and. (bv_mean .ge. vb_crit)) &
			)) then

				if (ice%bs(ii) .gt. s_ocn) then

					! temp gradient used in gravity drainage eqn. is calculated
					! as the average between gradients above and below, then
					! averaged across this and last timestep

					! Cox-Weeks gravity drainage 
					! -----------------------------------------------
					if (ii .eq. 1) then
!						if (z_snow .gt. 0) then
!							T_grad = (ice%snow%t(z_snow) - ice%t(ii+1)) / &
!							(c_5*(ice%th(ii+1) + ice%th(ii) + ice%snow%t(z_snow)))
!						else														 
						T_grad = (ice%t(ii) - ice%t(ii+1)) / &
							(c_5*(ice%th(ii+1) + ice%th(ii)))
!						endif
					elseif (ii .eq. int_z) then
						T_grad = (ice%t(ii-1) - t_ocn) / &
							(c_5*ice%th(ii-1) + ice%th(ii))
					else
						T_grad = (ice%t(ii-1) - ice%t(ii+1)) / &
							(c_5*(ice%th(ii-1) + ice%th(ii+1)) + ice%th(ii))
					endif

					ice%tgrad(ii) = T_grad

					!dsdt1 = -1.*T_grad*c_01* & ! 0.01 factor is to convert from 1/m to 1/cm
					!(1.68E-5 - 3.37E-7*bv_mean)*dtt_s ! dsdt in ppt/s, so multiplying by sec., the minus sign gives a negative dsdt

					dsdt1 = 4.2e-6_dbl_kind*dtt_s*T_grad* &
						(max(bv_mean*c_001 - 0.05_dbl_kind,1.0e-5_dbl_kind))**(1.2) ! Petrich et al. desal rate

					! record instantaneous dsdt rate (psu/s) 
					!ice%dsdt(ii) = -1.*T_grad*c_01* & ! 0.01 factor is to convert from 1/m to 1/cm
					!	 (1.68E-5 - 3.37E-7*bv_mean)

					! find latent heat associated with cox/week gravity drainage
					! --------------------------------------------											
					if (dsdt1 .lt. c0) then
						! get volume of brine frozen, as % of total salinity * current brine volume * thickness * Lf
						! not worrying about density or brine salinity - this is simply the total fresh ice needed
						lfheat(ii) = abs(dsdt1)/s_new(ii)*bv_mean*c_001*ice%th(ii)*Lf*(c1-bb_f)	 
					endif

					! estimate conductive heat flux into layer
					! --------------------------------------------											
					if (ii .eq. 1) then
						if (z_snow .gt. 1) then
							tmp1 = (ice%t(1) - ice%snow%t(1))/dth(z_snow+1)
						else
							tmp1 = (ice%t(1) - ice%snow%ts)/dth(z_snow+1)
						endif
						tmp2 = (ice%t(1) - ice%t(2))/dth(z_snow+2)
					else
						tmp1 = (ice%t(ii) - ice%t(ii-1))/dth(z_snow+ii)
						if (ii .eq. ice%z) then
							tmp2 = (ice%t(ii) - t_ocn)/dth(z_snow+ii+1)
						else
							tmp2 = (ice%t(ii) - ice%t(ii+1))/dth(z_snow+ii+1)
						endif							
					endif
					F0 = tmp1*ki(z_snow+ii) + tmp2*ki(z_snow+ii+1) - ed_w_ice(ii)

					! F0 calc from temp gradient 
					! --------------------------------------------											
					tmp1 = F0/(c1-bb_f) ! overcompensate F0 to simulate bubble fraction
					dhdt = sia2_dhdt_from_flux(tmp1)
							
					! record instaneous flux-based desal vars
					ice%dhdt_conv(ii) = dhdt
					ice%f0(ii) = F0

					! pre-calculate SLDM desalination, store in dsdt3
					! --------------------------------------------											

					! assume F0 is used to freeze brine, with stable salinity type desalintion
					! calculated using keff
					
					if (dhdt .gt. c0) then
						! keff from cox/weeks
						!if (dhdt .gt. 3.6e-5) then
						!		 keff=0.26/(0.26+0.74*exp(-7243.0*dhdt))
						!else
						!		 if (dhdt .lt. 2.0e-6) then
						!				keff=0.12
						!		 else
						!				keff=0.8925+0.0568*log(dhdt)
						!		 endif
						!endif											

						! keff from petrich
						keff = sia2_keff(dhdt)
						!print *,'keffs: ',keff,tmp4


						! record convective desal rate
						! ice%dsdt3(ii) = -1.*dhdt*c_01*ice%bs(ii)*(1.-keff)/ice%th(ii)
						 ice%dsdt3(ii) = F0/Lf*ice%bs(ii)/IceD/ice%th(ii)

					endif

					! decide how to operate - if salinity is not stable,
					! hold temperature constant and do and convective desalination
					! otherwise use cox-weeks not-too-convective desal
!											 if (bv_mean .ge. bv_conv .and. dhdt .gt. 0.) then
					if (bv_mean .ge. bv_conv .and. dhdt .gt. c0) then

						if (s_new(ii) .gt. keff*ice%bs(ii)) then
						
							sldm(ii) = 1	! yes, using SLDM
							
							! use all heat flux to freeze ice
							lfheat(ii) = F0
							
							! Find SLDM desalination	 
							dsdt3 = -1.*dtt_s*F0/Lf*ice%bs(ii)/IceD/ice%th(ii)/(c1-bb_f)															

							! scale dsdt according to bv^3
!							tmp1 = (ice%bv(ii)-bv_conv)		! difference
!							tmp2 = (250.-bv_conv)				! scale
							dsdt3 = min(dsdt1,dsdt3)	 ! find bigger one
																				 ! assuming bv = 0.5 is wide open
!							if (tmp1 .lt. tmp2) then
!								dsdt3 = dsdt1 + (dsdt3-dsdt1)*(tmp1/tmp2)**3
!							endif
							
							! take larger amount of two desal values (more negative)
							dsdt1 = min(dsdt1,dsdt3)		

						endif

					endif													 

					! prevent increasing salinity
					dsdt1 = min(c1,dsdt1)

					! bound salinity at sal_min, don't gravity drain upward...
					if ((s_new(ii) + dsdt1) .lt. min_sal) then
						dsdt1 = min_sal - s_new(ii)
					endif

					! record desal for use later during environment re-calculation
					dsdt_out(ii) = dsdt1

				else
					
					! no gravity desal
					ice%dsdt(ii) = c0
					ice%dhdt_conv(ii) = c0
					ice%dsdt3(ii) = c0
					ice%f0(ii) = c0
					ice%tgrad(ii) = c0

				endif	 ! end of bs gradient check for gravity drainage

			endif ! end of check for open brine channels

		enddo		 ! end of desal calculation loop

	end subroutine sia2_desal
	
! **********************************************************************


! **********************************************************************
! function sia2_vb_open 
! Purpose: find higest layer that is open to brine convection (i.e. 
! brine volume greater than bv_crit)
! ----------------------------------------------------------------------
	
	integer(kind=int_kind) pure function sia2_bv_open(ice_bv,nz,bv_crit)
	
		implicit none
		
		! Function Arguments
		! --------------------------------------------------------------------
 		real(kind=dbl_kind), dimension(nz), intent(in) :: &
 			ice_bv  				! ice brine volume
 		integer(kind=int_kind), intent(in) :: &
 			nz			  				! number of ice layers
 		real(kind=dbl_kind), intent(in) :: &
 			bv_crit  				! critical brine volume below which there is no brine convection

		! Internal Variables
		! --------------------------------------------------------------------
 		integer(kind=int_kind) :: &
 			i
		
						
		! Function Code
		! --------------------------------------------------------------------
		sia2_bv_open = 0
	   do i=1,nz
			if (ice_bv(i) .lt. bv_crit) then
				sia2_bv_open = i
			endif
    enddo

	end function sia2_bv_open
	
! **********************************************************************


! **********************************************************************
! FUNCTION: sia2_dhdt_from_flux (Saenz and Arrigo 2012) 
! returns the rate of vertical sea ice formation from the total heat flux
! into the layer (cm/s), using a regression of the Petrich et al. (2006) 
! interpretation of the Nakawo/Signa 197? stable salinity vs ice growth
! rate data (Petrich calculation was a revision of one originally made 
! by Cox and Weeks (1988) 
! ----------------------------------------------------------------------	
	real(kind=dbl_kind) pure function sia2_dhdt_from_flux(F0)
	
		implicit none
	
	! Function Arguments
	! --------------------------------------------------------------------
		real(kind=dbl_kind), intent(in) :: &
			F0            	! Ice layer heat flux (W/m-2) extracted from layer (positive is freezing)
			
	! Function Code
	! --------------------------------------------------------------------
		if (F0 .ge. 210.0_dbl_kind) then
			sia2_dhdt_from_flux = 2.e-4_dbl_kind
		elseif (F0 .le. 57.316_dbl_kind) then
			sia2_dhdt_from_flux = 1.449e-9_dbl_kind*F0**2 &
				+ 3.524e-7_dbl_kind*F0 &
				+ 5.212e-8_dbl_kind
		else
			sia2_dhdt_from_flux = 1.942e-13_dbl_kind*F0**4 &
				- 6.e-11_dbl_kind*F0**3 &
				+ 7.387e-9_dbl_kind*F0**2 &
				+ 1.755e-7_dbl_kind*F0	! (cm/s) dhdt vs. F0 regression at T=-1.88 freezing temp, petrich keff
		endif
	  	  
	end function sia2_dhdt_from_flux

! **********************************************************************

! **********************************************************************
! FUNCTION: sia2_keff
! returns the 'stable salinity' of sea ice at a particular growth rate
! according to Petrich et al. 2006.
! ----------------------------------------------------------------------	
	real(kind=dbl_kind) pure function sia2_keff(dhdt_cm_s)
	
		implicit none
	
	! Function Arguments
	! --------------------------------------------------------------------
		real(kind=dbl_kind), intent(in) :: &
			dhdt_cm_s            	! rate of vertical sea ice formation (cm/s)
			
	! Function Code
	! --------------------------------------------------------------------
		sia2_keff = &
			0.19_dbl_kind*(dhdt_cm_s*7.4074e4_dbl_kind)**(0.46_dbl_kind)
		sia2_keff = max(sia2_keff,0.12_dbl_kind)
	 
		 	  
	end function sia2_keff

! **********************************************************************

end module sia2_desalination


