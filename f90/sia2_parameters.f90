! ======================================================================
! Global Parameters
! ======================================================================

module sia2_parameters

	use sia2_constants, only: log_kind,int_kind,real_kind,dbl_kind
	
	implicit none

	public


	! ------------------------------------------------------------
	! run options
	! ------------------------------------------------------------
	integer(kind=int_kind), parameter :: &
		rstart 				=	0	    ,	&	!	1) use retart files to resume simulation 2) start new simulation	
		adv_on 				=	1	    ,	&	!	advect ice according to ice motion product	
		use_gpu  			=	1	    ,	&	!	perform delta-eddington light calculations using GPU card	
		stn_only 			=	0	    ,	&	!	only compute stations, ignoring model domain (mdh1,mdv1, etc.) below	
		atm_f 				=	2	    ,	&	!	atmospheric forcing source (0=NCEP,1=ECMWF)	
		atmo 					=	1	    ,	&	!	atmospheric boundary calculation type (1-CICE 0=muntersteiner old/kottmeier 2003)	
		woa 					=	1	    ,	&	!	toggle for use of specified forcing (0) or world ocean atlas forcing (1)	
		ohf 					=	2	    ,	&	!	ocean heat flux is specified (0) or derived from world ocean atlas (1) or generated from internal climatology (2)	
		validation 		=	0	    ,	&	!	toggle for modeling of (1) validation station, (0) or specified model domain	
		monte_carlo 	=	0	    ,	&	!	toggle for running purmutations of validation based on modified snowfall,temperature	
		override_ic 	=	0	    ,	&	!	override the built-in dependence on satellite ice conectration, assuming ice the whole run	
		use_ponds 		=	1	    ,	&	!	1=don't disappear surface melted ice, 0=melted ice/snow dissappears	
		use_drain 		=	1	    ,	&	!	1=drain surface ice if above bv_conv and turn into snow, 0=don't	
		mdh1 					=	202	  , & !	model domain - horizontal upper left (RossSeaBox = 372/744, Full = 202/404)  400  300	
		mdv1 					=	194	  , & !	model domain - vertical upper left (RossSeaBox = 424/848, Full = 194/388)  260 280	
		mdh2 					=	517	  , & !	model domain - horizontal lower right (RossSeaBox = 374/748, Full = 517/1035)	
		mdv2 					=	525	   	  !	model domain - vertical lower right (RossSeaBox = 426/752, Full = 525/1051)	

	! ------------------------------------------------------------
	! start and end timing
	! ------------------------------------------------------------
	real(kind=dbl_kind), parameter :: &
		begin_j_day 	=	71	    ,        & !	start date (decimal julian day)	
		begin_year 		=	1997	    ,        & !	start year	
		end_j_day 		=	182	    ,        & !	end date (decimal julian day)	
		end_year 			=	1998	             !	end year	

	! ------------------------------------------------------------
	! output options
	! ------------------------------------------------------------
	integer(kind=int_kind), parameter :: &
		write_f 			=	24	    ,        & !	write frequency (hours)	
		write_disk 		=	1	    ,        & !	actually write to disk - may not want to if debugging (1=yes, 0=no)	
		wr_stations 	=	0	             !	write out station netcdf files - costly disk access, maybe want to supress	

	! ------------------------------------------------------------
	! grid resolution
	! ------------------------------------------------------------

	integer(kind=int_kind), parameter :: &
		n_max_floes 					= 40						, & ! maximum number of ice categories
		z_max_ice 						= 42						, & ! maximum ice layers
		z_max_snow						= 26						, & ! maximum snow layers
		z_max_pack 						= z_max_ice + z_max_snow, &
		z_int_min 						=	1	    				, & !	defines minimum number of model layers	
		z_sk 									=	2	    				, & !	number of layers in skeletal layer	
		pl_max 								= 30								! maximum platelet layers
	real(kind=dbl_kind), parameter :: &
		h_max 								= 10.0_dbl_kind , & ! 1maximum ice thickness (m)
		pl_th 								= 0.02_dbl_kind , & ! platelet layer thickness (m)
		bot_th 								= 0.2_dbl_kind  ,	&	! thickness of bottom section considered bottom ice (m from bottom)
		sk_th_min 						=	0.001_dbl_kind, & !	minimum height of individual skeletal layer	
		sk_th_max 						=	0.01_dbl_kind	, & !	maximum height of individual skeletal layer	
		z_th_min 							=	0.02_dbl_kind	 ,	& !	minimum height of individual layer	
		z_th_max 							=	0.3937376_dbl_kind,	& !	maximum height of individual layer	
		z_th_fr 							=	0.01_dbl_kind	,	& !	default new consolidated frazil ice thickness
		z_th_crit 						=	0.02_dbl_kind	,	& !	critical height of individual layer, below which a smaller physics timestep mut be used (m)
		h_crit                = z_th_min*(z_max_ice-z_sk), & ! height where grid shifts to accordion style
    th_min                = z_th_min*z_int_min

	! ------------------------------------------------------------
	! adjustable parameters
	! ------------------------------------------------------------

		! timing
		real(kind=dbl_kind), save :: &
			dt 						=	1	    ,        & !	length of time step (hours)	
			dt_sub_1 			=	0.2	    ,        & !	fraction of main timestep for ice physics - normal is 0.05	
			dt_sub_2 			=	0.05	    ,        & !	fraction of main timestep for fast-changing ice physics - normal is 0.005	
			dt_sub_3 			=	0.05	    ,        & !	fraction of main timestep for very-fast-changing ice physics	
			dt_sub_1_tol 	=	0.5	    ,        & !	minimum temperature tolerance for use of dt_sub_1 (deg C)	
			dt_sub_2_tol 	=	1.3	             !	minimum temperature tolerance for use of dt_sub_2 (dec C)	

		! ocean
		real(kind=dbl_kind), save :: &
			Fw 						=	7.0	    ,        & !	Oceanic heat flux from water (W/m^2) (1.195 cal/m^2/sec * 4.1868 W*s/cal = 5 W/m^2 in Arrigo, others?)	
			Sw 						=	34.1	    ,        & !	salinity of seawater (psu) - normal = 34.95 (33.33)	
			Sd 						=	1027693	    ,        & !	density of seawater  (g/m^3)	
			Tw 						=	-1.8	    ,        & !	seawater temperature (degC)	
			sw_NH4 				=	0.	    ,        & !	seawater NH4 concentration, used if climatologies not used (µMolar)	
			sw_NO3 				=	31	    ,        & !	seawater N03 concentration, used if climatologies not used (µMolar)	
			sw_SiO4 			=	80.	    ,        & !	seawater SiO2 concentration, used if climatologies not used (µMolar)	
			sw_PO4 				=	2.1	    ,        & !	seawater PO4 concentration, used if climatologies not used (µMolar)	
			sw_poc 				=	0.	    ,        & !	seawater detritus concentration (g/m^3)	
			alg_wc 				=	35	             !	water column algal concentration - used when freezing/flooding/creating new ice (mgC/m^3)	

		! snow & ice
		integer(kind=int_kind), parameter :: &
			ic_n 									= 5							, & ! ice thickness categories (between 1 and 10)
			sc_n 									= 2 						, & ! lognormal snow thickness categories (currently, either 1 or 2) 
			la 										= 1 						, & ! lognormal light adjustment categories (between 1 and 9)
			snow_model 		=	1	    ,	&	!	snow model type	
			icecon_mod 		=	2	    ,	&	!	modify satellite ice concentration (0=no modification, 1=90%->100%, 2=80%->100%)	
			grid_model 		=	2	    		!	determines how vertical ice grid is constructed (0 = equally spaced grid, 1 = higher resolution at top/bottom, 2 = variable number of layers)	
		real(kind=dbl_kind), parameter :: &
			ksnow 				=	0.33_dbl_kind,        & !	standard snow conductivity (W/m/K) (somewhere i found = 0.0741Cal/m/s/degC)	
			den_s_dry 		=	0.35_dbl_kind,        & !	density dry snow (g/cm^3)	
			den_s_wet 		=	0.35_dbl_kind,        & !	density wet snow (g/cm^3)	
			den_s_switch 	=	-2.0_dbl_kind,        & !	temperature above which wet snow falls, below which dry snow falls	
			bb_f 					=	0.03_dbl_kind,        & !	bubble fraction in sea ice - causes density to be decreased by this amount	
			eps0_s 				= 0.97_dbl_kind,     & ! long-wave emissivity of snow surface - 
			eps0_i 				= 0.99_dbl_kind,     & ! long-wave emissivity of ice surface - 
			eps_snow 			=	eps0_s	    , 		& !	long wave emissivity of snow (%)	
			eps_ice 			=	eps0_i	    , 		& !	long wave emissivity of ice (%)	
			atan_max 			= 1.557407724654902_dbl_kind,	&	! = tan(1), used for atan function multiplier			
      atan_c_i = atan_max/0.5_dbl_kind,  &	! 0.5 m is the cutoff for the atan function that describes albedo
			chw 					=	0.006_dbl_kind, 	& !	ocean-ice heat transfer coefficient
      mu_w_min 			= 0.05_dbl_kind,   	&	! minimum friction velocity between ice and ocean, for use in ocean heat flux calc
      Fm_a_switch 	= 0.90_dbl_kind,   	&	! grid cell area above which mixed layer frazil is applied to bottom of ice
      mean_floe_diameter 	= 30.0_dbl_kind,   	&	
      alpha_lateral_melt 	= 0.66_dbl_kind,   		& 
      melt_f_denom = mean_floe_diameter*alpha_lateral_melt 

		! ice desalination and fluid transfer
		integer(kind=int_kind), parameter :: &
			desal 				=	7	             !	desalination type (0 = normal, 1 = with brine replacement, 2 = dtt brine flux)	
		real(kind=dbl_kind), save :: &
			bv_conv 			=	200.	   ,        & !	critical brine volume at which desal switches from direct convection to relative conv.	
			f_sk 					=	0.5	    ,        & !	fraction of skeletal layer open to convection	
			fb 						=	0.0511	    ,        & !	fraction of brine tube layer that is brine tubes	
			fi 						=	0.5	    ,        & !	minimum fractional sea ice coverage neccessary for model use	
			vb_crit 			=	50	    ,        & !	brine volume above which gravity drainage occurs (ppt)	
			conv_max 			=	1.74E-05	    ,        & !	maximum convective flux for use in desalination (cm^3/cm^2/s)	
			dbvdt_scale 	=	1.0	             !	scaling factor for desal dilution			

		! irradiance
		real(kind=dbl_kind), parameter :: &
			alb_s_dry 		=	0.98	    ,        & !	standard albedo dry snow (%)	
			alb_s_wet 		=	0.88	    ,        & !	standard albedo wet snow (%)	
			alb_i_dry 		=	0.58	    ,        & !	standard albedo of sea ice (%)	
			alb_i_wet 		=	0.505	    ,        & !	standard albedo of sea ice (%)	
			h_snowpatch 	=	0.02	    ,        & !	contant that determines the bare ice/snow covered ice ratio for a given snowdepth	
			par_to_swd 		=	2	    ,        & !	energy conversion factor between PAR (400-700nm) irradiance to total shortwave downward irradiance (250-4000nm) (energy/energy units) - 2.034 calc clear sky	
			a_ice_ir 			=	7.18	             !	weighted-mean ice absorption coefficient for near-IR (700-4000nm)	
		integer(kind=int_kind), parameter :: &
			par_cf 				=	0	             !	use climatology could fraction to correct par (0=off (cf already included in par), 1=yes)	

		! biology
		integer(kind=int_kind), parameter :: &
			alg_mig 			=	1	             !	algae migration (0=algae stay put verticaly, 1=algae move with their respective layers while ice grows (quasi-movement))
		real(kind=dbl_kind), parameter :: &
			Ek_max 				=	18	    ,        & !	Spectral photoaclimation parameter (microEin*m-2*s-1)	
			A 						=	1.4	    ,        & !	parameter of light utilization equation (dimensionless)	
			B 						=	0.12	    ,        & !	parameter of light utilization equation (dimensionless)	
			rg 						=	0.0631	    ,        & !	growth rate constant rg (1/degC) - from Epply et al 1972	
			G0 						=	0.81	    ,        & !	growth rate @ zero dec C (1/day) - from Epply at al 1972	
			xp 						=	0.01	    ,        & !	phytoplankton death/grazing rate (1/day)	
			remin 				=	0.03	    ,        & !	rate of poc/detritus remineralization (1/day)	
			remin_f 			=	1	    ,        & !	fraction poc remineralization to available N,P	
			c_chl 				=	35	    ,        & !	Carbon:Chlorophyl ratio (grams/gram)	
			c_n 					=	7	    ,        & !	Carbon:Nitrogen Ratio (moles/mole)	
			c_si 					=	4	    ,        & !	Carbon:Silicon Ratio (moles/mole)	
			c_p 					=	106	    ,        & !	Carbon:Phosphorus Ratio (moles/mole)	
			Ks_NH4 				=	1	    ,        & !	Half-saturation rate constant (microMolar)	
			Ks_NO3 				=	1	    ,        & !	Half-saturation rate constant (microMolar)	
			Ks_SiO4 			=	60	    ,        & !	Half-saturation rate constant (microMolar)	
			Ks_PO4 				=	0.1	    ,        & !	Half-saturation rate constant (microMolar)	
			alg_dtt 			=	1	    ,        & !	algal model calculation frequency (0=once per time step, 1=same frequency as sub-dt ice physics)	
			alg_mig_crit 	=	1.5	    ,        & !	maximum growth rate under which algae maintain position (cm/day)	
			min_alg 			=	3.5	             !	minimum microalgal concentration (mgC/m^3)	

		! initialization
		integer(kind=int_kind), parameter :: &
			iit 					=	1	    ,        & !	intial ice temp (0 = use linear gradient w/ airtemp as inital ice temp, 1 = use water temp)	
			iis 					=	0	    ,        & !	initial ice salinity (0 = use s_cont constant salinity 1 = use standard 1st year ice "C" curve, 2 = multiyear ice curve)	
			iin 					=	4	    ,        & !	initial ice nutrients (0 = full nutrients upon ice creation, 1 = linear depleted nutrients, 2=very little nutrients, 4=seawater nutrients in brine)	
			iif 					=	0	             !	initial snow flooding (m) - used only in validation mode	
		real(kind=dbl_kind), save :: &
			s_const 			=	9	    ,        & !	salinity constant to use if iis paramater is set to 0	
			n_f 					=	0.3	             !	nutrient fraction - used to scale inital forcing nutirent concentrations if iin set to 2	

		! other parameters
		integer(kind=int_kind), save :: &
			ts_is_at 			=	0	    ,        & !	short-circuit atmospheric heat transfer by fixing surface temp to air temp (0=off 1=on)	
			kevin_light 	=	0	    ,        & !	use Kevin's light model instead of gregg&carder 1990 (0=off 1=on)	
			use_pl 				=	0	    ,        & !	use a platelet layer - only works if validaiton is on (0=off 1=on)	
			ncep_f 				=	6	    ,        & !	ncep/doe II forcing frequency (6 = every 6 hours, 24 = daily)	
			pr_on 				=	0	    ,        & !	on/off toggle for causing precipitation to go directly into snow ice formation, if appropriate	
			use_mdiff 		=	1	    ,        & !	use molecular diffusion to supply nutrients to ice bottom (0=off 1=on)	
			no_flooding 	=	0	    ,        & !	prevents flooding if non-zero	
			flood_brine 	=	0	    ,        & !	(1) pushed brine upward or (2) floods directly with seawater	
			woa_depth 		=	3	    ,        & !	level of World Ocean Altas data to use as forcing (1=surface,2=10m,3=20m,4=30m,5=50m,6=75m,7=100m)	
			max_it 				=	100	    ,        & !	Newton-Raphson maximum iterations	
			snow_ridging 	=	0	    ,	&	!	1=conserve snow mass during ridging,0=don't ridge up snow (w/ mass loss)	
			snow_in_gaps 	=	1	             !	determines whether snow is inserted into ridging gaps, or just seawater	
		real(kind=dbl_kind), save :: &
			min_sal 			=	0.1	    ,        & !	minimum bulk ice salinity (ppt)	
			nr_tol 				=	0.005	    ,        & !	Newton-Raphson method tolerance for surface temp (T0)	
			gl_max_f 			=	0.075	    ,        & !	fraction of minimum layer height that is the trigger for boundary adjustment	
			fl_max_f 			=	0.075	    ,        & !	fraction of minimum layer height that is the trigger for boundary adjustment for surface flooding	
			T_ref 				=	-75	    ,        & !	temp from which all heat intergrals are calculated (degC) (hopefully the lower than the lowest ice temp)	
			temp_tol 			=	0.0000001	    ,        & !	temp difference tolerance for heat methods - if the difference in temp get smaller than this, take corrective action	
			fl_crit 			=	-0.02	    ,        & !	minimum freeboard required to flood ice surface (m)	
			fl_ratio 			=	0.5	    ,        & !	ratio of ice/water in flooded snow	
			da_f 					=	0.4	    ,        & !	distribution adjustment factor - effectively changes the s.d. of the log-normal distribution - value between 1 and 0	
			snow_min 			=	0.01	    ,        & !	minimum modeled snow depth (m) - below this value, the model assumes no snow for thermodynamic/light/albedo putposes	
			gc_offset 		=	0	    ,        & !	offset in hours for indexing gc_par_ice file (240*24=5760 for arrigo 1982 sim)	
			snow_fudge 		=	1.0	    ,        & !	snow depth fudge factor for multiplicatively adjusting the snow depth for light purposes	
			a_factor 			=	0.3	    ,        & !	wavelength independent snow absorption for delta eddington light routine	
			p_factor 			=	1	    ,        & !	ice strength scalar, used to compare ice strength to ice momentum during advection refinement	
			snow_rd_lim 	=	24	    ,        & !	number that determines when snow depth distribution is re-distributed	
			cv_void_f 		=	0.3	    ,        & !	convergence void fraction (fraction)	
			cvf_switch 		=	0.6	    ,        & !	ice height below which cv_void_f is affected by cvf_thin (m)
			cvf_thin			=	0.33333333	    ,        & !		multiplier of cv_void_f when ice is less than cvf_switch (fraction)
			ohf_skew 			=	0	    ,        & !	ocean heat flux skew.  This number of added to the ocean heat flux (but it is always kept above 0)	
			at_skew 			=	0	    ,        & !	air temperature flux skew.  This number of added to the air temp	
			snow_skew 		=	0	             !	snow depth skew - fractional snow depth/precipitation multiplier (0.5 = 50% increase)	

	! ------------------------------------------------------------
	! tracers
	! ------------------------------------------------------------

		! external/2D tracers
		integer(kind=int_kind), parameter :: &
			n_2t									= 4  								! number of 2D tracers
		integer(kind=int_kind), parameter :: &
			age_i									= 1, 							& !	
			ridged_i							= 2, 							& !
			snow_dist_i						= 3, 							& !
			snow_rd_i							= 4                 !

		! internal/3D tracers
		integer(kind=int_kind), parameter :: &
			n_3t 									= 5  								! number of 3D tracers		
		integer(kind=int_kind), parameter :: &
			no3_i									= 1, 							& !
			nh4_i									= 2, 							& !
			po4_i									= 3, 							& !
			sio3_i								= 4, 							& !
			det_i									= 5, 							& !
			smalg_i								= 6									! 


	! ------------------------------------------------------------
	! ice fluxes - sign is out of ice
	! ------------------------------------------------------------
		integer(kind=int_kind), parameter :: &
			n_flx									= 26  						! number of ice fluxes
		integer(kind=int_kind), parameter :: &
			w_smelt 						= 1, 							& !
			w_bmelt 						= 2, 							& !
			w_bfreeze 					= 3, 							& !
			w_flood 						= 4, 							& !
			w_frfreeze 					= 5, 							& !
			w_latmelt 					= 6, 							& !
			w_latfreeze 				= 7, 							& !
			w_snowmelt 					= 8, 							& !
			w_desal 						= 9, 							& !
			s_smelt 						= 10, 							& !
			s_bmelt 						= 11, 							& !
			s_bfreeze 					= 12, 							& !
			s_flood 						= 13, 							& !
			s_frfreeze 					= 14, 							& !
			s_latmelt 					= 15, 							& !
			s_latfreeze 				= 16, 							& !
			s_desal 						= 17, 							& !
			q_ssmelt 						= 18, 							& !
    	q_mlmelt 						= 19, 							& !
    	q_bcon 							= 20, 							& !
    	q_latmelt 					= 21,	 							& !
    	q_totalfreeze 			= 22,	 							& !
    	q_totalmelt 				= 23,	 							& !
    	v_frazil      			= 24,	 							& !
    	v_congelation 			= 25,	 							& !
    	v_snowice     			= 26	 							  !


	! ------------------------------------------------------------
	! ice thickness change tracking variables
	! ------------------------------------------------------------
		integer(kind=int_kind), parameter :: &
			n_dh									= 17  						! number of ice thickness changes
		integer(kind=int_kind), parameter :: &
			sn_precip 					= 1, 							& !
			sn_melt 						= 2, 							& !
			sn_subl 						= 3, 							& !
			sn_flood						= 4,							&	!
			ice_b_melt_ml 			= 5, 							& !
			ice_b_melt_con 			= 6, 							& !
			ice_b_grow_ml 			= 7, 							& !
			ice_b_grow_con 			= 8, 							& !
			ice_s_melt 					= 9, 							& !
			ice_s_subl 					= 10, 							& !
			ice_s_flood 				= 11, 							& !
			ice_s_drain 				= 12, 							& !
			ice_s_melt_flood 		= 13, 							& !
			ice_b_dhdt 					= 14, 							&	!
			ice_b_sal 					= 15,								&	!
			ice_s_sal 					= 16,								&	!
			sn_rain   					= 17 									!
		
	! ------------------------------------------------------------
	! input data grid indexes
	! ------------------------------------------------------------          
	integer(kind=int_kind), parameter :: &
			mp_x									= 192						, & ! NCEP/DOE II x bound
			mp_y									= 94						, & ! NCEP/DOE II y bound
			ec_x									= 144						, & ! ECMWF 40-yr Reanalysis x bound
			ec_y									= 73						, & ! ECMWF 40-yr Reanalysis y bound
			eci_x									= 240						, & ! ECMWF Interim Reanalysis x bound
			eci_y									= 121						, & ! ECMWF Interim Reanalysis y bound
			dg_x									= 360						, & ! 1 degree grid x bound
			dg_y									= 180						, & ! 1 degree grid y bound
			wavl									= 31						, & ! number of shortwave irradiance bins
			logn_n 								= 9							    ! lognormal distribution bins - not sure if this is used?

		! ------------------------------------------------------------
		! model domain grid indexes based on southern hemisphere EASE grid
		! ------------------------------------------------------------          
		integer(kind=int_kind), parameter :: &
			grid_v								= 721						, & !
			grid_h								= 721						, & !
			h1_ncep_i							= 202						, & !
			h2_ncep_i							= 518						, & !
			v1_ncep_i							= 194						, & !
			v2_ncep_i							= 526						  !

		! ------------------------------------------------------------
		! Enumerations
		! ------------------------------------------------------------
		integer(kind=int_kind), parameter :: &
			ncep_at								= 1					, & !
			ncep_p 								= 2					, & !
			ncep_h 								= 3					, & !
			ncep_fc								= 4					, & !
			ncep_u10							= 5					, & !
			ncep_v10							= 6					, & !
			ncep_pr								= 7							!
		integer(kind=int_kind), parameter :: &
			woa_t									= 1					, & !
			woa_s									= 2					, & !
			woa_n									= 3					, & !
			woa_p									= 4					, & !
			woa_si								= 5						  !

		! ------------------------------------------------------------
		! All kinds of parameters
		! ------------------------------------------------------------
    real(kind=dbl_kind), parameter :: &
    	ki_min 								= 0.563_dbl_kind			, & !  minimum sea ice thermal conductivity (W/m/K)
			pond_f_perc 					= 0.3e0_dbl_kind	    , & ! fraction of ponded/melt water that pushed down through brine network
			af_min 								= 5.0e-9_dbl_kind     , & ! minimum areal fraction of ice in a grid cell
      Ce 										= 2.1e-3_dbl_kind      , & ! coefficient of turbulent latent heat transfer (water vapor?)
      Ch 										= 2.0e-3_dbl_kind      , & ! coefficient of turbulent sensible heat transfer
      cp 										= 1005.0_dbl_kind       ! J/kg/K 

		! ------------------------------------------------------------
		! Model globals/parameters assigned during run or from
		! constant.txt run file
		! ------------------------------------------------------------
	
		real(kind=dbl_kind), allocatable :: &
				ida_multiplier(:),        &	!
				lda_multiplier(:),        &	!
				sda_multiplier(:)	        	!
		logical(kind=log_kind), save :: &
				leap_year,       					&	!
				start,      							&	!
				do_write,       					&	!
				z_odd												!
		real(kind=dbl_kind), save, dimension (301) :: &
				aice, &
				aph, &
				awater, &
				awater_tc, &
				awater_sc, &
				awater_v, &
				aozone, &
				rdsnow, &
				rwsnow, &
				surfacelight
		character (LEN=14) :: &
				datestr
		integer(2), save, pointer :: &
				SSMI_grid_int2(:,:), &
				icevec_grid_int2(:,:,:), &
				Ed_int2(:,:,:), &
				mp_grid_int2(:,:), &
				ec_grid_int2(:,:), &
				eci_grid_int2(:,:)
		real(kind=dbl_kind), save, pointer :: &
				kds_wet(:), &
				kds_dry(:), &
				lambda(:), &
				quanta_to_watts(:)              
		real(kind=dbl_kind), save, dimension (130) :: &
				mc_prod
		real(kind=dbl_kind), save, dimension (ic_n) :: &
				ic_h_max, &
				ic_h_med, &
				ic_h_min              

		integer(kind=int_kind), save :: &
			steps, last_day, last_hour, last_year, icevec_index,skdump, &
			dtt_1,dtt_2,dtt_3,n_dtt_1,n_dtt_2,n_dtt_3,snow_loaded, &
			mdh, mdv, last6hour,last12hour,n_stations,snowd_index, &
			hour24,tcells,sda_n, last3hour,ida_norm_n,lda_n,lda_norm_n, &
			f_index,f_index_next,i_temp,pur_clock,icecon_index

		real(kind=dbl_kind), save :: &
			pur_stepper,sk_h,bt_h,next_write_hour,aph_max,dt_s, &
			b_flux_max_tot,dt_years,dt_days, &
			dtt_s_1,dtt_s_2,dtt_s_3,cur_month,ida_d,lda_d,ad_denom, &
			total_flooded,a_ice_ir_cos
			
end module sia2_parameters
	