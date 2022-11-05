module sia2_types

  use sia2_constants
  use sia2_parameters
  
  implicit none

  public

  ! ------------------------------------------------------------
  ! Static array types (forcing & projection data)
  ! ------------------------------------------------------------

  ! Projection type - holds projection conversion information for various projections
  type, public :: proj_type
      integer(kind=int_kind), dimension(grid_h,grid_v) :: &
          mdh, &          ! model domain horizontal grid index
          mdv, &          ! model domain vertical grid index
          x_mp, &         ! NCEP/mercator grid horizontal(longitude) index
          y_mp, &         ! NCEP/mercator grid vertical(latitude) index
          x_ec, &         ! ECMWF 2.5 grid horizontal(longitude) index
          y_ec, &         ! ECMWF 2.5 grid vertical(latitude) index
          x_eci, &        ! ECMWF 1.5 grid horizontal(longitude) index
          y_eci, &        ! ECMWF 1.5 grid vertical(latitude) index
          x_dg, &         ! 1 degree grid horizontal(longitude) index
          y_dg, &         ! 1 degree grid vertical(latitude) index
          mi, &           ! grid index to modeled variables
          mask            ! ocean mask - 1 = ocean, 0 = excluded from analysis (land or coastline)
      real(kind=dbl_kind), dimension(grid_h,grid_v) :: &
          lat, &          ! central latitude of model grid cell
          lon             ! central longitude of model grid cell
  end type proj_type

  ! NCEP type - holds an entire set of NCEP forcing information
  type, public :: ncep_type
      ! air temp, 2m air pressure, specific humidity, cloud fraction, 
      ! u direction wind speed, v direciton wind speed, precipitation rate
      real(kind=dbl_kind), dimension(mp_x,mp_y) :: &
          at, &           ! 2m air temperature (kelvin)
          p, &            ! surface air pressure (mb)
          h, &            ! specific humidity ()
          fc, &           ! fraction cloud cover (%)
          u10, &          ! 10m horizontal wind speed vector (m/s) 
          v10, &          ! 10m vertical wind speed vector (m/s) 
          pr              ! precipitation rate (kg/m^2/s)
  end type ncep_type

  ! ed_type - holds a set of global spectral radiation based on the
  ! NCEP/DOE grid
  type, public :: ed_type
      real(kind=dbl_kind), dimension(mp_x,mp_y,wavl) :: &
          Ed              ! µEin/m^2/s
  end type ed_type

  ! mp_f_type - Mercator projection forcing type, holds all NCEP forcing data
  ! using the ~2 degree NCEP grid
  type, public :: mp_f_type  
      type (ncep_type) :: &
          ncep            ! current NCEP forcing dataset
      type (ncep_type) :: &
          ncep_next       ! next available NCEP forcing dataset
      type (ncep_type) :: &
          ncep_interp     ! NCEP forcing interpolated for current time step
      type (ed_type), dimension(2) :: &
          Edir, &         ! Direct irradience (µEin/m^2/s)
          Edif            ! Diffuse irradience (µEin/m^2/s)
  end type mp_f_type

  ! ECMWF type - holds an entire set of ECMWF forcing information
  type, public :: ecmwf_type
      ! air temp, 2m air pressure, specific humidity, cloud fraction, 
      ! u direction wind speed, v direciton wind speed, precipitation rate
      real(kind=dbl_kind), dimension(ec_x,ec_y) :: &
          at, &           ! 2m air temperature (kelvin)
          p, &            ! surface air pressure (Pa)
          dpt, &          ! dew point temperature (kelvin)
          fc, &           ! fraction cloud cover (%)
          u10, &          ! 10m horizontal wind speed vector (m/s) 
          v10, &          ! 10m vertical wind speed vector (m/s) 
          pr, &           ! total precipitation (m)
          ssr             ! surface solar radiation (W/m^2)
  end type ecmwf_type

  ! ECMWF Interim type - holds an entire set of ECMWF forcing information
  type, public :: ecmwf_int_type
      ! air temp, 2m air pressure, specific humidity, cloud fraction, 
      ! u direction wind speed, v direciton wind speed, precipitation rate
      real(kind=dbl_kind), dimension(eci_x,eci_y) :: &
          at, &           ! 2m air temperature (kelvin)
          p, &            ! surface air pressure (Pa)
          dpt, &          ! dew point temperature (kelvin)
          fc, &           ! fraction cloud cover (%)
          u10, &          ! 10m horizontal wind speed vector (m/s) 
          v10, &          ! 10m vertical wind speed vector (m/s) 
          pr, &           ! total precipitation (m)
          ssr             ! surface solar radiation (W/m^2)
  end type ecmwf_int_type

  ! ec_f_type - ECMWF projection forcing type, holds all ECMWF forcing data
  ! using the ECMWF grid
  type, public :: ec_f_type  
      type (ecmwf_type) :: &
          ecmwf           ! current NCEP forcing dataset
      type (ecmwf_type) :: &
          ecmwf_next      ! next available NCEP forcing dataset
      type (ecmwf_type) :: &
          ecmwf_interp    ! NCEP forcing interpolated for current time step
  end type ec_f_type

  ! eci_f_type - ECMWF Interim projection forcing type, holds all ECMWF Interim forcing data
  ! using the ECMWF grid
  type, public :: eci_f_type  
      type (ecmwf_int_type) :: &
          ecmwf           ! current NCEP forcing dataset
      type (ecmwf_int_type) :: &
          ecmwf_next      ! next available NCEP forcing dataset
      type (ecmwf_int_type) :: &
          ecmwf_interp    ! NCEP forcing interpolated for current time step
  end type eci_f_type

  ! woa_type - holds 1 degrees World Ocean Atlas data for a specific depth,
  ! which is determined by setting the woa_depth constant
  type, public :: woa_type
      real, dimension(dg_x,dg_y) :: &
          t, &            ! World Ocean Atlas water temp (dec C)
          s, &            ! World Ocean Atlas salinity (psu)
          n, &            ! World Ocean Atlas nitrate concentration (µMol)
          p, &            ! World Ocean Atlas phosphate concentration (µMol)
          si              ! World Ocean Atlas silica concentration (µMol)
  end type woa_type

  ! degree grid forcing type - holds forcing datasets available in a
  ! 1 degree grid
  type, public :: dg_f_type
      type (woa_type) :: &
          woa             ! current WOA forcing dataset
      type (woa_type) :: &
          woa_next        ! next available WOA forcing dataset
      type (woa_type) :: &
          woa_interp      ! WOA forcing interpolated for current time step
  end type dg_f_type

  ! EASE projection forcing type - holds forcing variables that are
  ! projected using ease grid
  type, public :: ease_f_type 
      real(kind=dbl_kind), dimension(grid_h,grid_v) :: &
          icecon, &       ! ssm/i ice concentration for current timestep (fraction)
          icecon_next     ! pre-loaded ssm/i ice concentration for next step (fraction)
      real(kind=dbl_kind), dimension(grid_h,grid_v) :: &
          sh, &           ! ssm/i snow depth for current timestep (cm)
          sh_next         ! pre-loaded ssm/i snow depth for next step (cm)
      real(kind=dbl_kind), dimension(grid_h,grid_v,7) :: &
          sh_past         ! array of past snow depths for validaiton of snow depth algorithm
      real(kind=dbl_kind), dimension(grid_h,grid_v,5) :: &
          icecon_past         ! array of past snow depths for validaiton of snow depth algorithm
      real(kind=dbl_kind), dimension(grid_h,grid_v,3) :: &
          ice_vec, &       ! ease grid ice vectors (u,v,flag)
          ice_vec_next     ! pre-loaded ease grid ice vectors for next step(u,v,flag)
      real(kind=dbl_kind), dimension(grid_h,grid_v,3,5) :: &
          icevec_past
  end type ease_f_type

  type, public :: validation_f_type ! validation station forcing type
      real(kind=dbl_kind) :: &
          at, &
          sh, &
          ih, &
          p, &
          lat, &
          lon, &
          iceh1, &
          ws, &
          rhum, &
          fc, &
          pr
      real(kind=dbl_kind) :: &
          fw, &
          t, &
          s, &
          no3, &
          nh4, &
          po4, &
          sioh4, &
          swd, &
          lwd
      integer(kind=int_kind) :: &
          time, &
          time_offset, &
          x_mp, &
          y_mp, &
          x_ec, &
          y_ec, &
          x_eci, &
          y_eci, &
          x_dg, &
          y_dg, &
          length, &
          begin_year, &
          mdh, &
          mdv
  end type validation_f_type

  ! ------------------------------------------------------------
  ! Dynamically allocated array types (modeled data)
  ! ------------------------------------------------------------

  ! Interpolated/corrected forcing data is stored for one modeled grid point
  ! is stored in this type
  type, public :: forcing_type 
      real(kind=dbl_kind) :: &
          at, &
          p, &
          h, &
          fc, &
          u10, &
          v10, &
          pr, &
          ws, &
          rhum, &
          swd, &
          lwd, &
          fw, &
          t, &
          s, &
          d, &
          no3, &
          nh4, &
          po4, &
          sioh4, &
          poc, &
          sh_interp, &
          sd_interp, &
          sh_new, &
          sh_interp_next, &
          sd_interp_next, &
          ic_interp, &
          ic_interp_next, &
          sh_interp_last, &
          sd_interp_last, &
          ivu_interp, &
          ivv_interp, &
          ivu_interp_next, &
          ivv_interp_next
      real(kind=dbl_kind), dimension(n_3t) :: &
          tr3
  end type forcing_type

  type, public :: meta_type
      integer(kind=int_kind) :: &
          status, &
          grid_h, &
          grid_v, &
          x_mp, &
          y_mp, &
          x_ec, &
          y_ec, &
          x_eci, &
          y_eci, &
          x_dg, &
          y_dg
      real(kind=dbl_kind) :: &
          melt_loss, &
          adv_loss, &
          md_loss, &
          cong_growth, &
          snow_growth, &
          adv_gain, &
          md_gain, &
          a_convg, &
          a_new, &
          a_drop, &
          lat, &
          lon, &
          pr_clim, &
          pr_ssmi, &
          pxl_h_offset, &
          pxl_v_offset, &
          prod_sum_int, &
          prod_sum_bot, &
          bm_lost, &
          tlim, &
          llim, &         ! light limitation fraction (dimensionless)
          nlim, &         ! nitrogen limitation fraction (dimensionless)
          plim, &         ! phosphorus limitation fraction (dimensionless)
          silim, &        ! silica limitation fraction (dimensionless)
          slim, &        ! salinity limitation fraction (dimensionless)
          salt_flux, &
          h2o_flux, &
          p_wgt_int, &
          p_wgt_af_int, &
          p_wgt_bot, &
          p_wgt_af_bot, &
          pwid_sum_int, &
          pwid_sum_af_int, &
          pwsd_sum_int, &
          pwsd_sum_af_int, &
          pwt_sum_int, &
          pwt_sum_af_int, &
          pwid_sum_bot, &
          pwid_sum_af_bot, &
          pwsd_sum_bot, &
          pwsd_sum_af_bot, &
          pwt_sum_bot, &
          pwt_sum_af_bot
          
      integer(kind=int_kind), dimension(8) :: &
          ia, &           ! i_ease index of 8 surround cells
          ja, &           ! j_ease index of 8 surrounding cells
          mia             ! mi index of 8 surroundind cells
  end type meta_type

  type, public :: snow_type
      integer(kind=int_kind) :: &
          z
      real(kind=dbl_kind) :: &
          depth, &
          d_small, &
          ts
      real(kind=dbl_kind), dimension(z_max_snow) :: &
          t, &
          d, &
          heat, &
          th, &
          melt, &
          ed_w
  end type snow_type

  type, public :: pond_type
      real(kind=dbl_kind) :: &
          t, &
          s, &
          d, &
          heat, &
          th, &
          perc, &
          smalg, &        ! brine-based algal concentration (mgC/m^3)
          poc, &          ! particulate organic carbon (detritus) (mgC/m^3)
          no3, &          ! ice layer brine NO3 concentration (µMol)
          nh4, &          ! ice layer brine NH4 concentration (µMol)
          po4, &          ! ice layer brine PO4 concentration (µMol)
          sioh4           ! ice layer brine SiOH4 concentration (µMol)                  
      real(kind=dbl_kind), dimension(n_3t) :: &
          tr3
  end type pond_type


  type, public :: ice_type
  integer(kind=int_kind) :: &
      z               ! number of active layers
  integer(kind=int_kind), dimension(z_max_ice) :: &
          bced, &         ! brine channeled - indicated whether layer has good connection with underlying seawater
          drained         ! layer is drained of brine, affecting optical properties and thermal conductivity and capactitance
      real(kind=dbl_kind) :: &
          ts, &           ! surface temperature of top layer (ice OR snow) (degC)
          fbh, &          ! freeboard height of ice pack (m)
          pr_sum, &       ! precipitation sum - record total precipitation in cell (kg/m^2)
          c_gl_bv, &      ! congelation growth layer brine volume - used to determine density new ice growth
          sh_offset, &    ! offset for calculating thue snow height following flooding (m) - used with validation data
          sh_prev, &      ! snow height from previous time step (m) - used to calculate incremental snow height change during timestep 
          af, &           ! area fraction - records fraction of grid cell that this ice type represents
          Ed0, &
          Ed1, &
          PAR_bot, &
          PAR_bot_pl, &
          ed0_nir_dir, &
          ed0_nir_dif, &
          age, &
          ridged, &
          snow_dist, &
          snow_rd
          
      real(kind=dbl_kind), dimension(z_max_ice) :: &
          t, &            ! ice layer temperature (degC)
          s, &            ! ice layer bulk salinity (psu)
          d, &            ! ice layer density (g/m^2)
          bs, &           ! ice layer brine salinity (psu)
          bd, &           ! ice layer brine density (g/m^2)
          bv, &           ! ice layer brine volume (ppt)
          llim, &         ! light limitation fraction (dimensionless)
          nlim, &         ! nitrogen limitation fraction (dimensionless)
          plim, &         ! phosphorus limitation fraction (dimensionless)
          silim, &        ! silica limitation fraction (dimensionless)
          slim, &         ! salinity limitation fraction (dimensionless)
          heat, &         ! ice layer weirdo-enthalpy (J/g)
          th, &           ! ice layer thickness (m)
          id, &           ! ice layer bottom depth (m)
          smalg, &        ! brine-based algal concentration (mgC/m^3)
          prod, &         ! time-step productivity (mgC/grid cell)
          Ik1, &          ! time-step Ik-prime (µEin/m^2/s)
          dsdt, &         ! desal rate (psu/timestep)
          fbv, &          ! brine volume flux m^3/m^2/timestep
          dhdt_conv, &       ! storage for timestep-derived dhdt from last timestep, since desal occurs in reponse to heat flux
          f0, &
          dsdt3, &
          tgrad, &
          ed_w
      real(kind=dbl_kind), dimension(z_max_ice+1) :: &
          poc, &          ! particulate organic carbon (detritus) (mgC/m^3)
          no3, &          ! ice layer brine NO3 concentration (µMol) - instantaneous concentration
          nh4, &          ! ice layer brine NH4 concentration (µMol) - instantaneous concentration
          po4, &          ! ice layer brine PO4 concentration (µMol) - instantaneous concentration
          sioh4, &        ! ice layer brine SiOH4 concentration (µMol) - instantaneous concentration
          no3_mean, &          ! ice layer brine NO3 concentration (µMol) - conc. available to algae, over time
          nh4_mean, &          ! ice layer brine NH4 concentration (µMol) - conc. available to algae, over time
          po4_mean, &          ! ice layer brine PO4 concentration (µMol) - conc. available to algae, over time
          sioh4_mean           ! ice layer brine SiOH4 concentration (µMol) - conc. available to algae, over time
      real(kind=dbl_kind), dimension(n_2t) :: &
          tr2
      real(kind=dbl_kind), dimension(n_3t,z_max_ice) :: &
          tr3
      real(kind=dbl_kind), dimension(n_flx) :: &
          flux
      type (snow_type) :: &
          snow         ! snow_type - an ice grid cell may have more than one snow_type
      type (pond_type) :: &
          pond            ! pond_type - melt pond variables
      real(kind=dbl_kind), dimension(z_max_ice,la) :: &
          pur             ! Photosynthetically useable radiation by depth and wavelength (µEin/m^s/s/wavelength)
  end type ice_type


  type, public :: ice_pack_type
  
    real(kind=dbl_kind) :: &
      t_mean, &             ! mean surface temperature
      s_mean, &             ! mean bulk salinity (psu)
      d_mean, &             ! mean ice density
      th_mean, &            ! mean ice thickness
      q_mean, &             ! mean ice enthalpy (J/g)
      ts_mean, &            ! mean ice or snow surface temperature
      t1_mean, &            ! mean snow-ice interface temperature
      smalg_mean, &         ! algal concentration
      smalg_sum, &          ! total algal biomass  
      af_total, &           ! area fraction - records fraction of grid cell that this ice type represents
      ff_mean, &            ! area fraction of ice with surface flooding
      t_snow_mean, &        ! mean snow temperature
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
    real(kind=dbl_kind), dimension(n_2t) :: &
      tr2_mean
    real(kind=dbl_kind), dimension(n_3t) :: &
      tr3_mean
    type(ice_type), dimension(n_max_floes) :: &
      ice

  end type ice_pack_type


  type, public :: platelet_type
      real(kind=dbl_kind) :: &
          t, & 
          s, & 
          d, & 
          bs, & 
          bd, & 
          bv 
      real(kind=dbl_kind), dimension(pl_max) :: &
          no3, &
          nh4, &
          po4, &
          sioh4, &
          Ik1, &
          smalg, &
          poc, &
          pur
  end type platelet_type                  

  type, public :: snow_adv_type
      integer(kind=int_kind) :: &
          z
      real(kind=dbl_kind) :: &
          sh_prev, &
          depth           ! new advected snow depth, before ridging
      real(kind=dbl_kind), dimension(z_max_snow) :: &
          th_new
  end type snow_adv_type

  type, public :: ice_adv_type
      integer(kind=int_kind) :: &
          z               ! new number of ice layers resulting from advection
      real(kind=dbl_kind) :: &
          af, &           ! areal fraction of advected ice in this category, not yet ridged, could be > 1!
          depth, &           ! new advected internal ice depth depth (before ridging)
          t_mean
      real(kind=dbl_kind), dimension(z_max_ice) :: &
          th_new, &       ! new ice layer thicknesses resulting from advection (m)
          id_new, &       ! new ice layer depths resulting from advection (m)
          th_corr, &      ! corrected ice layer thicknesses resulting from boundaries during advection (m)
          th_debug, &     ! used to record particular ice thicknesses during advection for debugging purposes (m)
          th_snow_new, &  ! new snow layer thicknesses resulting from advection (m)
          th_old          ! previous ice layer thicknesses (before advection) (m)
  end type ice_adv_type

  type, public :: adv_type              
      integer(kind=int_kind) :: &
          v_mod           
      real(kind=dbl_kind) :: &
          adv_sum, &      ! sum of in and out advection fractions
          out, &          ! areal percentage of cell advected out
          out1, &         ! convergence-optimzed areal percentage of cell advected out
          icecon_corr, &  ! fractional correction of icecon due to next ice concentration forcing (dimensionless)
          a, &            ! ice-covered area of cell (m^2)
          af, &            ! ice-covered fraction of cell (m^2)
          a_new, &        ! ice-covered area of cell after advection(m^2)
          a_out, &        ! ice-covered area advected out of cell(m^2)
          a_in, &         ! ice-covered area advected into cell(m^2)
          p, &            ! ice strength (N)
          cvm0, &         ! convergence minimization parameter, original
          cvm1, &         ! convergence minimization parameter, currently optimized
          v_scale, &      ! velocity scalar, used to change velocity using shared memory
          mass, &         ! mean ice pack mass (kg/m^3)
          a_open, &       ! pre-advection area open water
          af_open, &      ! pre-advection fraction of open water
          af_new, &       ! area fraction of new ice, used to ad new ice during redistribution of ice categories
          a_convg, &      ! area of input converging
          f_convg, &      ! fraction of input converging
          a_melt, &
          a_allowed, &
          id_mean, &
          a_drop
      real(kind=dbl_kind), dimension(8)  :: &
          in, &           ! convergence-optimzed areal percentage inputs
          in1, &          ! areal percentage of bordering cell advected in - before advection
          in1a            ! area of advection in                         
!       real(kind=dbl_kind), dimension(0:8)  :: &
!         a_convg, &      ! area of input converging
!         f_convg         ! fraction of input converging
      real(kind=dbl_kind), dimension(ic_n)  :: &
          a_rdg           ! final total area that goes ends up ridging according to rf (km^2) 
      real(kind=dbl_kind), dimension(ic_n,ic_n)  :: &
          rf              ! ridging factor - areal reduction scalar for ridged ice
      type (ice_adv_type), dimension(sc_n,ic_n) :: &
          ice            ! snow_adv_type - holds advection variables for snow layers
      type (snow_adv_type), dimension(sc_n,ic_n) :: &
          snow            ! snow_adv_type - holds advection variables for snow layers
  end type adv_type

  type, public :: stn_write_type
      logical :: &
          valid           ! turns on/off writing for station
      integer(kind=int_kind) :: &
          step            ! index to the current writeout step
      character (LEN=80) :: &
          fname           ! string holds output filename for station
      real(kind=dbl_kind) :: &
          time            ! write out model time
      real(kind=dbl_kind), dimension(18) :: &
          s               ! array to hold single-value station variables
      real(kind=dbl_kind), dimension(29,z_max_ice) :: &
          z               ! array to hold z-dimension station variables
      real(kind=dbl_kind), dimension(5,z_max_ice+1) :: &
          z1              ! array to hold z1-dimension station variables
      real(kind=dbl_kind), dimension(3,wavl) :: &
          wavl            ! array to hold wavl dimension station variables
      real(kind=dbl_kind), dimension(3) :: &
          icevec          ! array to station icevec data
  end type stn_write_type

  
  contains
  
end module sia2_types
