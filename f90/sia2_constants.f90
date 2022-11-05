module sia2_constants

	implicit none

	public

	! ------------------------------------------------------------
	! Elemental data types/sizes
	! ------------------------------------------------------------

	integer, parameter :: &
			int_kind							= 4						, & ! integers are 4 bytes
      log_kind 							= kind(.true.), & ! logical kind
			real_kind							= 4						, & ! floats are 4 bytes in KPP
			dbl_kind							= 4								! doubles are set to 4 for KPP

	! ------------------------------------------------------------
	! Empirical and Physical Constants
	! ------------------------------------------------------------

	real(kind=dbl_kind), parameter :: &
			cw										=	3.96_dbl_kind		, &	! heat capacity of water at 0degC (J g-1 K-1)
			ci0										=	2.113_dbl_kind		, &	! constant in ice heat capacity equation (W m-1 K-1)
			IceD									=	918000_dbl_kind			, &	! density of pure ice (g m-3)
			mu										=	-0.054_dbl_kind   ,	&	! linear liquidus constant for seawater
			kelvin0								= 273.15_dbl_kind	 	  , &  !	0degC in K
      qq1 									= 1.16378e7_dbl_kind  , & ! q1 constant in turbulent latentent heat flux calc, CICE v. 4
      qq2 									= 5897.8_dbl_kind     , &  ! q2 constant in turbulent latentent heat flux calc, CICE v. 4
      Lv 										= 2501.0_dbl_kind , & ! kJ/kg latent heat of vaporization of water (2260)
      Lf 										= 334.0_dbl_kind , & ! kJ/kg latent heat of fusion of water (334)
      R_air 								= 287.058_dbl_kind  , & ! Universal gas constant for dry air - J/kg/K; pressure must be in Pascals for: density = pressure/RT   
      R_h2o 								= 461.495_dbl_kind  , &! Universal gas constant for water vapor - J/kg/K; pressure must be in Pascals for: density = pressure/RT
      heat_snow0 						= 0.2309_dbl_kind*kelvin0 + 0.0034_dbl_kind*kelvin0**2

		real(kind=dbl_kind), parameter :: pi=3.141592_dbl_kind
		real(kind=dbl_kind), parameter :: steph_boltz=5.66e-8_dbl_kind !W/m^2/K^4
		real(kind=dbl_kind), parameter :: gC_mC = 12.01_dbl_kind !gramsC/molesC
		real(kind=dbl_kind), parameter :: mC_gC = 1.0_dbl_kind/12.01_dbl_kind !molesC/gramC
		real(kind=dbl_kind), parameter :: cell_side = 25.067525_dbl_kind ! km
		real(kind=dbl_kind), parameter :: cell_area = cell_side**2 ! km^2
		real(kind=dbl_kind), parameter :: grav = 9.80616_dbl_kind ! gravity, m/s^2
	
		integer(kind=int_kind), parameter :: &
			bins 					=	31	      !	num of wavelength bins	

		real(kind=dbl_kind), parameter :: fe_A = 2.7798e-6_dbl_kind   ! hPa/K^4, turbutent heat flux coefficient
		real(kind=dbl_kind), parameter :: fe_B = -2.6913e-3_dbl_kind  ! hPa/K^3, turbutent heat flux coefficient
		real(kind=dbl_kind), parameter :: fe_C = 0.97920_dbl_kind     ! hPa/K^2, turbutent heat flux coefficient
		real(kind=dbl_kind), parameter :: fe_D = -158.64_dbl_kind     ! hPa/K, turbutent heat flux coefficient
		real(kind=dbl_kind), parameter :: fe_E = 9653.2_dbl_kind      ! hPa, turbutent heat flux coefficient

		real(kind=dbl_kind), parameter :: ss0 = 8.0_dbl_kind     ! 1st year ice salinilty standard curve value (top)
		real(kind=dbl_kind), parameter :: ss1 = 6.3_dbl_kind     ! 1st year ice salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ss2 = 5.6_dbl_kind     ! 1st year ice salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ss3 = 5.3_dbl_kind    ! 1st year ice salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ss4 = 5.2_dbl_kind     ! 1st year ice salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ss5 = 5.1_dbl_kind     ! 1st year ice salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ss6 = 4.9_dbl_kind     ! 1st year ice salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ss7 = 4.8_dbl_kind     ! 1st year ice salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ss8 = 4.8_dbl_kind     ! 1st year ice salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ss9 = 6.2_dbl_kind     ! 1st year ice salinilty standard curve value (bottom)

		real(kind=dbl_kind), parameter :: ssm0 = 0.1_dbl_kind     ! multi-year salinilty standard curve value (top)
		real(kind=dbl_kind), parameter :: ssm1 = 0.2_dbl_kind     ! multi-year salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ssm2 = 0.2_dbl_kind     ! multi-year salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ssm3 = 0.6_dbl_kind     ! multi-year salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ssm4 = 1.9_dbl_kind     ! multi-year salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ssm5 = 3.1_dbl_kind     ! multi-year salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ssm6 = 3.3_dbl_kind     ! multi-year salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ssm7 = 3.4_dbl_kind     ! multi-year salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ssm8 = 3.9_dbl_kind     ! multi-year salinilty standard curve value (mid)
		real(kind=dbl_kind), parameter :: ssm9 = 6.2_dbl_kind     ! multi-year salinilty standard curve value (bottom)

		real(kind=dbl_kind), parameter :: a0 = -0.12073_dbl_kind    ! salinity limiting poly coefficient 0 (dimensionless)
		real(kind=dbl_kind), parameter :: a1 = 0.07097_dbl_kind     ! salinity limiting poly coefficient 1 (dimensionless)
		real(kind=dbl_kind), parameter :: a2 = -0.00133_dbl_kind    ! salinity limiting poly coefficient 2 (dimensionless)
		real(kind=dbl_kind), parameter :: a3 = 6.3427e-6_dbl_kind  ! salinity limiting poly coefficient 3 (dimensionless)

		real(kind=dbl_kind), parameter :: a_star = 0.05_dbl_kind  ! mean value of areal fraction participating in ridging
		real(kind=dbl_kind), parameter :: r_a_star = 20.0_dbl_kind  ! 1/a_star
		real(kind=dbl_kind), parameter :: adv_mu = 4.0_dbl_kind  !  used to calc lambda constant in ice category redistribution function

		! seawater equation of state coefficient alpha
		real(kind=real_kind), parameter :: d_a(7) = (/ -1.36471e-1_real_kind, &
			4.68181e-2_real_kind, 8.07004e-1_real_kind, -7.45353e-3_real_kind, &
			-2.94418e-3_real_kind, 3.43570e-5_real_kind, 3.48658e-5_real_kind/)
		! seawater equation of state coefficient beta
		real(kind=real_kind), parameter :: d_b(7) = (/ 5.06423e-1_real_kind, &
			-3.57109e-3_real_kind, -8.76148e-4_real_kind, 5.25243e-5_real_kind, &
			1.57976e-5_real_kind, -3.46686e-7_real_kind, -1.68764e-7_real_kind/)
		! seawater equation of state coefficient gamma
		real(kind=real_kind), parameter :: d_g(7) = (/ -5.52640e-4_real_kind, &
			4.88584e-6_real_kind, 9.96027e-7_real_kind, -7.25139e-8_real_kind, &
			-3.98736e-9_real_kind, 4.00631e-10_real_kind, 8.26368e-11_real_kind/)

		real(kind=dbl_kind), parameter :: logn_multiplier9(9) = &
			(/0.102_dbl_kind, 0.272_dbl_kind, 0.427_dbl_kind, 0.532_dbl_kind, &
			0.721_dbl_kind, 0.952_dbl_kind, 1.31_dbl_kind, 1.74_dbl_kind, &
			3.31_dbl_kind/)

		real(kind=dbl_kind), parameter :: logn_multiplier6(6) = &
			(/0.145_dbl_kind, 0.385_dbl_kind, 0.585_dbl_kind, 0.860_dbl_kind, &
			1.345_dbl_kind, 2.70_dbl_kind/)
		real(kind=dbl_kind), parameter :: logn_multiplier5(5) = &
			(/1.82e-1_dbl_kind, 4.08e-1_dbl_kind, 6.96e-1_dbl_kind, &
			1.14_dbl_kind, 2.59_dbl_kind/)

		real(kind=dbl_kind), target :: gauss_bins_3(4) = (/0.0_dbl_kind, &
			1.3717e-01_dbl_kind, 3.0777e-01_dbl_kind, 1.0_dbl_kind/)
		real(kind=dbl_kind), target :: gauss_bins_4(5) = (/0.0_dbl_kind, &
			1.0153e-01_dbl_kind, 2.1483e-01_dbl_kind, 3.6633e-01_dbl_kind, &
			1.0_dbl_kind/)
		real(kind=dbl_kind), target :: gauss_bins_5(6) = (/0.0_dbl_kind, &
			8.0522e-02_dbl_kind, 1.6677e-01_dbl_kind, 2.6766e-01_dbl_kind, &
			4.0738e-01_dbl_kind, 1.0_dbl_kind/)
		real(kind=dbl_kind), target :: gauss_bins_6(7) = (/0.0_dbl_kind, &
			6.7155e-02_dbl_kind, 1.3749e-01_dbl_kind, 2.1515e-01_dbl_kind, &
			3.0840e-01_dbl_kind, 4.4080e-01_dbl_kind, 1.0_dbl_kind/)
		real(kind=dbl_kind), target :: gauss_bins_7(8) = (/0.0_dbl_kind, &
			5.7288e-02_dbl_kind, 1.1649e-01_dbl_kind, 1.8014e-01_dbl_kind, &
			2.5207e-01_dbl_kind, 3.3991e-01_dbl_kind, 4.6626e-01_dbl_kind, &
			1.0_dbl_kind/)
		real(kind=dbl_kind), target :: gauss_bins_8(9) = (/0.0_dbl_kind, &
			5.0286e-02_dbl_kind, 1.0185e-01_dbl_kind, 1.5595e-01_dbl_kind, &
			2.1515e-01_dbl_kind, 2.8294e-01_dbl_kind, 3.6696e-01_dbl_kind, &
			4.8950e-01_dbl_kind, 1.0_dbl_kind/)
		real(kind=dbl_kind), target :: gauss_bins_9(10) = (/0.0_dbl_kind, &
			4.4558e-02_dbl_kind, 9.0070e-02_dbl_kind, 1.3749e-01_dbl_kind, &
			1.8810e-01_dbl_kind, 2.4411e-01_dbl_kind, 3.0872e-01_dbl_kind, &
			3.8956e-01_dbl_kind, 5.0859e-01_dbl_kind, 1.0_dbl_kind/)
		real(kind=dbl_kind), target :: gauss_bins_10(11) = (/0.0_dbl_kind, &
			4.0102e-02_dbl_kind, 8.0840e-02_dbl_kind, 1.2285e-01_dbl_kind, &
			1.6709e-01_dbl_kind, 2.1483e-01_dbl_kind, 2.6798e-01_dbl_kind, &
			3.3004e-01_dbl_kind, 4.0802e-01_dbl_kind, 5.2355e-01_dbl_kind, &
			1.0_dbl_kind/)

		real(kind=dbl_kind), parameter :: pur_c = 0.0246059207996924_dbl_kind
		real(kind=dbl_kind), parameter :: one_over_pur_c = 40.6406250000000_dbl_kind          
		real(kind=dbl_kind), parameter :: pur_c_2 = 0.0244140625000000_dbl_kind
		real(kind=dbl_kind), parameter :: one_over_pur_c_2 = 40.96_dbl_kind
!          integer (kind=1), parameter :: pur_0 = -127
		integer (kind=2), parameter :: pur_0 = -32767
		real(kind=dbl_kind), parameter :: af_min = 5.0e-9_dbl_kind
		
		character*20, save :: version_string = 'beta - initial      '

		real(kind=dbl_kind), parameter :: c9999 = 9999.0_dbl_kind
		real(kind=dbl_kind), parameter :: c_5 = 0.5_dbl_kind
		real(kind=dbl_kind), parameter :: cn_5 = -0.5_dbl_kind      
		real(kind=dbl_kind), parameter :: c_1 = 0.1_dbl_kind
		real(kind=dbl_kind), parameter :: c_01 = 0.01_dbl_kind
		real(kind=dbl_kind), parameter :: c_001 = 0.001_dbl_kind
		real(kind=dbl_kind), parameter :: c1_5 = 3.0_dbl_kind/2.0_dbl_kind
		real(kind=dbl_kind), parameter :: c_333 = 2.0_dbl_kind/3.0_dbl_kind
		real(kind=dbl_kind), parameter :: c_111 = 1.0_dbl_kind/9.0_dbl_kind
    real(kind=dbl_kind), parameter :: c_1e3 = 1.0e-3_dbl_kind       
    real(kind=dbl_kind), parameter :: c_1e6 = 1.0e-6_dbl_kind       

		real(kind=dbl_kind), parameter :: d0_ = 0.0_dbl_kind
		real(kind=dbl_kind), parameter :: c0 = 0.0_dbl_kind
    real(kind=dbl_kind), parameter :: c1 = 1.0_dbl_kind    
    real(kind=dbl_kind), parameter :: c2 = 2.0_dbl_kind   
    real(kind=dbl_kind), parameter :: c3 = 3.0_dbl_kind           
    real(kind=dbl_kind), parameter :: c4 = 4.0_dbl_kind       
    real(kind=dbl_kind), parameter :: c5 = 5.0_dbl_kind       
    real(kind=dbl_kind), parameter :: c6 = 6.0_dbl_kind       
    real(kind=dbl_kind), parameter :: c7 = 7.0_dbl_kind       
    real(kind=dbl_kind), parameter :: c8 = 8.0_dbl_kind       
    real(kind=dbl_kind), parameter :: c9 = 9.0_dbl_kind       
    real(kind=dbl_kind), parameter :: c10 = 10.0_dbl_kind       
    real(kind=dbl_kind), parameter :: c100 = 100.0_dbl_kind       
    real(kind=dbl_kind), parameter :: c1000 = 1000.0_dbl_kind       
    real(kind=dbl_kind), parameter :: c1e3 = 1.e3_dbl_kind       
    real(kind=dbl_kind), parameter :: c1e5 = 1.e5_dbl_kind       
    real(kind=dbl_kind), parameter :: c1e6 = 1.e6_dbl_kind       

    real(kind=dbl_kind), parameter :: nc1 = -1.0_dbl_kind    
			
end module sia2_constants


