module sia2_common
	
	use sia2_parameters

	public
	
		real(kind=dbl_kind), allocatable :: ida_multiplier(:)
		real(kind=dbl_kind), allocatable :: lda_multiplier(:)
		real(kind=dbl_kind), allocatable :: sda_multiplier(:)

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

end module sia2_common