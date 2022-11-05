list of shit for new ice heat flux

call sia2_env_atmo

call sia2_conductivity
	
	<setup newton-raphson iteration>

	call sia2_heat_solver
	
	<check/repeat>
	
<update boundary fluxes/ice grwowth rates>

<update ice state>

<find ice boundary adjustments for grid changes>

call sia2_boundaries
		r_depth = r_depth + f_depth
<make sure we are not growing too large/small h_min/h_max, deal with loss of biomass. ice etc>

call sia2_new_grid (for ice)

call sia2_regrid_ice

<deal with skeletal shit>

call sia2_regrid_snow
