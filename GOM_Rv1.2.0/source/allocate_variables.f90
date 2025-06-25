!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine allocates variables (not all but some of them), 
!! which were defined in mod_global_variables.f90, and initialize those variables.
!! Other variables will be allocated when (or where) it is required.
!! 
subroutine allocate_variables
	use mod_global_variables
	
	! local variables =========================================================!
	implicit none
	! end of local variables ==================================================!
	
	!=== C0 ~ C19: maxnod, maxele, maxlayer, maxface ==========================!
   allocate(top_layer_at_node(maxnod),										&	! (maxnod)
   &			bottom_layer_at_node(maxnod),									&	! (maxnod)
   &			wetdry_node(maxnod))												! (maxnod)
   top_layer_at_node 	= 0
   bottom_layer_at_node = 0
   wetdry_node 			= 0
	
   allocate(Etide_species_coef_at_node(0:2,maxnod),					&	! (0:2,maxnod)
	&			Etide_species_coef_at_element(0:2,maxele))					! (0:2,maxele)
   Etide_species_coef_at_node 	= 0.0_dp
	Etide_species_coef_at_element = 0.0_dp
		
   allocate(ob_element_flag(maxele),										&	! (maxele)
   &			Qb_element_flag(maxele),										&	! (maxele)
   &			top_layer_at_element(0:maxele), 								&	! (0:maxele)
   &			bottom_layer_at_element(0:maxele))								! (0:maxele)
   ob_element_flag 			= 0
   Qb_element_flag			= 0
   top_layer_at_element 	= 0
   bottom_layer_at_element = 0
   
   
	
   
   allocate(eta_cell(maxele),													&	! (maxele)
   &			eta_cell_new(maxele),											&	! (maxele)
   &			eta_node(maxnod))														! (maxnod)
   eta_cell 		= 0.0_dp
   eta_cell_new 	= 0.0_dp
   eta_node 		= 0.0_dp


   
	allocate(num_sub_elm_iteration(maxlayer,maxnod))						! (maxlayer,maxnod)
	num_sub_elm_iteration = 0
	
   allocate(dz_node(maxlayer,maxnod),							&	! (maxlayer,maxnod)
   &			dz_node_new(maxlayer,maxnod),						&	! (maxlayer,maxnod)
	&			dzhalf_node(0:maxlayer,maxnod),					&	! (0:maxlayer,maxnod), actual size = (1:maxlayer-1,maxnod)
	&			dzhalf_node_new(0:maxlayer,maxnod),				&	! (0:maxlayer,maxnod), actual size = (1:maxlayer-1,maxnod)
	&			u_node(0:maxlayer,maxnod),							&	! (0:maxlayer,maxnod)
	&			v_node(0:maxlayer,maxnod),							&	! (0:maxlayer,maxnod)
	&			w_node(0:maxlayer,maxnod))								! (0:maxlayer,maxnod)
	dz_node 				= 0.0_dp
	dz_node_new 		= 0.0_dp
	dzhalf_node 		= 0.0_dp
	dzhalf_node_new 	= 0.0_dp
	u_node 				= 0.0_dp
	v_node 				= 0.0_dp
	w_node 				= 0.0_dp
	
	
   allocate(dz_cell(maxlayer,maxele),						&	! (maxlayer,maxele)
   &			dzhalf_cell(0:maxlayer,maxele),				&	! (0:maxlayer,maxele), actual size = (1:maxlayer-1,maxele)
   &			dz_cell_new(maxlayer,maxele),					&	! (maxlayer,maxele)
   &			dzhalf_cell_new(0:maxlayer,maxele),			&	! (0:maxlayer,maxele), actual size = (1:maxlayer-1,maxele)
   &			wn_cell(0:maxlayer,maxele),					&	! (0:maxlayer,maxele)
   &			wn_cell_new(0:maxlayer,maxele))					! (0:maxlayer,maxele)
   dz_cell 		= 0.0_dp
   dzhalf_cell = 0.0_dp
   dz_cell_new = 0.0_dp
   dzhalf_cell_new = 0.0_dp
   wn_cell 		= 0.0_dp
   wn_cell_new 	= 0.0_dp
   
   allocate(delta_z(maxlayer),												&	! (maxlayer), delta_z is the globalized vertical layer depth (see the difference with dz_face)
   &			z_level(0:maxlayer))													! (maxlayer)
   delta_z = 0.0_dp
   z_level = 0.0_dp

	allocate(top_layer_at_face(maxface),									& 	! (maxface)
	&			bottom_layer_at_face(maxface),								&	! (maxface)	
	&			boundary_type_of_face(maxface),								&	! (maxface)
	&			isflowside(maxface),												&	! (maxface)
	&			isflowside2(maxface),											&	! (maxface)
	&			isflowside3(maxface),											&	! (maxface)
	&			coriolis_factor(maxface))											! (maxface)
	top_layer_at_face 		= 0
	bottom_layer_at_face 	= 0
	boundary_type_of_face 	= 0
	isflowside 					= 0
	isflowside2 				= 0
	isflowside3 				= 0
	coriolis_factor 			= 0.0_dp
 
   ! 2D variables
   allocate(dz_face(maxlayer,maxface),					&	! (maxlayer,maxface), dz_face is vertical layer depth at each side (this, not delta_z, will be used for actual calculation)
   &			dz_face_new(maxlayer,maxface),			&	! (maxlayer,maxface)
   &			dzhalf_face(0:maxlayer,maxface),			&	! (0:maxlayer,maxface)
   &			dzhalf_face_new(0:maxlayer,maxface),	&	! (0:maxlayer,maxface)
   &			un_face(maxlayer,maxface),					&	! (maxlayer,maxface), normal velocity at face
   &			vn_face(maxlayer,maxface),					&	! (maxlayer,maxface), tangent velocity at face
   &			un_face_new(maxlayer,maxface),			&	! (maxlayer,maxface), normal velocity at face
   &			vn_face_new(maxlayer,maxface),			&	! (maxlayer,maxface), tangent velocity at face
   &			Av(0:maxlayer,maxface),						&	! (0:maxlayer,maxface)
   &			Kv(0:maxlayer,maxele),						&	! (0:maxlayer,maxele)
!   &			Kv_f(0:maxlayer,maxface),					&	! (0:maxlayer,maxface)
!   &			Kv_n(0:maxlayer,maxnod),					&	! (0:maxlayer,maxnod)
   &			Kh(maxlayer,maxface))							! (maxlayer,maxface)
   dz_face 				= 0.0_dp
   dz_face_new			= 0.0_dp
   dzhalf_face_new	= 0.0_dp
   dzhalf_face 		= 0.0_dp
   un_face 				= 0.0_dp
   vn_face 				= 0.0_dp
   un_face_new 		= 0.0_dp
   vn_face_new			= 0.0_dp
   Av 					= 0.0_dp
   Kv 					= 0.0_dp
!   Kv_f					= 0.0_dp
!   Kv_n					= 0.0_dp
	Kh						= 0.0_dp
	
   allocate(bed_elev(maxele),	&	! (maxele)
   &			lon_cell(maxele), &	! (maxele)
   &			lat_cell(maxele))		! (maxele)
   bed_elev = 0.0_dp
   lon_cell = 0.0_dp
   lat_cell = 0.0_dp
	

   allocate(AinvDeltaZ1(maxlayer,maxface), 	&
   &			AinvDeltaZ2(maxlayer,maxface),	&
   &			AinvG1(maxlayer,maxface), 			&
   &			AinvG2(maxlayer,maxface))
   AinvDeltaZ1 = 0.0_dp
   AinvDeltaZ2 = 0.0_dp
   AinvG1 		= 0.0_dp
   AinvG2 		= 0.0_dp
   
   
   allocate(un_ELM(maxlayer,maxface), &
   &			vn_ELM(maxlayer,maxface))
   un_ELM = 0.0_dp
   vn_ELM = 0.0_dp

   allocate(wind_stress_normal(maxface), &
   &			wind_stress_tangnt(maxface))
   wind_stress_normal = 0.0_dp
   wind_stress_tangnt = 0.0_dp



   ! transport model variables -----------------------------------------------!
   ! jw, check variable size again. For example (maxlayer+1)
   allocate(air_temperature1(maxnod), 					&
   &			air_temperature2(maxnod),					&
   &			srad(maxnod), 									&
   &			hradu(maxnod), 								&
   &			hradd(maxnod), 								&
   &			shum1(maxnod), 								&
   &			shum2(maxnod), 								&
   &			sflux(maxnod), 								&
   &			fluxsu(maxnod), 								&
   &			fluxlu(maxnod)) 
   air_temperature1 	= 0.0_dp
   air_temperature2	= 0.0_dp
   srad					= 0.0_dp
   hradu					= 0.0_dp
   hradd					= 0.0_dp
   shum1					= 0.0_dp
   shum2					= 0.0_dp
   sflux					= 0.0_dp
   fluxsu				= 0.0_dp
   fluxlu				= 0.0_dp
   
   allocate(velo_u_transport(0:maxlayer+1,maxnod), &	
   &			velo_v_transport(0:maxlayer+1,maxnod), &
   &			velo_w_transport(0:maxlayer+1,maxnod))
	velo_u_transport 	= 0.0_dp
	velo_v_transport	= 0.0_dp
	velo_w_transport	= 0.0_dp



   allocate(rho_face(maxlayer,maxface), 					&
   &			rho_cell(maxlayer,maxele), 					&
   &			rho_node(maxlayer,maxnod))					
   rho_face	= 0.0_dp
   rho_cell	= 0.0_dp
   rho_node	= 0.0_dp


	! bottom_friction_variables -----------------------------------------------!
	allocate(manning(maxface), 						&
	&			bf_height(maxface), 						&
	&			bottom_roughness(maxface),				&	
	&			bottom_drag_coefficient(maxface), 	&
   &			Gamma_B(maxface),							&
   &			Gamma_T(maxface),							&
   &			Cdb(maxface))
   manning						= 0.0_dp
   bf_height					= 0.0_dp
	bottom_roughness 			= 0.0_dp
	bottom_drag_coefficient = 0.0_dp
   Gamma_B 						= 0.0_dp
   Gamma_T						= 0.0_dp ! this is surface friction term
   Cdb 							= 0.0_dp

  	allocate(sponge_relax(maxnod), &
  	&			etaic(maxele))
  	sponge_relax 	= 0.0_dp
  	etaic 			= 0.0_dp
   
   
	! C7_1 Time series of transport variables ---------------------------------!
	allocate(salt_ser_data_num(max_salt_ser))
	salt_ser_data_num = 0
		
	allocate(salt_ser_time(max_salt_data_num,max_salt_ser), &
	&			salt_ser_salt(max_salt_data_num,max_salt_ser))
	salt_ser_time = 0.0_dp
	salt_ser_salt = 0.0_dp
	
	allocate(temp_ser_data_num(max_temp_ser))
	temp_ser_data_num = 0
		
	allocate(temp_ser_time(max_temp_data_num,max_temp_ser), &
	&			temp_ser_temp(max_temp_data_num,max_temp_ser))
	temp_ser_time = 0.0_dp
	temp_ser_temp = 0.0_dp

	allocate(temp_ob(max_ob_element))
	temp_ob = 0.0_dp
	
	! air_ser.inp variables
	allocate(air_ser_time(max_air_data_num), 		&
	&			T_air(max_air_data_num), 				&
	&			T_dew(max_air_data_num),				&
	&			air_wind_speed(max_air_data_num),	&
	&			solar_radiation(max_air_data_num),	&
	&			rain(max_air_data_num), 				&
	&			evaporation(max_air_data_num), 		&
	&			cloud(max_air_data_num))
	air_ser_time 		= 0.0_dp
	T_air 				= 0.0_dp
	T_dew 				= 0.0_dp
	air_wind_speed 	= 0.0_dp
	solar_radiation 	= 0.0_dp
	rain 					= 0.0_dp
	evaporation 		= 0.0_dp
	cloud 				= 0.0_dp
	
	allocate(phi_n(maxele), &
	&			phi_sw(maxele))
	phi_n 	= 0.0_dp
	phi_sw 	= 0.0_dp	
	
	allocate(phi_sz(maxlayer,maxele))
	phi_sz	= 0.0_dp
	   
	!=== C20 ~ C29 ============================================================!
	! open boundary related variables -----------------------------------------!
	allocate(ob_nodes(max_ob_element,2),			&
	&			ob_eta_type(max_ob_element),			&
	&			harmonic_ser_id(max_ob_element),		&
	&			eta_ser_id(max_ob_element), 			&
	&			salt_ser_id(max_ob_element),			&
	&			temp_ser_id(max_ob_element))
	ob_nodes 			= 0
	ob_eta_type 		= 0
	harmonic_ser_id 	= 0
	eta_ser_id 			= 0
	salt_ser_id			= 0
	temp_ser_id			= 0
		
   allocate(ob_cell_id(max_ob_element),			&
   &			ob_face_id(max_ob_element))
   ob_cell_id = 0
   ob_face_id = 0
	
	! allocate harmonic_ser.inp related variables -----------------------------!
   allocate(no_tidal_constituent(num_harmonic_ser),									&
   &			tidal_phase_shift(num_harmonic_ser),										&
   &			tidal_amplitude(max_no_tidal_constituent,num_harmonic_ser), 		&
   &			tidal_phase(max_no_tidal_constituent,num_harmonic_ser), 				&
   &			tidal_period(max_no_tidal_constituent,num_harmonic_ser), 			&
   &			tidal_nodal_factor(max_no_tidal_constituent,num_harmonic_ser), 	&
   &			equilibrium_argument(max_no_tidal_constituent,num_harmonic_ser))
   no_tidal_constituent 				= 0
   tidal_phase_shift						= 0.0_dp
   tidal_amplitude 						= 0.0_dp
   tidal_phase 							= 0.0_dp
   tidal_period 							= 0.0_dp
   tidal_nodal_factor					= 0.0_dp
   equilibrium_argument					= 0.0_dp

	! allocate eta_ser.inp related variables ----------------------------------!	
	allocate(eta_ser_data_num(max_eta_ser))
	eta_ser_data_num = 0
	allocate(eta_ser_time(max_eta_data_num,max_eta_ser), &
	&			eta_ser_eta(max_eta_data_num,max_eta_ser))
	eta_ser_time = 0.0_dp
	eta_ser_eta = 0.0_dp
	
	! allocate Holland storm related variables --------------------------------!
	allocate(eta_from_Holland_at_ob(max_ob_element))
	eta_from_Holland_at_ob = 0.0_dp
	
! 	allocate(num_ob_element(max_ob_segment,max_ob_element), &
! 	&			no_ob_elements(max_ob_segment))
! 	num_ob_element = 0
! 	no_ob_elements = 0
   
! 	allocate(num_serial_ob_node(max_ob_segment), &
! 	&			ob_facenum(max_ob_segment,max_ob_face))
! 	num_serial_ob_node 	= 0
! 	ob_facenum 				= 0

! 	allocate(i_tidal_boundary_temperature_type(max_ob_segment), &
! 	&			i_tidal_boundary_salinity_type(max_ob_segment))
! 	i_tidal_boundary_temperature_type 	= 0
! 	i_tidal_boundary_salinity_type		= 0
	

!   allocate(tide_amplitude_at_ob(max_ob_segment,max_no_tidal_constituent,max_ob_element), &
!   &			tidal_phase_lag_degree(max_ob_segment,max_no_tidal_constituent,max_ob_element))
!   tide_amplitude_at_ob 	= 0.0_dp
!   tidal_phase_lag_degree 	= 0.0_dp
   


   allocate(tidal_constituent_name_at_ob(max_no_tidal_constituent), & 	! character
   &			tide_name(max_no_tidal_constituent))								! character
      



! 	allocate(tide_Q(max_ob_segment),	&
! 	&			tth(max_ob_segment), 		&
! 	&			sth(max_ob_segment), 		&
! 	&			ath(max_ob_segment))		
! 	tide_Q 	= 0.0_dp
! 	tth 		= 0.0_dp
! 	sth 		= 0.0_dp
! 	ath 		= 0.0_dp

	allocate(eta_ob(max_ob_element))
	eta_ob = 0.0_dp
	
	allocate(eta_at_ob_old(max_ob_element), &
	&			eta_at_ob_new(max_ob_element))
	eta_at_ob_old = 0.0_dp
	eta_at_ob_new = 0.0_dp

   allocate(Etide_species(max_Etide_species))
   Etide_species = 0
   
   allocate(Etide_amplitude(max_Etide_species)		, &
   &			Etide_frequency(max_Etide_species)		, &
   &			Etide_nodal_factor(max_Etide_species)	, &
   &			Etide_astro_arg_degree(max_Etide_species))
   Etide_amplitude 			= 0.0_dp
   Etide_frequency 			= 0.0_dp
   Etide_nodal_factor 		= 0.0_dp
   Etide_astro_arg_degree 	= 0.0_dp
	
	
	! river boundary variables ------------------------------------------------!
	allocate(Q_nodes(max_Q_bc,2), 	&
	&			Q_ser_id(max_Q_bc), 		&
	&			Q_boundary(max_Q_bc,2), &
	&			Q_salt_ser_id(max_Q_bc),&
	&			Q_temp_ser_id(max_Q_bc))
	Q_nodes 			= 0
	Q_ser_id 		= 0
	Q_boundary 		= 0
	Q_salt_ser_id 	= 0
	Q_temp_ser_id 	= 0
	

	allocate(q_data_num(max_Q_ser))
	q_data_num 	= 0

	allocate(q_ser_time(max_q_data_num,max_Q_ser), &
	&			q_ser_Q(max_q_data_num,max_Q_ser))
	q_ser_time 	= 0.0_dp
	q_ser_Q 		= 0.0_dp


	allocate(Qu_boundary(max_Q_bc), &
	&			Qv_boundary(max_Q_bc))
	Qu_boundary = 0.0
	Qv_boundary = 0.0

	
	allocate(Q_portion(max_Q_bc), &
	&			Q_add(max_Q_bc))
	Q_portion 	= 0.0_dp	
	Q_add 		= 0.0_dp
	
   allocate(u_boundary(num_ob_cell), &
   &			v_boundary(num_ob_cell))
   u_boundary = 0.0_dp
   v_boundary = 0.0_dp
   
	allocate(salt_at_obck(maxlayer,num_ob_cell), &
	&			temp_at_obck(maxlayer,num_ob_cell), &
	&			salt_at_Qbc(num_Qb_cell), 					&
	&			temp_at_Qbc(num_Qb_cell))
	salt_at_obck = 0.0_dp
	temp_at_obck = 0.0_dp
	salt_at_Qbc = 0.0_dp
	temp_at_Qbc = 0.0_dp
	
	! Withdrawl/return & Source/Sink boundary variables -----------------------!
	allocate(WR_nodes(max_WR_bc,2), 		&
	&			WR_layer(max_WR_bc),			&
	&			WR_Q_ser_id(max_WR_bc),		&
	&			WR_portion(max_WR_bc),		&
	&			WR_salt_ser_id(max_WR_bc),	&
	&			WR_temp_ser_id(max_WR_bc))
	WR_nodes = 0
	WR_layer = 0
	WR_Q_ser_id = 0
	WR_portion = 0.0_dp
	WR_salt_ser_id = 0
	WR_temp_ser_id = 0

	allocate(WR_boundary(max_WR_bc,2),	&
	&			isflowside4(maxface),		&
	&			WR_element_flag(maxele))
	WR_boundary = 0
	isflowside4 = 0
	WR_element_flag = 0

	allocate(WRu_boundary(max_WR_bc), &
	&			WRv_boundary(max_WR_bc))
	WRu_boundary = 0.0
	WRv_boundary = 0.0
	allocate(Q_add_WR(max_WR_bc))
	Q_add_WR		= 0.0_dp

	allocate(SS_cell(max_SS_bc), 			&
	&			SS_layer(max_SS_bc),			&
	&			SS_Q_ser_id(max_SS_bc),		&
	&			SS_portion(max_SS_bc),		&
	&			SS_salt_ser_id(max_SS_bc),	&
	&			SS_temp_ser_id(max_SS_bc))
	SS_cell = 0
	SS_layer = 0
	SS_Q_ser_id = 0
	SS_portion = 0.0_dp
	SS_salt_ser_id = 0
	SS_temp_ser_id = 0
	! allocate(SSw_boundary(max_SS_bc))
	! SSw_boundary = 0.0_dp
	allocate(Q_add_SS(max_SS_bc))
	Q_add_SS = 0.0_dp
	
	allocate(SS_element_flag(maxele))
	SS_element_flag = 0
	

	!=== C30 ~ C39 ============================================================!
	! Air/Sea Boundary Conditions ---------------------------------------------!
	allocate(windp_ser_data_num(max_windp_data_num))
	windp_ser_data_num = 0
	allocate(windp_station_node(max_windp_station))
	windp_station_node = 0

	allocate(windp_ser_time(max_windp_data_num,max_windp_station),	&
	&			windp_u(max_windp_data_num,max_windp_station),			&
	&			windp_v(max_windp_data_num,max_windp_station),			&
	&			windp_p(max_windp_data_num,max_windp_station))
	windp_ser_time = 0.0_dp
	windp_u = 0.0_dp
	windp_v = 0.0_dp
	windp_p = 0.0_dp
   
	allocate(wind_u_at_node(maxnod),			&
	&			wind_v_at_node(maxnod),			&
   &			wind_u_at_face(maxface),		&
   &			wind_v_at_face(maxface))
   wind_u_at_node = 0.0_dp
   wind_v_at_node = 0.0_dp
   wind_u_at_face = 0.0_dp
   wind_v_at_face = 0.0_dp

	allocate(airp_at_node(maxnod), &
	&			airp_at_face(maxface))
	airp_at_node = 0.0_dp
	airp_at_face = 0.0_dp


  	! air_p_variables ---------------------------------------------------------!
   allocate(air_p0(maxnod),													&	! (maxnod)
   &			air_p1(maxnod),													&	! (maxnod)
   &			air_p2(maxnod))														! (maxnod)
   air_p0 = 0.0_dp
   air_p1 = 0.0_dp
   air_p2 = 0.0_dp
   
   ! wind_variables ----------------------------------------------------------!
   allocate(wind_u0(maxnod), 													&	! (maxnod)
   &			wind_u1(maxnod), 													&	! (maxnod)
   &			wind_u2(maxnod),													&	! (maxnod)
   &			wind_v0(maxnod),													&	! (maxnod)
   &			wind_v1(maxnod), 													&	! (maxnod)
   &			wind_v2(maxnod))														! (maxnod)
   wind_u0 = 0.0_dp
   wind_v0 = 0.0_dp
   wind_u1 = 0.0_dp
   wind_u2 = 0.0_dp
   wind_v1 = 0.0_dp
   wind_v2 = 0.0_dp
	
   

	!=== C40 ~ C49 ============================================================!
	! Output Control Options --------------------------------------------------!
	allocate(eta_from_bottom(max(1,tser_station_num)))
	eta_from_bottom = 0.0_dp
	
   allocate(tser_station_cell(max_tser_station), &
   &			tser_station_node(max_tser_station))
   tser_station_cell = 0
   tser_station_node = 0
   
   allocate(tser_station_name(max_tser_station)) ! character
   allocate(max_eta_node(maxnod), 	&
   &			max_flood_time(maxnod), &
   &			flood_id(maxnod))
   max_eta_node 	= 0.0_dp
   max_flood_time 	= 0.0_dp
   flood_id 			= 0
   
   ! vertically averaged variables
   allocate(ubar_node(maxnod),												&
   &			vbar_node(maxnod),												&
   &			sbar_node(maxnod),												&
   &			tbar_node(maxnod),												&
   &			rbar_node(maxnod))
   ubar_node = 0.0_dp
   vbar_node = 0.0_dp
   sbar_node = 0.0_dp
   tbar_node = 0.0_dp
   rbar_node = 0.0_dp
	
	
	! define netcdf output variable names here
	netcdf_variable_name_cell(1) = 'eta_cell'
	netcdf_variable_name_cell(2) = 'wn_cell'
	netcdf_variable_name_cell(3) = 'salt_cell'
	netcdf_variable_name_cell(4) = 'temp_cell'
	netcdf_variable_name_cell(5) = 'rho_cell'
	netcdf_variable_name_cell(6) = 'Kv'

	! define netcdf output variable names here
	netcdf_variable_name_node(1) = 'eta_node'
	netcdf_variable_name_node(2) = 'u_node'
	netcdf_variable_name_node(3) = 'v_node'
	netcdf_variable_name_node(4) = 'w_node'
	netcdf_variable_name_node(5) = 'salt_node'
	netcdf_variable_name_node(6) = 'temp_node'
	netcdf_variable_name_node(7) = 'rho_node'
	netcdf_variable_name_node(8) = 'wind_u_at_node'
	netcdf_variable_name_node(9) = 'wind_v_at_node'
	netcdf_variable_name_node(10) = 'airp_at_node'
	netcdf_variable_name_node(11) = 'ubar_node'
	netcdf_variable_name_node(12) = 'vbar_node'
	netcdf_variable_name_node(13) = 'sbar_node'
	netcdf_variable_name_node(14) = 'tbar_node'
	netcdf_variable_name_node(15) = 'rbar_node'

	! define netcdf output variable names here
	netcdf_variable_name_face(1) = 'Av'
	netcdf_variable_name_face(2) = 'Kh'

	!===	END =================================================================!
	
	!=== Transport variables ==================================================!
	allocate(salt_cell(maxlayer,maxele),		&
	&			salt_cell_new(maxlayer,maxele),	&
	&			salt_node(maxlayer,maxnod),		&
	&			salt_face(maxlayer,maxface), 		&
	&			temp_cell(maxlayer,maxele),		&
	&			temp_cell_new(maxlayer,maxele),	&
	&			temp_node(maxlayer,maxnod),		&
	&			temp_face(maxlayer,maxface))
	salt_cell = 0.0_dp
	salt_cell_new = 0.0_dp
	salt_node = 0.0_dp
	salt_face = 0.0_dp
	
	! initialize to 20 degrees for safety:
	! this will be updated with either ref_temp (in C6-3) or from temp.ser
	temp_cell = 20.0_dp
	temp_cell_new = 20.0_dp
	temp_node = 20.0_dp
	temp_face = 20.0_dp

! 	allocate(con_cell(maxlayer,maxele,maxtran2),			&
! 	&			con_cell_new(maxlayer,maxele,maxtran2),	&
! 	&			con_node_new(maxlayer,maxnod,maxtran2),	&
! 	&			con_face_new(maxlayer,maxface,maxtran2))
! 	allocate(r_jk(maxlayer,maxface,maxtran2), &
! 	&			r_ik(0:maxlayer,maxele,maxtran2))
! 	con_cell = 0.0_dp
! 	con_cell_new = 0.0_dp
! 	con_node_new = 0.0_dp
! 	con_face_new = 0.0_dp
! 	r_jk = 0.0_dp ! set to upwind first
! 	r_ik = 0.0_dp ! set to upwind first
	
	!=== From subroutines =====================================================!
	! calculate_velocity_at_node.f90
   allocate(u_face_level(0:maxlayer,maxface), &
   &			v_face_level(0:maxlayer,maxface))
   u_face_level = 0.0_dp
   v_face_level = 0.0_dp





	allocate(tem0(maxlayer,maxnod), &
	&			sal0(maxlayer,maxnod), &
	&			q2(0:maxlayer,maxface), &
	&			xl(0:maxlayer,maxface), &
! 	&			velo_u_transport(0:maxlayer,maxnod), &
! 	&			velo_v_transport(0:maxlayer,maxnod), &
! 	&			velo_w_transport(0:maxlayer,maxnod), &
	&			srho(maxlayer,maxface), &
	&			erho(maxlayer,maxele), &
	&			prho(maxlayer,maxnod))

end subroutine allocate_variables