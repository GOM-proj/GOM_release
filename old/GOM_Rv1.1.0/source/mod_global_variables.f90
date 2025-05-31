!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! Collection of Global Variables
!!
module mod_global_variables
	implicit none
	!! define precision control parameters ===================================!!
	integer, parameter :: sp = kind(0.0) 				! jw
	integer, parameter :: dp = kind(0.d0) 				! jw
	! jw
	
	!! === cpu time variables ================================================!!
	real(dp),save :: start, finish						! jw
	real(dp),save :: elapsed_time = 0.0_dp 			! jw
	integer, save :: gom_start_year, gom_start_month, gom_start_day, gom_start_hour, gom_start_minute, gom_start_second			! jw
	integer, save :: gom_finish_year, gom_finish_month, gom_finish_day, gom_finish_hour, gom_finish_minute, gom_finish_second	! jw
	
	!! === Parameters ========================================================!!
   real(dp),parameter :: small_06	=	1.e-6
	
	! jw
	! jw
   real(dp),parameter :: pi = 3.141592653589793_dp
   real(dp),parameter :: deg2rad = pi/180.0_dp		! jw
   real(dp),parameter :: rad2deg = 180.0_dp/pi		! jw
   
   ! jw
	! jw
	! jw

	! jw
	
   !! === main.inp variables ================================================!!
	!!========================================================================!!
	!! C0 ~ C19: General Model Configuration Information							  !!
	!!========================================================================!!

	! jw
	character(len=200), save :: project_name
	
	!! C01	Define general parameters ---------------------------------------!!
 	real(dp),save :: gravity									! jw
   real(dp),save :: rho_o										! jw
   real(dp),save :: rho_a										! jw
   integer, save :: max_no_neighbor_node					! jw

	!! C1	Read grid information -----------------------------------------------!
	integer, save :: maxnod, maxele							! jw
	real(dp),save :: h_node_adjust							! jw
	integer, save :: ana_depth									! jw
	integer, save :: coordinate_system						! jw
	real(dp),save :: lon, lat									! jw
	integer, save :: voronoi									! jw
	integer, save :: node_mirr, cell_mirr
	
	! jw
	integer, save :: utm_projection_zone					! jw
	real(dp),save :: lon_mid, lat_mid						! jw
	real(dp),save :: standard_meridian						! jw
	
	!! C2 Read vertical layer information
	integer, save :: maxlayer									! jw
	real(dp),save :: MSL											! jw
	! jw
	
	!! C2_1
   real(dp),save,allocatable,dimension(:) ::	&
   &			delta_z,									&			! jw
   &			z_level												! jw
	
	! jw
   real(dp),save :: delta_z_min								! jw
	
	!! C3 Read time related variables -----------------------------------------!
	real(dp),save :: dt											! jw
	integer, save :: ndt											! jw
	integer, save :: data_start_year							! jw
	integer, save :: data_time_shift							! jw
   integer, save :: start_year, start_month, start_day, start_hour, start_minute, start_second	! jw

	! jw
   integer, save :: it = 0										! jw
   real(dp),save :: current_jday_1900, air_p_julian_day_1900_2, wind_julian_day_1900_2   
   
   ! jw
   ! jw
   ! jw
   ! jw
   ! jw
   real(dp),save :: julian_day , jday, jday_1900		
   real(dp),save :: elapsed_time_hr
   integer, save :: year, month, day, hour, minute
   integer, save :: reference_diff_days
   
	!! C4	Restart & screen show control ---------------------------------------!
	integer, save :: restart_in, restart_out, restart_freq, ishow, ishow_frequency

   !! C5_01 Momentum equation variables I - Propagation/Nonlinear advection ---!
   real(dp),save :: theta
   integer, save :: advection_flag
   real(dp),save :: adv_onoff_depth
   integer, save :: ELM_backtrace_flag,	ELM_sub_iter, ELM_min_iter, ELM_max_iter

	!! C5_02	Momentum equation variables II - Bottom friction
   integer, save :: bf_flag
   integer, save :: bf_varying
   real(dp),save,allocatable,dimension(:) :: manning		! jw
   real(dp),save,allocatable,dimension(:) :: bf_height	! jw
   real(dp),save :: von_Karman									! jw
	
	!! C5_03	Momentum equation variables III - Baroclinic gradient/Coriolis
   integer, save :: baroclinic_flag
   integer, save :: ana_density
   real(dp),save :: ref_salt, ref_temp					! jw
   real(dp),save :: dry_depth
   integer, save :: Coriolis_option
   real(dp),save :: Coriolis								! jw
	! jw
   real(dp),save,allocatable,dimension(:) :: 	&
   &			coriolis_factor								! jw
	
	!! C5_04 Momentum equation variables IV - Vertical & horizontal diffusion --!
   real(dp),save :: Ah_0									! jw
   real(dp),save :: Kh_0									! jw
   real(dp),save :: Smagorinsky_parameter
	real(dp),save,allocatable,dimension(:,:) :: &	! jw
	&			Kh													! jw
	
	real(dp),save :: Av_0									! jw
	real(dp),save :: Kv_0									! jw
   real(dp),save,allocatable,dimension(:,:) ::	&
   ! jw
   ! jw
   ! jw
   ! jw
   ! jw
   ! jw
   ! jw
   ! jw
   ! jw
   &			Av													! jw
   
	real(dp),save,allocatable,dimension(:,:) :: 	&
	&			Kv,											&	! jw
	&			Kv_f,											&	! jw
	&			Kv_n												! jw

	!! C6		Matrix solver option: OpneMP version of the Jacobi Preconditioned Conjugate Gradient method:
	integer, save :: max_iteration_pcg
	real(dp),save :: error_tolerance
	integer, save :: pcg_result_show



	!! C7		off/on tranport equations & selection of the solver -----------------!	
	integer,save :: transport_flag, transport_solver

	!! C7_1	TVD solver -------------------------------------------------------!
	integer,save :: trans_sub_iter, h_flux_limiter, v_flux_limiter
	
	!! C7_2	Turn off/on transport variables ----------------------------------!
	! jw
	! jw
	integer,save :: is_tran(2)=0
	integer,save :: num_tran_ser(2)=0
	integer,save :: tran_ser_shape(2)=0
	
	! jw
   integer,save :: maxtran = 2							! jw
	integer,save :: maxtran2								! jw
	integer,save,allocatable,dimension(:) :: &
	&			tran_id											! jw
	
	! jw
	integer,save :: is_salt, num_salt_ser
   integer,save :: salt_ser_shape						! jw
	integer,save :: max_salt_ser, max_salt_data_num
	integer,save,allocatable,dimension(:) :: 	&
	&			salt_ser_data_num								! jw
	
	! jw
	! jw
	
	real(dp),save,allocatable,dimension(:,:) :: 	&
	&			salt_ser_time,								&	! jw
	&			salt_ser_salt									! jw
	
	! jw
	real(dp),save,allocatable,dimension(:,:) :: salt_at_obck	! jw
	real(dp),save,allocatable,dimension(:) :: salt_at_Qbc 	! jw
	
	!! C7_3	Heat transport control variables 1 -------------------------------!
	integer,save :: heat_option, sol_swr
	real(dp),save:: fWz_a, fWz_b, fWz_c, wind_height
	
	!! C7_4 	Heat transport control variables 2 -------------------------------!
	real(dp),save:: light_extinction, sol_absorb, sed_water_exchange, T_sed, sed_temp_coeff
	
	! jw
	real(dp),save,allocatable,dimension(:) :: &
	&			phi_n,  									& 	! jw
	&			phi_sw										! jw
	real(dp),save,allocatable,dimension(:,:):: &
	&			phi_sz										! jw
		
	! jw
	integer,save :: max_air_data_num
	real(dp),save,allocatable,dimension(:) ::	&
	&			air_ser_time,							&	! jw
	&			T_air,									&	! jw
	&			T_dew,									&	! jw
	&			air_wind_speed,						&	! jw
	&			solar_radiation,						&	! jw
	&			cloud,									&	! jw
	&			rain,										&	! jw
	&			evaporation									! jw
	
	! jw
	real(dp),save,allocatable,dimension(:,:) :: temp_at_obck	! jw
	real(dp),save,allocatable,dimension(:) :: temp_at_Qbc 	! jw

	! jw
	integer,save :: is_temp, num_temp_ser
   integer,save :: temp_ser_shape						! jw
	integer,save :: max_temp_ser, max_temp_data_num
	integer,save,allocatable,dimension(:) :: 	&
	&			temp_ser_data_num								! jw
	
	real(dp),save,allocatable,dimension(:,:) :: 	&
	&			temp_ser_time,								&	! jw
	&			temp_ser_temp									! jw

   real(dp),save,allocatable,dimension(:) ::		&
   &			temp_ob											! jw

	
	!! C8		Model termination control ----------------------------------------!
	integer, save :: terminate_check_freq
	real(dp),save :: eta_min_terminate,	eta_max_terminate, uv_terminate, salt_terminate, temp_terminate
	
		      
	
   
	
	
	!!========================================================================!!
	!! C20 ~ C29, Tidal & River Boundary Conditions									  !!
	!!========================================================================!!
   ! jw
   integer, save :: num_ob_cell							! jw
   integer, save :: ob_info_option						! jw
   integer, save :: num_tide_bc							! jw
   integer, save :: num_harmonic_ser					! jw
   integer, save :: check_tide							! jw
   integer, save :: num_eta_bc							! jw
   integer, save :: num_eta_ser							! jw
   integer, save :: check_eta
   integer, save :: eta_ser_shape						! jw
   
   ! jw
   ! jw
   integer, save :: eta_ser_frequency
   integer, save :: max_eta_data_num, max_eta_ser
   real(dp),save,allocatable :: 						&
   &			eta_ser_time(:,:),						&	! jw
   &			eta_ser_eta(:,:)								! jw
	
	integer,save,allocatable ::						&
	&			eta_ser_data_num(:)							! jw
	
   real(dp),save,allocatable,dimension(:) ::		&
   &			eta_ob											! jw
   
   real(dp),save,allocatable,dimension(:) ::		&	! jw
   &			eta_at_ob_old,								&	! jw
   &			eta_at_ob_new									
   
   ! jw
	integer, save, allocatable, dimension(:,:) ::&
	&			ob_nodes 										! jw
   ! jw
   integer, save, allocatable, dimension(:) :: 	&
   &			ob_cell_id,									&	! jw
   &			ob_face_id										! jw
	
	integer, save, allocatable, dimension(:) :: 	&
	&			ob_eta_type, 								&	! jw
   &			harmonic_ser_id,							&	! jw
   &			eta_ser_id,									&	! jw
   &			salt_ser_id,								&	! jw
   &			temp_ser_id										! jw
   
   
   real(dp),save :: eta_ser_dt

   	
   integer, save :: max_ob_node							! jw
   integer, save :: max_ob_element						! jw
   integer, save :: max_ob_face							! jw

   integer, save,allocatable :: 						&
   &			num_ob_element(:,:), 					&	! jw
   &			no_ob_elements(:)								! jw
	
	! jw
   integer, save, allocatable, dimension(:) ::  &
   &			no_tidal_constituent							! jw
   integer, save :: max_no_tidal_constituent			! jw
   real(dp),save, allocatable, dimension(:) :: 	&
   &			tidal_phase_shift								! jw

   real(dp),save, allocatable, dimension(:,:)::	&
   &			tidal_amplitude,   						&	! jw
   &  		tidal_phase, 								&	! jw
   &			tidal_period,								&	! jw
   &			tidal_nodal_factor,						&	! jw
   &			equilibrium_argument							! jw


   ! jw
   integer, save :: tide_spinup, baroclinic_spinup	
   real(dp),save :: tide_spinup_period, baroclinic_spinup_period

   
   ! jw
   integer, save :: num_Qb_cell, Qbc_info_option, num_Q_ser, check_Q, Q_ser_shape
   real(dp),save, allocatable, dimension(:) :: 	&
   &			Qu_boundary, 								&	! jw
   &			Qv_boundary										! jw
   
   ! jw
   integer, save :: max_Q_bc, max_Q_ser
   
   ! jw
	integer, save, allocatable :: 					&
	&			Q_nodes(:,:), 								&	! jw
	&			Q_ser_id(:),								&	! jw
	&			Q_salt_ser_id(:),							&	! jw
	&			Q_temp_ser_id(:)								! jw
	real(dp),save, allocatable :: 					&
	&			Q_portion(:)									! jw
   
   ! jw
   integer, save, allocatable :: 					&
   &			Q_boundary(:,:)								! jw
   real(dp),save, allocatable :: 					&
   &			Q_add(:)											! jw
   real(dp),save, allocatable :: 					&
   &			q_ser_time(:,:), 							&	! jw
   &			q_ser_Q(:,:)									! jw
   integer, save, allocatable :: 					&
   &			q_data_num(:)									! jw
   
   integer, save :: q_interp_method
   integer, save :: face_id, element_id
	integer, save :: max_q_data_num 						! jw
	
	! jw
	! jw
	integer, save :: num_WR_cell

	! jw
	integer, save :: max_WR_bc
	integer, save, allocatable ::	&
	&			WR_boundary(:,:),		&	! jw
	&			isflowside4(:),		&	! jw
	&			WR_element_flag(:)		! jw
   real(dp),save, allocatable, dimension(:) :: 	&
   &			WRu_boundary, 								&	! jw
   &			WRv_boundary									! jw
   real(dp),save, allocatable :: 					&
   &			Q_add_WR(:)										! jw

	
	! jw
	integer, save, allocatable :: 	&
	&			WR_nodes(:,:),				&	! jw
	&			WR_layer(:),				&	! jw
	&			WR_Q_ser_id(:),			&	! jw
	&			WR_salt_ser_id(:),		&	! jw
	&			WR_temp_ser_id(:)				! jw
	real(dp),save, allocatable ::		&
	&			WR_portion(:)					! jw
	
		
	! jw
	integer, save :: num_SS_cell

	! jw
	integer, save :: max_SS_bc
	integer, save, allocatable :: SS_element_flag(:)	! jw
! jw
! jw
   real(dp),save, allocatable :: 						&
   &			Q_add_SS(:)											! jw

	! jw
	integer, save, allocatable ::		&
	&			SS_cell(:),					&	! jw
	&			SS_layer(:),				&	! jw
	&			SS_Q_ser_id(:),			&	! jw
	&			SS_salt_ser_id(:),		&	! jw
	&			SS_temp_ser_id(:)				! jw
	real(dp),save, allocatable ::		&
	&			SS_portion(:)					! jw
	
	
   ! jw
   integer, save :: no_Etide_species 
   real(dp),save :: Etide_cutoff_depth
   
   ! jw
   integer, save :: max_Etide_species
   
   ! jw
   integer, save, allocatable :: 					&
   &			Etide_species(:) 								! jw
   real(dp),save, allocatable, dimension(:) :: 	&
   &			Etide_amplitude, 							&	! jw
   &			Etide_frequency, 							&	! jw
   &			Etide_nodal_factor, 						&	! jw
   &			Etide_astro_arg_degree						! jw

	! jw
	real(dp),save,allocatable, dimension(:,:)::	&
	&			Etide_species_coef_at_node,			&	! jw
	&			Etide_species_coef_at_element				! jw
	   
!=============================================================================!
! jw
! jw
!=============================================================================!
   ! jw
   integer, save :: wind_flag, airp_flag, num_windp_ser, wind_formular, wind_spinup
   real(dp),save :: wind_spinup_period
	
	! jw
	! jw
	integer, save :: max_windp_data_num					! jw
	integer, save :: max_windp_station					! jw
	! jw
	integer, save,allocatable,dimension(:) ::		&
	&			windp_ser_data_num,						&	! jw
	&			windp_station_node							! jw
	real(dp),save,allocatable,dimension(:,:):: 	&
	&			windp_ser_time,							&	! jw
	&			windp_u,										&	! jw
	&			windp_v,										&	! jw
	&			windp_p											! jw
	! jw
	
	! jw
! jw
! jw
! jw
! jw
	
! jw
! jw
	
! jw
! jw
! jw
! jw
	! jw
	
	
   real(dp),save,allocatable,dimension(:)	::		&
   &			wind_u_at_node,							&	! jw
   &			wind_v_at_node,							&	! jw
   &			wind_u_at_face, 							&	! jw
   &			wind_v_at_face,							&	! jw
   &			wind_stress_normal,						&	! jw
   &			wind_stress_tangnt							! jw
	
	real(dp),save,allocatable,dimension(:) ::		&
	&			airp_at_node,								&	! jw
	&			airp_at_face									! jw
	   

   ! jw
   integer, save :: hurricane_flag															! jw
   integer, save :: hurricane_data_type, hurricane_interp_method					! jw
   integer, save :: hurricane_dt																! jw
	
	! jw
	integer, save :: hurricane_start_year, hurricane_start_month, hurricane_start_day
	integer, save :: hurricane_end_year, hurricane_end_month, hurricane_end_day
	
	! jw
	real(dp),save :: hurricane_start_jday, hurricane_end_jday
   
	! jw
   real(dp),save,allocatable,dimension(:) :: 	&
   &			wind_u0, 									&	! jw
   &			wind_u1, 									&	! jw
   &			wind_u2,										&	! jw
   &			wind_v0,										&	! jw
   &			wind_v1, 									&	! jw
   &			wind_v2,										&	! jw
   &			air_p0, 										&	! jw
   &			air_p1, 										&	! jw
   &			air_p2											! jw
	! jw
	! jw
   real(dp),save,allocatable,dimension(:) :: 	&	
   &			wind_u1_new, 								&	! jw
   &			wind_v1_new, 								&	! jw
   &			wind_u2_new, 								&	! jw
   &			wind_v2_new, 								&	! jw
   &			air_p1_new, 								&	! jw
   &			air_p2_new,									&	! jw
   &			shiftx, 										&	! jw
   &			shifty											! jw
	
	! jw
   integer, save :: &
   &			hurric_year_1, hurric_month_1, hurric_day_1, hurric_hour_1, hurric_minute_1, &
   &			hurric_year_2, hurric_month_2, hurric_day_2, hurric_hour_2, hurric_minute_2
   
   real(dp),save :: &
   &			hurric_x, hurric_y , hurric_latitude, 			&
   &  	   hurric_delta_pressure, hurric_mwr,				&
   &			hurric_x_1, hurric_y_1 , hurric_latitude_1,  &
   &  	   hurric_delta_pressure_1, hurric_mwr_1,			&
   &			hurric_x_2, hurric_y_2 , hurric_latitude_2,	&
   &  	   hurric_delta_pressure_2, hurric_mwr_2
   real(dp),save :: &
   &			hurric_julian_day_1     , hurric_julian_day_2,		&
   &  	   hurric_julian_day_1900_1, hurric_julian_day_1900_2	  
   real(dp),save :: hurric_coriolis, hurric_beta
   real(dp),save :: hurricane_time_1, hurricane_time_2
   real(dp),save :: hurricane_center_x, hurricane_center_y
	
	! jw
   integer, save :: hurricane_read_interval
	real(dp),save, allocatable :: regular_grid_xc(:), regular_grid_yc(:)
	integer, save, allocatable :: regular_grid_node_count(:,:), regular_grid(:,:,:)
	integer, save :: regular_grid_xi, regular_grid_yi ! jw
	real(dp),save:: regular_grid_half_dx, regular_grid_half_dy
	
	integer,save, allocatable :: hurricane_search_nodes(:,:), hurricane_search_count(:)

   ! jw
   integer, save :: holland_flag
   real(dp),save :: holland_start_jday, holland_end_jday

	! jw
   real(dp),save,allocatable,dimension(:) :: 	&
   &			eta_from_Holland_at_ob						! jw
   
!=============================================================================!
! jw
! jw
!=============================================================================!

   ! jw
   integer, save :: check_grid_2D, check_grid_2DO, check_grid_3D, check_grid_format, check_grid_info
   real(dp),save :: check_grid_unit_conv

   ! jw
	integer, save :: cell_connectivity_size_2D	! jw
	integer, save :: cell_connectivity_size_3D	! jw

   
   ! jw
   integer, save :: tser_station_num, tser_info_option, tser_hloc, tser_time, tser_frequency, tser_format
   
   ! jw
   real(dp),save :: tser_time_conv ! jw
   integer, save :: max_tser_station
	real(dp),save,allocatable,dimension(:) :: 	&
	&			eta_from_bottom								! jw
   
	! jw
   integer, save, allocatable, dimension(:) :: 	&
   &			tser_station_cell,						&	! jw
   &			tser_station_node								! jw
   character(25), save,allocatable, dimension(:) :: &
   &			tser_station_name								! jw

   ! jw
   integer, save :: tser_eta, tser_H, tser_u, tser_v, tser_salt, tser_temp, tser_airp
   
	! jw
	integer, save :: IS2D_switch, IS2D_time, IS2D_frequency, IS2D_start, IS2D_end, IS2D_format, IS2D_binary, IS2D_File_freq
	real(dp),save :: IS2D_unit_conv
	real(dp),save :: IS2D_time_conv 						! jw
	integer, save :: IS2D_File_num, IS2D_vtk_num
	character(len = 100), save :: IS2D_File_name
	integer, save :: IS2D_flood_map

	! jw
	real(dp),save,allocatable,dimension(:) :: max_eta_node, max_flood_time	! jw
	integer, save,allocatable,dimension(:) :: flood_id								! jw

	! jw
	integer, save :: IS2D_variable(6) = 0
	integer, save :: zone_num_2D							! jw

	
	! jw
	integer, save :: IS3D_surf_switch, IS3D_full_switch, IS3D_time, IS3D_frequency, IS3D_start, IS3D_end, &
	&					  IS3D_format, IS3D_grid_format, IS3D_binary, IS3D_File_freq
	real(dp),save :: IS3D_unit_conv
	real(dp),save :: IS3D_time_conv ! jw
	integer, save :: zone_num_3D_surf, zone_num_3D_full
	integer, save :: IS3D_File_num_surf, IS3D_File_num_full, IS3D_vtk_num
	character(len = 100), save :: IS3D_File_name_surf, IS3D_File_name_full

	! jw
	integer, save :: IS3D_variable(6) = 0
	
	! jw
	integer, save :: IS2D_dump_switch, IS2D_dump_binary, IS2D_dump_time, IS2D_dump_File_freq, &
	&					  IS2D_dump_frequency, IS2D_dump_start, IS2D_dump_end
	integer, save :: IS2D_dump_File_num
	character(len = 100), save :: IS2D_dump_File_name
	real(dp),save :: IS2D_dump_time_conv 						! jw

	! jw
	integer, save :: IS3D_dump_switch, IS3D_dump_binary, IS3D_dump_time, IS3D_dump_File_freq, &
	&					  IS3D_dump_frequency, IS3D_dump_start, IS3D_dump_end
	integer, save :: IS3D_dump_File_num
	character(len = 100), save :: IS3D_dump_File_name
	real(dp),save :: IS3D_dump_time_conv 						! jw
	
	! jw
	integer, save :: dia_advection, dia_momentum, dia_freesurface, dia_eta_at_ob
	integer, save :: dia_bottom_friction, dia_face_velocity, dia_node_velocity
   !=== End of main.inp variables ============================================!
   
   
	!==========================================================================!
   ! jw
   real(dp),save :: 	&
   &			wtime1, wtime2, wtratio, windx1_s, windx2_s, windy1_s, windy2_s,	&
   &      	spinup_function_wind
      

   ! jw
   ! jw
   real(dp),save,allocatable,dimension(:,:) ::	&
	! jw
	! jw
   &			un_ELM, 										&	! jw
   &			vn_ELM,										&	! jw
   &			un_face,										&	! jw
   &			vn_face,										&	! jw
   &			un_face_new,								&	! jw
   &			vn_face_new,								&	! jw
   ! jw
   ! jw
   &			wn_cell,										&	! jw
   &			wn_cell_new,								&	! jw
   ! jw
	&			u_node,										&	! jw
	&			v_node,										&	! jw
	&			w_node											! jw

	! jw
	! jw
	! jw
   real(dp),save,allocatable,dimension(:,:) :: 	&
   &			u_face_level,								&	! jw
   &			v_face_level									! jw
   
   real(dp),save,allocatable,dimension(:) :: 	&
   &			u_boundary, 								&	! jw
   &			v_boundary										! jw
   
   ! jw
   real(dp),save,allocatable,dimension(:,:) :: 	&
   &			velo_u_transport, 						&	! jw
   &			velo_v_transport, 						&	! jw
   &			velo_w_transport								! jw
   
   ! jw
   real(dp),save,allocatable,dimension(:) :: 	&
   &			ubar_node,									&	! jw
   &			vbar_node,									&	! jw
   &			sbar_node,									&	! jw
   &			tbar_node,									&	! jw
   &			rbar_node										! jw
   ! jw
   
   ! jw
   integer, save :: maxface								! jw
   real(dp),save :: xn_min, yn_min 						! jw
   integer, save :: start_end_node(4,4,3) 			! jw
   
   ! jw
   real(dp),save,allocatable,dimension(:) :: 	&
   ! jw
   &			x_node, 										&	! jw
   &			y_node, 										&	! jw
   &			lon_node,									&	! jw
   &			lat_node,									&	! jw
   &			h_node, 										&	! jw
   &			eta_node,									&	! jw
   ! jw
   &			x_cell,       								&	! jw
   &			y_cell,       								&	! jw
   &			lon_cell, 									&	! jw
   &			lat_cell,									&	! jw
	&			h_cell,										&	! jw
	&			bed_elev,									&	! jw
   &			eta_cell,									&	! jw
	&			eta_cell_new,								&	! jw
	! jw
   &			face_length,                 			&	! jw
   &  		x_face, 										&	! jw
   &			y_face,              					&	! jw
   &			h_face											! jw
	
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	real(dp),save,allocatable,dimension(:,:) ::	&
	! jw
	&			dz_node,										&	! jw
	&			dzhalf_node,								&	! jw
	&			dz_node_new,								&	! jw
	&			dzhalf_node_new,							&	! jw
   ! jw
   &			dz_cell,										&	! jw
   &			dzhalf_cell,								&	! jw
   &			dz_cell_new,								&	! jw
   &			dzhalf_cell_new,							&	! jw
   ! jw
   &			dz_face,										&	! jw
   &			dzhalf_face,								&	! jw
   &			dz_face_new,								&	! jw
   &			dzhalf_face_new								! jw
   																! jw
   																! jw
	! jw
   integer, save, allocatable, dimension(:) ::	&
   ! jw
   &			top_layer_at_node,						&	! jw
   &			bottom_layer_at_node,					&	! jw
	! jw
	&			top_layer_at_element, 					&	! jw
   &			bottom_layer_at_element,				&	! jw
	! jw
   &			top_layer_at_face, 						&	! jw
   &			bottom_layer_at_face							! jw

	! jw
	! jw
   integer, save, allocatable, dimension(:) ::	&
   &			adj_cells_at_node,						&	! jw
   &			adj_nodes_at_node,						&	! jw
   &			initial_wetdry_node,						&	! jw
   &			wetdry_node 									! jw
	integer, save, allocatable, dimension(:,:)::	&  
	&			adj_cellnum_at_node,	   				&	! jw
	&			adj_nodenum_at_node,   					&	! jw
	&			node_count_each_element						! jw
	! jw
   integer, save, allocatable, dimension(:,:):: &
   &			nodenum_at_cell,							& 	! jw
  	&			nodenum_at_cell_tec,						&	! jw
   &			adj_cellnum_at_cell,						&	! jw
	&			facenum_at_cell								! jw
	! jw
   integer, save,allocatable,dimension(:,:) :: 	&
   &			adj_cellnum_at_face, 					&	! jw
   &			nodenum_at_face								! jw

	! jw
	! jw
   real(dp),save,allocatable,dimension(:) :: 	&
   &			area   	                      			! jw
   integer, save, allocatable, dimension(:) ::	&
   &			tri_or_quad										! jw
	integer,save,allocatable,dimension(:,:) ::	&
   &			sign_in_outflow								! jw
	! jw
   real(dp),save,allocatable,dimension(:) :: 	&
   &			delta_j,   									&	! jw
   &			cos_theta,									&	! jw
   &			sin_theta,									&	! jw
   &			cos_theta2,									&	! jw
   &			sin_theta2										! jw
   ! jw

	! jw
	! jw
   integer, save, allocatable :: 					&
   &			num_serial_ob_node(:) 						! jw
! jw
   ! jw
   integer,save,allocatable,dimension(:) ::		&
   &			ob_element_flag,							&	! jw
   &			Qb_element_flag								! jw
   ! jw
   integer,save,allocatable,dimension(:) :: 		&
   &			boundary_type_of_face,					&	! jw
   &			isflowside,		 							&	! jw
   &			isflowside2,								&	! jw
   &			isflowside3										! jw


   integer, save :: i_sponge_layer_flag
   real(dp),save :: spinup_function_tide, spinup_function_baroclinic
      
   integer, save, allocatable, dimension(:) ::  &
   &			i_tidal_boundary_temperature_type, 	&	! jw
   &			i_tidal_boundary_salinity_type			! jw

   character(len = 10), save, allocatable, dimension(:) :: 	&
   &			tidal_constituent_name_at_ob,			&	! jw
   &  		tide_name										! jw

   real(dp),save, allocatable, dimension(:) :: 	&
   &			tide_Q, 										& 	! jw
   &			tth, 											& 	! jw
   &			sth, 											& 	! jw
   &			ath		  										! jw
	! jw


	! jw
   real(dp),save,allocatable,dimension(:,:) ::	&	! jw
   &			AinvDeltaZ1, 								&	! jw
   &			AinvDeltaZ2,								&	! jw
   &			AinvG1, 										&	! jw
   &			AinvG2											! jw

   

	
	! jw
	integer, save, allocatable, dimension(:,:)::	&  
	&			num_sub_elm_iteration						! jw
   
   real(dp),save :: denominator_min_for_matrix		! jw

	! jw
   real(dp),save,allocatable,dimension(:) :: 	&
   &			bottom_roughness, 						&	! jw
   &			bottom_drag_coefficient,				&	! jw
   &			Gamma_B,										&	! jw
   &			Gamma_T,										&	! jw
   &			Cdb												! jw
   
   ! jw
   ! jw
   ! jw
   real(dp),save,allocatable,dimension(:,:) :: 	&
   &			salt_cell,									&	! jw
   &			salt_cell_new,								&	! jw
   &			salt_node,									&	! jw
   &			salt_face,									&	! jw
   &			temp_cell,									&	! jw
   &			temp_cell_new,								&	! jw
   &			temp_node,									&	! jw
   &			temp_face										! jw
   
! jw
! jw
! jw
! jw
! jw
! jw
! jw

	
   real(dp),save,allocatable,dimension(:) :: 	&
   &			air_temperature1, 						&	! jw
   &			air_temperature2,							&	! jw
   &			srad, 										&	! jw
   &			hradu, 										&	! jw
   &			hradd, 										&	! jw
   &			shum1, 										&	! jw
   &			shum2, 										&	! jw
   &			sflux, 										&	! jw
   &			fluxsu, 										&	! jw
   &			fluxlu 											! jw
   
      
   ! jw
   ! jw
   real(dp),save,allocatable,dimension(:,:) :: 	&	! jw
   &			rho_face,									&	! jw
   &			rho_cell,									&	! jw
   &			rho_node											! jw
   
   real(dp),save :: salinity_min, salinity_max, temperature_min, temperature_max
   integer, save :: i_density_flag


   character(len=2), save :: i_turbulence_model_name, stability_function 
   integer, save :: i_turbulence_flag
   real(dp),save :: qd, qd2, vd, td
      
   real(dp),save :: schk, schpsi
   real(dp),save :: cmiu0, cpsi1, cpsi2, rpub, rmub, rnub, psimin, eps_min, bgdiff, h1_pp, h2_pp,   &
   &       	tdmin_pp, vdmax_pp2, vdmin_pp2, vdmax_pp1, vdmin_pp1 , q2min 

   real(dp),save,allocatable,dimension(:) :: 	&
   &			sponge_relax, 								& 	! jw
   &			etaic												! jw
   real(dp),save :: ttt, qq 
	
	! jw
   integer, save :: 										&
   &        heat_model_flag, 							&
   &        initial_temperature_field_flag, initial_salinity_field_flag,      &
   &        i_transport_model_flag
   integer, save :: ifile
	! jw
	!==========================================================================!
	


	! jw
! jw
! jw
! jw
! jw
! jw
	real(dp),save,allocatable,dimension(:,:) ::	&
	&			tem0,											&	! jw
	&			sal0,											&	! jw
	&			q2,											&	! jw
	&			xl,											&	! jw
! jw
! jw
! jw
	&			srho,											&	! jw
	&			erho,											&	! jw
	&			prho												! jw
	
end module mod_global_variables

