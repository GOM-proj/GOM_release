!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
module mod_file_definition
	implicit none
	! declare input & output folders of your project ==========================!
 	character(len=100), parameter :: inp_folder = "./input"
 	character(len=100), parameter :: out_folder = "./output"
	! =========================================================================!
	
	!==========================================================================!
	! 11 	~  99: input files
	! 101 ~ 999: output files
	!==========================================================================!
	integer, parameter :: 					&
	&	pw_main_inp			 		= 11, 	&
	&	pw_node_inp 				= 12, 	&
	&	pw_cell_inp					= 13, 	&
	&	pw_hurricane_ser			= 14, 	&		! wind& pressure (fort.22) data input file
	&	pw_hurricane_center		= 15,		&		! hurricane center (eye) position data, fort.23
	&	pw_eta_ser					= 16,		&		! eta_ser1.inp or eta_ser2.inp
! 	&	pw_eta_ser_interp			= 17,		&		! eta_ser_interp.inp
	&	pw_q_ser						= 18,		&		! q_ser.inp
	&	pw_harmonic_ser			= 19,		&		! harmonic_ser.inp
	&	pw_windp_ser				= 20,		&		! windp_ser.inp; wind & pressure time series
! 	&	pw_windp_model				= 21,		&		! windp_model.inp; wind & pressure from wind model
	&	pw_salt_init				= 22,		&		! salt_init.inp; initial salinity at cell center
	&	pw_temp_init				= 23,		&		! temp_init.inp; initial temperature at cell center; not yet included
	&	pw_bottom_roughness_inp	= 24,		&		! bottom_roughness.inp; spatially varying bottom roughness
	&	pw_salt_ser					= 25,		&		! salt_ser.inp; salinity time series
	&	pw_temp_ser					= 26,		&		! temp_ser.inp; temperature time series
	&	pw_Obc_info					= 27,		&		! Obc_info.inp; open boundary element information
	&	pw_Qbc_info					= 28,		&		! Qbc_info.inp; River boundary element information
	&	pw_tser_station_info		= 29,		&		! tser_station_info.inp; time series output station information
	&	pw_air_ser					= 30,		&		! air_ser.inp; atmospheric information time series
	&	pw_restart_inp				= 99,		&		! restart_***.inp
	! From here, output control ===============================================!
	&	pw_node_mirr				= 101,	&
	&	pw_main_mirr				= 102,	&
	&	pw_cell_mirr				= 103, 	&
	&	pw_voronoi_center			= 106,	&		! voronoi center information
	&	pw_voronoi_face_center	= 107,	&		! voronoi face center information
	&	pw_eta_tser_from_msl_tec= 111,	&		! eta time series from mean sea level (msl)
	&	pw_H_tser_tec				= 112,	&		! total water depth time sereis
	&	pw_u_tser_tec				= 113,	&		! u time series
	&	pw_v_tser_tec				= 114,	&		! v time series
	&	pw_salt_tser_tec			= 115,	&		! salt time series
	&	pw_temp_tser_tec			= 116, 	&		! temp time series
	&	pw_airp_tser_tec			= 117,	&		! air pressure time series
	&	pw_eta_tser_from_msl_vtk= 121,	&		! eta time series from mean sea level (msl)
	&	pw_H_tser_vtk				= 122,	&		! total water depth time sereis
	&	pw_u_tser_vtk				= 123,	&		! u time series
	&	pw_v_tser_vtk				= 124,	&		! v time series
	&	pw_salt_tser_vtk			= 125,	&		! salt time series
	&	pw_temp_tser_vtk			= 126, 	&		! temp time series
	&	pw_airp_tser_vtk			= 127,	&		! air pressure time series
! 	&	pw_wind_u_tser				= 114,	&		! wind u time series
! 	&	pw_wind_v_tser				= 115,	&		! wind v time series
! 	&	pw_pressure_tser			= 116,	&		! pressure time series, substitued with airp_tser
! 	&	pw_water_depth_3D			= 117,	&		! 
! 	&	pw_velocity_field			= 118,	&		! vertically averaged velocities
! 	&	pw_all_field_data			= 119,	&
	&	pw_check_Q					= 130,	&		! check_Q.dat
	&	pw_check_grid_2D_tec		= 141, 	&
	&	pw_check_grid_2DO_tec	= 142,	&
	&	pw_check_grid_3D_tec		= 143,	&
	&	pw_check_grid_2D_vtk		= 151, 	&
	&	pw_check_grid_2DO_vtk	= 152,	&
	&	pw_check_grid_3D_vtk		= 153,	&
	&	pw_check_grid_info		= 160,	&
	&	pw_tec2D						= 171, 	&
	&	pw_tec2D_binary			= 172, 	&
	&	pw_tec3D_surf				= 173, 	&
	&	pw_tec3D_full				= 174, 	&
	&	pw_vtk2D						= 180,	&
	&	pw_vtk3D_surf				= 181,	&
	&	pw_vtk3D_full				= 182,	&	
	&	pw_flood_map_tec			= 191, 	&
	&	pw_flood_map_vtk			= 192, 	&
	&	pw_dump2D					= 201,	&
	&	pw_dump3D					= 202,	&
	&	pw_restart_out				= 250,	&
	&	pw_netcdf_out				= 351,	&
	&	pw_netcdf_grid				= 352,	&
	! diagnostic files --------------------------------------------------------!
	&	pw_dia_geometry			= 300,	&		! check geometry in set_geometry.f90
	&	pw_dia_advection			= 301, 	&		! diagnostic file for 'solve_nonlinear_advection.f90'
	&	pw_dia_momentum			= 302,	&		! diagnostic file for 'solve_momentum_equation.f90'
	&	pw_dia_freesurface		= 303,	&		! diagnostic file for 'solve_free_surface_equation.f90'
	&	pw_dia_eta_at_ob			= 304,	&		! diagnostic file in	 'solve_free_surface_equation.f90'
	&	pw_dia_bottom_friction 	= 305,	&		! diagnostic file for 'calculate_bottom_friction.f90'
	&	pw_dia_face_velocity_uv	= 306,	&		! diagnostic file for 'calculate_horizontal_velocities.f90'
	&	pw_dia_face_velocity_w	= 307,	&		! diagnostic file for 'calculate_vertical_velocities.f90'
	&	pw_dia_node_velocity		= 308,	&		! diagnostic file for 'calculate_velocity_at_node.f90'
	! log files ---------------------------------------------------------------!
	&	pw_run_log					= 900				! locate run_log at the end
		
	
	character(len = 100), parameter :: 	&
	&	id_main_inp 				= trim(inp_folder)//'/main.inp',						&
	&	id_node_inp 				= trim(inp_folder)//'/node.inp',						&
	&	id_cell_inp					= trim(inp_folder)//'/cell.inp',						&
	&	id_hurricane_ser			= trim(inp_folder)//'/hurricane_ser.inp',			&	! wind & pressure data input file
	&	id_hurricane_center		= trim(inp_folder)//'/hurricane_center.inp',		&	
	&	id_eta_ser1					= trim(inp_folder)//'/eta_ser1.inp',				&
	&	id_eta_ser2					= trim(inp_folder)//'/eta_ser2.inp',				&
! 	&	id_eta_ser_interp			= trim(inp_folder)//'/eta_ser_interp.inp',		&
	&	id_q_ser1					= trim(inp_folder)//'/q_ser1.inp',					&
	&	id_q_ser2					= trim(inp_folder)//'/q_ser2.inp',					&
	&	id_harmonic_ser			= trim(inp_folder)//'/harmonic_ser.inp',			&
	&	id_windp_ser				= trim(inp_folder)//'/windp_ser.inp',				&
! 	&	id_windp_model				= trim(inp_folder)//'/windp_model.inp',			&
	&	id_salt_init				= trim(inp_folder)//'/salt_init.inp',				&
	&	id_temp_init				= trim(inp_folder)//'/temp_init.inp',				&
	&	id_bottom_roughness_inp	= trim(inp_folder)//'/bottom_roughness.inp',		&
	&	id_salt_ser1				= trim(inp_folder)//'/salt_ser1.inp',				&
	&	id_salt_ser2				= trim(inp_folder)//'/salt_ser2.inp',				&
	&	id_temp_ser1				= trim(inp_folder)//'/temp_ser1.inp',				&
	&	id_temp_ser2				= trim(inp_folder)//'/temp_ser2.inp',				&
	&	id_Obc_info					= trim(inp_folder)//'/Obc_info.inp',				&
	&	id_Qbc_info					= trim(inp_folder)//'/Qbc_info.inp',				&
	&	id_tser_station_info		= trim(inp_folder)//'/tser_station_info.inp',	&
	&	id_air_ser					= trim(inp_folder)//'/air_ser.inp',					&
	&	id_restart_inp				= trim(inp_folder)//'/restart.inp',					&
	! From here, output control ===============================================!
	&	id_node_mirr				= trim(out_folder)//'/node_mirr.dat',				&
	&	id_main_mirr				= trim(out_folder)//'/main_mirr.dat',				&
	&	id_cell_mirr				= trim(out_folder)//'/cell_mirr.dat',				&
	&	id_voronoi_center 		= trim(out_folder)//'/voronoi_center.dat',		&
	&	id_voronoi_face_center	= trim(out_folder)//'/voronoi_face_center.dat',	&
	&	id_eta_tser_from_msl_tec= trim(out_folder)//'/eta_tser_from_msl.dat',	&	! .dat or .txt
	&	id_H_tser_tec				= trim(out_folder)//'/H_tser.dat',					&	! .dat or .txt
	&	id_u_tser_tec				= trim(out_folder)//'/u_tser.dat',					&	! .dat or .txt
	&	id_v_tser_tec				= trim(out_folder)//'/v_tser.dat',					&	! .dat or .txt
	&	id_salt_tser_tec			= trim(out_folder)//'/salt_tser.dat',				&	! .dat or .txt
	&	id_temp_tser_tec			= trim(out_folder)//'/temp_tser.dat',				&	! .dat or .txt
	&	id_airp_tser_tec			= trim(out_folder)//'/airp_tser.dat',				&	! .dat or .txt
	&	id_eta_tser_from_msl_vtk= trim(out_folder)//'/eta_tser_from_msl.txt',	&	! .dat or .txt
	&	id_H_tser_vtk				= trim(out_folder)//'/H_tser.txt',					&	! .dat or .txt
	&	id_u_tser_vtk				= trim(out_folder)//'/u_tser.txt',					&	! .dat or .txt
	&	id_v_tser_vtk				= trim(out_folder)//'/v_tser.txt',					&	! .dat or .txt
	&	id_salt_tser_vtk			= trim(out_folder)//'/salt_tser.txt',				&	! .dat or .txt
	&	id_temp_tser_vtk			= trim(out_folder)//'/temp_tser.txt',				&	! .dat or .txt
	&	id_airp_tser_vtk			= trim(out_folder)//'/airp_tser.txt',				&	! .dat or .txt
! 	&	id_wind_u_tser				= trim(out_folder)//'/wind_u_tser',					&	! .dat or .txt
! 	&	id_wind_v_tser				= trim(out_folder)//'/wind_v_tser',					&	! .dat or .txt
! 	&	id_pressure_tser			= trim(out_folder)//'/pressure_tser',				&	! .dat or .txt, substitued with airp_tser
! 	&	id_water_depth_3D			= trim(out_folder)//'/water_depth_3D.dat',		&
! 	&	id_velocity_field			= trim(out_folder)//'/velocity_field.dat',		&
! 	&	id_all_field_data			= trim(out_folder)//'/all_field_data.dat',		&
	&	id_check_Q					= trim(out_folder)//'/check_Q.dat',					&
	&	id_check_grid_2D_tec		= trim(out_folder)//'/check_grid_2D.dat',			&
	&	id_check_grid_2DO_tec	= trim(out_folder)//'/check_grid_2DO.dat',		&
	&	id_check_grid_3D_tec		= trim(out_folder)//'/check_grid_3D.dat',			&
	&	id_check_grid_2D_vtk		= trim(out_folder)//'/check_grid_2D.vtk',			&
	&	id_check_grid_2DO_vtk	= trim(out_folder)//'/check_grid_2DO.vtk',		&
	&	id_check_grid_3D_vtk		= trim(out_folder)//'/check_grid_3D.vtk',			&
	&	id_check_grid_info		= trim(out_folder)//'/check_grid_info.dat',		&
	&	id_tec2D						= trim(out_folder)//'/tec2D_',						&
	&	id_tec2D_binary			= trim(out_folder)//'/tec2D_binary.dat',			&	
	&	id_tec3D_surf 				= trim(out_folder)//'/tec3D_surf_',					&
	&	id_tec3D_full				= trim(out_folder)//'/tec3D_full_', 				&	
	&	id_vtk2D						= trim(out_folder)//'/vtk2D_',						&
	&	id_vtk3D_surf 				= trim(out_folder)//'/vtk3D_surf_',					&
	&	id_vtk3D_full				= trim(out_folder)//'/vtk3D_full_', 				&	
	&	id_flood_map_tec			= trim(out_folder)//'/tec2D_flood_map.dat',		&
	&	id_flood_map_vtk			= trim(out_folder)//'/vtk2D_flood_map.vtk',		&
	&	id_dump2D					= trim(out_folder)//'/dump2D_',						&
	&	id_dump3D					= trim(out_folder)//'/dump3D_',						&
	&	id_restart_out				= trim(out_folder)//'/restart_',						&
	&	id_netcdf_out				= trim(out_folder)//'/netcdf_out_',					&
	&	id_netcdf_grid				= trim(out_folder)//'/netcdf_grid.nc',				&
	! diagnostic files --------------------------------------------------------!
	&	id_dia_geometry			= trim(out_folder)//'/dia_check_geometry.dat',			&
	&	id_dia_advection			= trim(out_folder)//'/dia_nonlinear_advection.dat',	&
	&	id_dia_momentum			= trim(out_folder)//'/dia_momentum_equation.dat',		&
	&	id_dia_freesurface		= trim(out_folder)//'/dia_freesurface_equation.dat',	&
	&	id_dia_eta_at_ob			= trim(out_folder)//'/dia_eta_at_ob.dat',					&
	&	id_dia_bottom_friction	= trim(out_folder)//'/dia_bottom_friction.dat',			&
	&	id_dia_face_velocity_uv	= trim(out_folder)//'/dia_face_velocity_uv.dat',		&
	&	id_dia_face_velocity_w	= trim(out_folder)//'/dia_face_velocity_w.dat',			&
	&	id_dia_node_velocity		= trim(out_folder)//'/dia_node_velocity.dat',			&
	! log files ---------------------------------------------------------------!
	&	id_run_log 					= trim(out_folder)//'/run.log'							! locate run_log at the end
end module mod_file_definition