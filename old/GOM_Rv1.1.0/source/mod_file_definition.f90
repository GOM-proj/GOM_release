!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
module mod_file_definition
	implicit none
	! jw
 	character(len=100), parameter :: inp_folder = "./input"
 	character(len=100), parameter :: out_folder = "./output"
	! jw
	
	!==========================================================================!
	! jw
	! jw
	!==========================================================================!
	integer, parameter :: 					&
	&	pw_main_inp			 		= 11, 	&
	&	pw_node_inp 				= 12, 	&
	&	pw_cell_inp					= 13, 	&
	&	pw_hurricane_ser			= 14, 	&		! jw
	&	pw_hurricane_center		= 15,		&		! jw
	&	pw_eta_ser					= 16,		&		! jw
! jw
	&	pw_q_ser						= 18,		&		! jw
	&	pw_harmonic_ser			= 19,		&		! jw
	&	pw_windp_ser				= 20,		&		! jw
! jw
	&	pw_salt_init				= 22,		&		! jw
	&	pw_temp_init				= 23,		&		! jw
	&	pw_bottom_roughness_inp	= 24,		&		! jw
	&	pw_salt_ser					= 25,		&		! jw
	&	pw_temp_ser					= 26,		&		! jw
	&	pw_Obc_info					= 27,		&		! jw
	&	pw_Qbc_info					= 28,		&		! jw
	&	pw_tser_station_info		= 29,		&		! jw
	&	pw_air_ser					= 30,		&		! jw
	&	pw_restart_inp				= 99,		&		! jw
	! jw
	&	pw_node_mirr				= 101,	&
	&	pw_main_mirr				= 102,	&
	&	pw_cell_mirr				= 103, 	&
	&	pw_voronoi_center			= 106,	&		! jw
	&	pw_voronoi_face_center	= 107,	&		! jw
	&	pw_eta_tser_from_msl_tec= 111,	&		! jw
	&	pw_H_tser_tec				= 112,	&		! jw
	&	pw_u_tser_tec				= 113,	&		! jw
	&	pw_v_tser_tec				= 114,	&		! jw
	&	pw_salt_tser_tec			= 115,	&		! jw
	&	pw_temp_tser_tec			= 116, 	&		! jw
	&	pw_airp_tser_tec			= 117,	&		! jw
	&	pw_eta_tser_from_msl_vtk= 121,	&		! jw
	&	pw_H_tser_vtk				= 122,	&		! jw
	&	pw_u_tser_vtk				= 123,	&		! jw
	&	pw_v_tser_vtk				= 124,	&		! jw
	&	pw_salt_tser_vtk			= 125,	&		! jw
	&	pw_temp_tser_vtk			= 126, 	&		! jw
	&	pw_airp_tser_vtk			= 127,	&		! jw
! jw
! jw
! jw
! jw
! jw
! jw
	&	pw_check_Q					= 130,	&		! jw
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
	! jw
	&	pw_dia_geometry			= 300,	&		! jw
	&	pw_dia_advection			= 301, 	&		! jw
	&	pw_dia_momentum			= 302,	&		! jw
	&	pw_dia_freesurface		= 303,	&		! jw
	&	pw_dia_eta_at_ob			= 304,	&		! jw
	&	pw_dia_bottom_friction 	= 305,	&		! jw
	&	pw_dia_face_velocity_uv	= 306,	&		! jw
	&	pw_dia_face_velocity_w	= 307,	&		! jw
	&	pw_dia_node_velocity		= 308,	&		! jw
	! jw
	&	pw_run_log					= 900				! jw
		
	
	character(len = 100), parameter :: 	&
	&	id_main_inp 				= trim(inp_folder)//'/main.inp',						&
	&	id_node_inp 				= trim(inp_folder)//'/node.inp',						&
	&	id_cell_inp					= trim(inp_folder)//'/cell.inp',						&
	&	id_hurricane_ser			= trim(inp_folder)//'/hurricane_ser.inp',			&	! jw
	&	id_hurricane_center		= trim(inp_folder)//'/hurricane_center.inp',		&	
	&	id_eta_ser1					= trim(inp_folder)//'/eta_ser1.inp',				&
	&	id_eta_ser2					= trim(inp_folder)//'/eta_ser2.inp',				&
! jw
	&	id_q_ser1					= trim(inp_folder)//'/q_ser1.inp',					&
	&	id_q_ser2					= trim(inp_folder)//'/q_ser2.inp',					&
	&	id_harmonic_ser			= trim(inp_folder)//'/harmonic_ser.inp',			&
	&	id_windp_ser				= trim(inp_folder)//'/windp_ser.inp',				&
! jw
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
	! jw
	&	id_node_mirr				= trim(out_folder)//'/node_mirr.dat',				&
	&	id_main_mirr				= trim(out_folder)//'/main_mirr.dat',				&
	&	id_cell_mirr				= trim(out_folder)//'/cell_mirr.dat',				&
	&	id_voronoi_center 		= trim(out_folder)//'/voronoi_center.dat',		&
	&	id_voronoi_face_center	= trim(out_folder)//'/voronoi_face_center.dat',	&
	&	id_eta_tser_from_msl_tec= trim(out_folder)//'/eta_tser_from_msl.dat',	&	! jw
	&	id_H_tser_tec				= trim(out_folder)//'/H_tser.dat',					&	! jw
	&	id_u_tser_tec				= trim(out_folder)//'/u_tser.dat',					&	! jw
	&	id_v_tser_tec				= trim(out_folder)//'/v_tser.dat',					&	! jw
	&	id_salt_tser_tec			= trim(out_folder)//'/salt_tser.dat',				&	! jw
	&	id_temp_tser_tec			= trim(out_folder)//'/temp_tser.dat',				&	! jw
	&	id_airp_tser_tec			= trim(out_folder)//'/airp_tser.dat',				&	! jw
	&	id_eta_tser_from_msl_vtk= trim(out_folder)//'/eta_tser_from_msl.txt',	&	! jw
	&	id_H_tser_vtk				= trim(out_folder)//'/H_tser.txt',					&	! jw
	&	id_u_tser_vtk				= trim(out_folder)//'/u_tser.txt',					&	! jw
	&	id_v_tser_vtk				= trim(out_folder)//'/v_tser.txt',					&	! jw
	&	id_salt_tser_vtk			= trim(out_folder)//'/salt_tser.txt',				&	! jw
	&	id_temp_tser_vtk			= trim(out_folder)//'/temp_tser.txt',				&	! jw
	&	id_airp_tser_vtk			= trim(out_folder)//'/airp_tser.txt',				&	! jw
! jw
! jw
! jw
! jw
! jw
! jw
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
	! jw
	&	id_dia_geometry			= trim(out_folder)//'/dia_check_geometry.dat',			&
	&	id_dia_advection			= trim(out_folder)//'/dia_nonlinear_advection.dat',	&
	&	id_dia_momentum			= trim(out_folder)//'/dia_momentum_equation.dat',		&
	&	id_dia_freesurface		= trim(out_folder)//'/dia_freesurface_equation.dat',	&
	&	id_dia_eta_at_ob			= trim(out_folder)//'/dia_eta_at_ob.dat',					&
	&	id_dia_bottom_friction	= trim(out_folder)//'/dia_bottom_friction.dat',			&
	&	id_dia_face_velocity_uv	= trim(out_folder)//'/dia_face_velocity_uv.dat',		&
	&	id_dia_face_velocity_w	= trim(out_folder)//'/dia_face_velocity_w.dat',			&
	&	id_dia_node_velocity		= trim(out_folder)//'/dia_node_velocity.dat',			&
	! jw
	&	id_run_log 					= trim(out_folder)//'/run.log'							! jw
end module mod_file_definition