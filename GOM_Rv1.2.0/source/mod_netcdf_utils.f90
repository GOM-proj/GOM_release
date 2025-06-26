!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================!
!!
!! Define netcdf-related variables
!!
module mod_netcdf_utils
	use netcdf
	implicit none
	
	integer :: ncid
	
	! jw
	integer :: node_dimid, cell_dimid, face_dimid, layer_dimid, level_dimid ! jw
	integer :: four_dimid, five_dimid
	
	! jw
	integer :: time_dimid 
	
	! jw
	! jw
	integer :: varid_eta_node, varid_u_node, varid_v_node, varid_w_node, varid_salt_node, varid_temp_node, varid_rho_node
	integer :: varid_wind_u_at_node, varid_wind_v_at_node, varid_airp_at_node
	integer :: varid_ubar_node, varid_vbar_node, varid_sbar_node, varid_tbar_node, varid_rbar_node
 	integer :: varid_blayer_node, varid_tlayer_node

	! jw
	integer :: varid_eta_cell, varid_wn_cell, varid_salt_cell, varid_temp_cell, varid_rho_cell, varid_Kv
 	integer :: varid_blayer_cell, varid_tlayer_cell

 	
 	! jw
 	integer :: varid_Av, varid_Kh
 	integer :: varid_blayer_face, varid_tlayer_face
 	
 	! jw
 	integer :: varid_it, varid_elapsed_time_hr, varid_elapsed_time_day, varid_julian_day, varid_local_time
	! jw
end module mod_netcdf_utils
