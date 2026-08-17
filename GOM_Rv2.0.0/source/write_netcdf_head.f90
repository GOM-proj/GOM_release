!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_netcdf_head
	use mod_global_variables
	use mod_file_definition
	use mod_netcdf_utils
	use netcdf
	implicit none
	
	integer :: i

	! jw
		
	! jw
	! jw
	call netcdf_check(nf90_def_dim(ncid, "maxtime", netcdf_maxtime, time_dimid))
	call netcdf_check(nf90_def_dim(ncid, "maxnod", maxnod, node_dimid))
	call netcdf_check(nf90_def_dim(ncid, "maxele", maxele, cell_dimid))
	call netcdf_check(nf90_def_dim(ncid, "maxface", maxface, face_dimid))
	call netcdf_check(nf90_def_dim(ncid, "maxlayer", maxlayer, layer_dimid))
	call netcdf_check(nf90_def_dim(ncid, "maxlevel", maxlayer+1, level_dimid))
	call netcdf_check(nf90_def_dim(ncid, "five", 5, five_dimid)) ! jw
	
	! jw
	! jw
	call netcdf_check(nf90_def_var(ncid,"it",nf90_int,(/time_dimid/),varid_it))
	call netcdf_check(nf90_put_att(ncid,varid_it,"long_name","model iteration (elapsed) time step"))
	
 	call netcdf_check(nf90_def_var(ncid,"elapsed_time_hr",nf90_double,(/time_dimid/),varid_elapsed_time_hr))
 	call netcdf_check(nf90_put_att(ncid,varid_elapsed_time_hr,"long_name","model elapsed time in [hour]"))
		
 	call netcdf_check(nf90_def_var(ncid,"elapsed_time_day",nf90_double,(/time_dimid/),varid_elapsed_time_day))
 	call netcdf_check(nf90_put_att(ncid,varid_elapsed_time_day,"long_name","model elapsed time in [day]"))

 	call netcdf_check(nf90_def_var(ncid,"julian_day",nf90_double,(/time_dimid/),varid_julian_day))
 	call netcdf_check(nf90_put_att(ncid,varid_julian_day,&
 	& 	"long_name","model elapsed julian day from the given start year; see C3 in main.inp"))
	
 	call netcdf_check(nf90_def_var(ncid,"local_time",nf90_int,(/five_dimid,time_dimid/),varid_local_time))
 	call netcdf_check(nf90_put_att(ncid,varid_local_time,&
 	& 	"long_name","model simulation local time in [yyyy,mm,dd,hh,mm]"))

	if(maxval(netcdf_variable_node) > 0) then
		! jw
	 	call netcdf_check(nf90_def_var(ncid,"blayer_at_node",nf90_int,(/node_dimid,time_dimid/),varid_blayer_node))
	 	call netcdf_check(nf90_put_att(ncid,varid_blayer_node,"long_name","bottom layer number at node"))
		! jw
	 	call netcdf_check(nf90_def_var(ncid,"tlayer_at_node",nf90_int,(/node_dimid,time_dimid/),varid_tlayer_node))
	 	call netcdf_check(nf90_put_att(ncid,varid_tlayer_node,"long_name","top layer number at node"))
	end if
		
	if(maxval(netcdf_variable_cell) > 0) then
		! jw
	 	call netcdf_check(nf90_def_var(ncid,"blayer_at_cell",nf90_int,(/cell_dimid,time_dimid/),varid_blayer_cell))
	 	call netcdf_check(nf90_put_att(ncid,varid_blayer_cell,"long_name","bottom layer number at cell"))
		! jw
	 	call netcdf_check(nf90_def_var(ncid,"tlayer_at_cell",nf90_int,(/cell_dimid,time_dimid/),varid_tlayer_cell))
	 	call netcdf_check(nf90_put_att(ncid,varid_tlayer_cell,"long_name","top layer number at cell"))
	end if
	
	if(maxval(netcdf_variable_face) > 0) then
		! jw
	 	call netcdf_check(nf90_def_var(ncid,"blayer_at_face",nf90_int,(/face_dimid,time_dimid/),varid_blayer_face))
	 	call netcdf_check(nf90_put_att(ncid,varid_blayer_face,"long_name","bottom layer number at face"))
		! jw
	 	call netcdf_check(nf90_def_var(ncid,"tlayer_at_face",nf90_int,(/face_dimid,time_dimid/),varid_tlayer_face))
	 	call netcdf_check(nf90_put_att(ncid,varid_tlayer_face,"long_name","top layer number at face"))
	end if
	
		
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
	! jw
	! jw
	! jw
	! jw

	! jw
	! jw
	! jw

	! jw
	if(maxval(netcdf_variable_node) > 0) then
		do i=1,15
			if(netcdf_variable_node(i) == 1) then
				if(i==1) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/node_dimid,time_dimid/),varid_eta_node)) ! jw
					call netcdf_check(nf90_put_att(ncid,varid_eta_node,"units","[m]"))
					call netcdf_check(nf90_put_att(ncid,varid_eta_node,"long_name",&
					&	"water surface elevation at node [m]"))
				else if(i==2) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/level_dimid,node_dimid,time_dimid/),varid_u_node))
					call netcdf_check(nf90_put_att(ncid,varid_u_node,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_u_node,"long_name",&
					&	"velocity u (true east/west) at node (at each vertical level) [m/s]"))				
				else if(i==3) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/level_dimid,node_dimid,time_dimid/),varid_v_node))
					call netcdf_check(nf90_put_att(ncid,varid_v_node,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_v_node,"long_name",&
					&	"velocity v (true north/south) at node (at each vertical level) [m/s]"))				
				else if(i==4) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/level_dimid,node_dimid,time_dimid/),varid_w_node))
					call netcdf_check(nf90_put_att(ncid,varid_w_node,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_w_node,"long_name",&
					&	"velocity w at node (at each vertical level) [m/s]"))
				else if(i==5) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/layer_dimid,node_dimid,time_dimid/),varid_salt_node))
					call netcdf_check(nf90_put_att(ncid,varid_salt_node,"units","[psu]"))
					call netcdf_check(nf90_put_att(ncid,varid_salt_node,"long_name",&
					&	"salinity at horizontal nodal point at each veritcal layer [psu]"))				
				else if(i==6) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/layer_dimid,node_dimid,time_dimid/),varid_temp_node))
					call netcdf_check(nf90_put_att(ncid,varid_temp_node,"units","[Celcius]"))
					call netcdf_check(nf90_put_att(ncid,varid_temp_node,"long_name","temperature at node [C]"))				
				else if(i==7) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/layer_dimid,node_dimid,time_dimid/),varid_rho_node))
					call netcdf_check(nf90_put_att(ncid,varid_rho_node,"units","[kg/m3]"))
					call netcdf_check(nf90_put_att(ncid,varid_rho_node,"long_name","rho at node [kg/m3]"))				
				else if(i==8) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/node_dimid,time_dimid/),varid_wind_u_at_node))
					call netcdf_check(nf90_put_att(ncid,varid_wind_u_at_node,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_wind_u_at_node,"long_name",&
					&	"wind speed in x-direction at each node [m/s]"))				
				else if(i==9) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/node_dimid,time_dimid/),varid_wind_v_at_node))
					call netcdf_check(nf90_put_att(ncid,varid_wind_v_at_node,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_wind_v_at_node,"long_name",&
					&	"wind speed in y-direction at each node [m/s]"))				
				else if(i==10) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/node_dimid,time_dimid/),varid_airp_at_node))
					call netcdf_check(nf90_put_att(ncid,varid_airp_at_node,"units","[Millibars]"))
					call netcdf_check(nf90_put_att(ncid,varid_airp_at_node,"long_name","air pressure at node [Millibars]"))				
				else if(i==11) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/node_dimid,time_dimid/),varid_ubar_node))
					call netcdf_check(nf90_put_att(ncid,varid_ubar_node,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_ubar_node,"long_name",&
					&	"vertically averaged velocity u (true east/wets) at node [m/s]"))				
				else if(i==12) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/node_dimid,time_dimid/),varid_vbar_node))
					call netcdf_check(nf90_put_att(ncid,varid_vbar_node,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_vbar_node,"long_name",&
					&	"vertically averaged velocity v (true north/south) at node [m/s]"))				
				else if(i==13) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/node_dimid,time_dimid/),varid_sbar_node))
					call netcdf_check(nf90_put_att(ncid,varid_sbar_node,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_sbar_node,"long_name",&
					&	"vertically averaged salinity at node [psu]"))				
				else if(i==14) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/node_dimid,time_dimid/),varid_tbar_node))
					call netcdf_check(nf90_put_att(ncid,varid_tbar_node,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_tbar_node,"long_name",&
					&	"vertically averaged temperature at node [C]"))				
				else if(i==15) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_node(i),nf90_float,&
					&	(/node_dimid,time_dimid/),varid_rbar_node))
					call netcdf_check(nf90_put_att(ncid,varid_rbar_node,"units","[kg/m3]"))
					call netcdf_check(nf90_put_att(ncid,varid_rbar_node,"long_name",&
					&	"vertically averaged water density (rho) at node [kg/m3]"))				
				end if
			end if
		end do
	end if

	! jw
	if(maxval(netcdf_variable_cell) > 0) then
		do i=1,6
			if(netcdf_variable_cell(i) == 1) then
				if(i==1) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_cell(i),nf90_float,&
					&	(/cell_dimid,time_dimid/),varid_eta_cell)) ! jw
					call netcdf_check(nf90_put_att(ncid,varid_eta_cell,"units","[m]"))
					call netcdf_check(nf90_put_att(ncid,varid_eta_cell,"long_name",&
					&	"water surface elevation at cell center [m], [cell]"))
				else if(i==2) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_cell(i),nf90_float,&
					&	(/level_dimid,cell_dimid,time_dimid/),varid_wn_cell))
					call netcdf_check(nf90_put_att(ncid,varid_wn_cell,"units","[m/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_wn_cell,"long_name","velocity w at prism face [m/s], [cell,level]"))				
				else if(i==3) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_cell(i),nf90_float,&
					&	(/layer_dimid,cell_dimid,time_dimid/),varid_salt_cell))
					call netcdf_check(nf90_put_att(ncid,varid_salt_cell,"units","[psu]"))
					call netcdf_check(nf90_put_att(ncid,varid_salt_cell,"long_name","salinity at cell center [psu], [cell,layer]"))				
				else if(i==4) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_cell(i),nf90_float,&
					&	(/layer_dimid,cell_dimid,time_dimid/),varid_temp_cell))
					call netcdf_check(nf90_put_att(ncid,varid_temp_cell,"units","[Celcius]"))
					call netcdf_check(nf90_put_att(ncid,varid_temp_cell,"long_name","temperature at cell center [C], [cell,layer]"))
				else if(i==5) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_cell(i),nf90_float,&
					&	(/layer_dimid,cell_dimid,time_dimid/),varid_rho_cell))
					call netcdf_check(nf90_put_att(ncid,varid_rho_cell,"units","[Kg/m3]"))
					call netcdf_check(nf90_put_att(ncid,varid_rho_cell,"long_name","Density at element center [kg/m3], [cell,layer]"))				
				else if(i==6) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_cell(i),nf90_float,&
					&	(/level_dimid,cell_dimid,time_dimid/),varid_Kv))
					call netcdf_check(nf90_put_att(ncid,varid_Kv,"units","[m2/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_Kv,"long_name","vertical eddy diffusivity at cell [m2/s], [cell,level]"))				
				end if
			end if
		end do
	end if

	! jw
	if(maxval(netcdf_variable_face) > 0) then
		do i=1,2
			if(netcdf_variable_face(i) == 1) then
				if(i==1) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_face(i),nf90_float,&
					&	(/level_dimid,face_dimid,time_dimid/),varid_Av))
					call netcdf_check(nf90_put_att(ncid,varid_Av,"units","[m2/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_Av,"long_name","vertical eddy viscosity [m2/s]"))				
				else if(i==2) then
					! jw
					call netcdf_check(nf90_def_var(ncid,netcdf_variable_name_face(i),nf90_float,&
					&	(/layer_dimid,face_dimid,time_dimid/),varid_Kh))
					call netcdf_check(nf90_put_att(ncid,varid_Kh,"units","[m2/s]"))
					call netcdf_check(nf90_put_att(ncid,varid_Kh,"long_name","Horizontal diffusivity [m2/s]"))				
				end if
			end if
		end do
	end if
	! jw
	
	call netcdf_check(nf90_enddef(ncid))	
end subroutine write_netcdf_head