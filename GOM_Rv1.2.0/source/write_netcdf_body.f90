!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_netcdf_body
	use mod_global_variables
	use mod_file_definition
	use mod_netcdf_utils
	use netcdf
	implicit none
	
	integer :: i
	
	! jw
	integer :: it_array(1), local_time_array(5)
	real(dp):: elapsed_time_hr_array(1), elapsed_time_day_array(1), julian_day_array(1)
	! jw
	
	! jw
	! jw
	it_array(1) = it
	elapsed_time_hr_array(1) = elapsed_time_hr
	elapsed_time_day_array(1) = elapsed_time/86400.0_dp
	julian_day_array(1) = julian_day+1

 	local_time_array(1) = year
 	local_time_array(2) = month
 	local_time_array(3) = day
 	local_time_array(4) = hour
 	local_time_array(5) = minute

	call netcdf_check(nf90_put_var(ncid, varid_it, it_array, start=(/netcdf_buff3/), count=(/1/)))
 	call netcdf_check(nf90_put_var(ncid, varid_elapsed_time_hr, elapsed_time_hr_array, start=(/netcdf_buff3/), count=(/1/)))
 	call netcdf_check(nf90_put_var(ncid, varid_elapsed_time_day, elapsed_time_day_array, start=(/netcdf_buff3/), count=(/1/)))
 	call netcdf_check(nf90_put_var(ncid, varid_julian_day, julian_day_array, start=(/netcdf_buff3/), count=(/1/)))
 	call netcdf_check(nf90_put_var(ncid, varid_local_time, local_time_array, start=(/1,netcdf_buff3/), count=(/5,1/)))
	
	! jw
	if(maxval(netcdf_variable_node) > 0) then
	 	call netcdf_check(nf90_put_var(ncid, varid_blayer_node, bottom_layer_at_node, start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
	 	call netcdf_check(nf90_put_var(ncid, varid_tlayer_node, top_layer_at_node, start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
	end if
	if(maxval(netcdf_variable_cell) > 0) then
	 	call netcdf_check(nf90_put_var(ncid, varid_blayer_cell, bottom_layer_at_element, start=(/1,netcdf_buff3/), count=(/maxele,1/)))
	 	call netcdf_check(nf90_put_var(ncid, varid_tlayer_cell, top_layer_at_element, start=(/1,netcdf_buff3/), count=(/maxele,1/)))
	end if
	if(maxval(netcdf_variable_face) > 0) then
	 	call netcdf_check(nf90_put_var(ncid, varid_blayer_face, bottom_layer_at_face, start=(/1,netcdf_buff3/), count=(/maxface,1/)))
	 	call netcdf_check(nf90_put_var(ncid, varid_tlayer_face, top_layer_at_face, start=(/1,netcdf_buff3/), count=(/maxface,1/)))
	end if
	
	! jw

	! jw
	if(maxval(netcdf_variable_node) > 0) then
		do i=1,15
			if(netcdf_variable_node(i) == 1) then
				if(i==1) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_eta_node, eta_node, &
					&	start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
				else if(i==2) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_u_node, u_node, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer+1,maxnod,1/)))
				else if(i==3) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_v_node, v_node, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer+1,maxnod,1/)))
				else if(i==4) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_w_node, w_node, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer+1,maxnod,1/)))
				else if(i==5) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_salt_node, salt_node, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer,maxnod,1/)))
				else if(i==6) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_temp_node, temp_node, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer,maxnod,1/)))
				else if(i==7) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_rho_node, rho_node, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer,maxnod,1/)))
				else if(i==8) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_wind_u_at_node, wind_u_at_node, &
					&	start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
				else if(i==9) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_wind_v_at_node, wind_v_at_node, &
					&	start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
				else if(i==10) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_airp_at_node, airp_at_node, &
					&	start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
				else if(i==11) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_ubar_node, ubar_node, &
					&	start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
				else if(i==12) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_vbar_node, vbar_node, &
					&	start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
				else if(i==13) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_sbar_node, sbar_node, &
					&	start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
				else if(i==14) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_tbar_node, tbar_node, &
					&	start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
				else if(i==15) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_rbar_node, rbar_node, &
					&	start=(/1,netcdf_buff3/), count=(/maxnod,1/)))
				end if
			end if
		end do
	end if

	! jw
	if(maxval(netcdf_variable_cell) > 0) then
		do i=1,6
			if(netcdf_variable_cell(i) == 1) then
				if(i==1) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_eta_cell, eta_cell_new, &
					&	start=(/1,netcdf_buff3/), count=(/maxele,1/)))
				else if(i==2) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_wn_cell, wn_cell_new, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer+1,maxele,1/)))
				else if(i==3) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_salt_cell, salt_cell_new, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer,maxele,1/)))
				else if(i==4) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_temp_cell, temp_cell_new, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer,maxele,1/)))
				else if(i==5) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_rho_cell, rho_cell, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer,maxele,1/)))
				else if(i==6) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_Kv, Kv, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer+1,maxele,1/)))
				end if
			end if
		end do
	end if

	! jw
	if(maxval(netcdf_variable_face) > 0) then
		do i=1,2
			if(netcdf_variable_face(i) == 1) then
				if(i==1) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_Av, Av, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer+1,maxface,1/)))
				else if(i==2) then ! jw
					call netcdf_check(nf90_put_var(ncid, varid_Kh, Kh, &
					&	start=(/1,1,netcdf_buff3/), count=(/maxlayer,maxface,1/)))
				end if
			end if
		end do
	end if

end subroutine write_netcdf_body