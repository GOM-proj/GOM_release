!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Note: Main body data are identical for Tecplot (***.dat) and Paraview (*.txt):
!! 	Actually, Tecplot doesn't require to separate data with commas, but Paraview does.
!!		Thus, I include commas to separate data for general purpose.
!! 
subroutine write_tser_out_files
   use mod_global_variables
   use mod_file_definition
   implicit none
   
   integer :: k,l
   integer :: is, ie, nd
   real(dp),dimension(maxlayer,tser_station_num) :: u_cell, v_cell
   real(dp),dimension(tser_station_num) :: ubar, vbar, sbar, tbar
   integer ,dimension(tser_station_num) :: num_vertical_layer   
   character(len=40) :: format1, format2
   ! End of local variables ==================================================!

	! prepare for vertically averaed values
	ubar = 0.0_dp
	vbar = 0.0_dp
	sbar = 0.0_dp
	tbar = 0.0_dp
	u_cell = 0.0_dp
	v_cell = 0.0_dp
	num_vertical_layer = 1

	write(format1,'(A,I0,A)') '(E15.5E4, ',tser_station_num,'(A,E15.5E4))'
	write(format2,'(A,I0,A)') '(E15.5E4, ',(maxlayer+1)*tser_station_num,'(A,E15.5E4))'
   	
	! write time series files =================================================!
	! main body of time series files is identical between Tecplot and VTK files
	if(tser_eta == 1) then ! water elevation from msl =========================!
	   if(tser_hloc == 1) then
  			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_eta_tser_from_msl_tec,format1) julian_day*tser_time_conv, &
			   &		(',', eta_cell(tser_station_cell(is)), is=1,tser_station_num)
			end if
  			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_eta_tser_from_msl_vtk,format1) julian_day*tser_time_conv, &
			   &		(',', eta_cell(tser_station_cell(is)), is=1,tser_station_num)
			end if			
		else if(tser_hloc == 2) then
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_eta_tser_from_msl_tec,format1) julian_day*tser_time_conv, &
			   &		(',', eta_node(tser_station_node(is)), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_eta_tser_from_msl_vtk,format1) julian_day*tser_time_conv, &
			   &		(',', eta_node(tser_station_node(is)), is=1,tser_station_num)
			end if			
		end if		
	end if
	
	if(tser_H == 1) then ! total water depth ==================================!
		if(tser_hloc == 1) then
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_H_tser_tec,format1) julian_day*tser_time_conv, &
			   &		(',', h_cell(tser_station_cell(is)) + eta_cell(tser_station_cell(is)), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_H_tser_vtk,format1) julian_day*tser_time_conv, &
			   &		(',', h_cell(tser_station_cell(is)) + eta_cell(tser_station_cell(is)), is=1,tser_station_num)
			end if			
		else if(tser_hloc == 2) then
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_H_tser_tec,format1) julian_day*tser_time_conv, &
			   &		(',', h_node(tser_station_node(is)) + eta_node(tser_station_node(is)), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_H_tser_vtk,format1) julian_day*tser_time_conv, &
			   &		(',', h_node(tser_station_node(is)) + eta_node(tser_station_node(is)), is=1,tser_station_num)
			end if			
		end if
	end if
	
	
	! Note:
	! 		u,v,w,salt,temp variables are defined at vertical layers
	! 		Vertically averaged values at node are already calculated in the main program 
	! 		(e.g., ubar_node, vbar_node, sbar_node, tbar_node)
	! 		thus, they don't need to be recalculated here.
 	if(tser_u == 1) then ! u ==================================================!
		if(tser_hloc == 1) then ! at cell
			! calculate velocities at cell center (prism) first
			do is=1,tser_station_num
				ie = tser_station_cell(is)
				
				! note: u_cell is a local variable, but u_node is a global variable
				ubar(is) = 0.0_dp
				u_cell(:,is) = 0.0_dp 

				! If it is a dry element, skip rest part
		      if(top_layer_at_element(ie) == 0) then
		      	cycle ! do nothing and go to the next element
		      end if
				
				num_vertical_layer(is) = top_layer_at_element(ie) - bottom_layer_at_element(ie) + 1
				
				! calcualte value at cell center (prism center)
				do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
					do l=1,tri_or_quad(ie)
						nd = nodenum_at_cell(l,ie)
   					u_cell(k,is) = u_cell(k,is) + u_node(k-1,nd) + u_node(k,nd)
   				end do
   				u_cell(k,is) = u_cell(k,is)/(tri_or_quad(ie)*2) ! prism average
				end do
				
				! calculate vertically averaged value
   			do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
   				ubar(is) = ubar(is) + u_cell(k,is)
   			end do
   			ubar(is) = ubar(is)/num_vertical_layer(is)
			end do
			
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_u_tser_tec,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', ubar(is), (',', u_cell(k,is), k=1,maxlayer), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_u_tser_vtk,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', ubar(is), (',', u_cell(k,is), k=1,maxlayer), is=1,tser_station_num)
			end if			
		else if(tser_hloc == 2) then ! at node
			! ubar_node is already calculated in update_variables.f90
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_u_tser_tec,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', ubar_node(tser_station_node(is)), &
			   &	(',', u_node(k,tser_station_node(is)), k=1,maxlayer), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_u_tser_vtk,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', ubar_node(tser_station_node(is)), &
			   &	(',', u_node(k,tser_station_node(is)), k=1,maxlayer), is=1,tser_station_num)
			end if			
		end if		
	end if
	
	if(tser_v == 1) then	! v ==================================================!
		if(tser_hloc == 1) then ! at cell
			! calculate velocities at cell center (prism) first
			do is=1,tser_station_num
				ie = tser_station_cell(is)

				! note: v_cell is a local variable, but v_node is a global variable
				vbar(is) = 0.0_dp
				v_cell(:,is) = 0.0_dp 
				
				! If it is a dry element, skip rest part
		      if(top_layer_at_element(ie) == 0) then
		      	cycle ! do nothing and go to the next element
		      end if

				num_vertical_layer(is) = top_layer_at_element(ie) - bottom_layer_at_element(ie) + 1
				
				! calcualte value at cell center (prism center)
				do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
					do l=1,tri_or_quad(ie)
						nd = nodenum_at_cell(l,ie)
   					v_cell(k,is) = v_cell(k,is) + v_node(k-1,nd) + v_node(k,nd)
   				end do
   				v_cell(k,is) = v_cell(k,is)/(tri_or_quad(ie)*2) ! prism average
				end do
				
				! calculate vertically averaged value
   			do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
   				vbar(is) = vbar(is) + v_cell(k,is)
   			end do
   			vbar(is) = vbar(is)/num_vertical_layer(is)
			end do
			
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_v_tser_tec,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', vbar(is), (',', v_cell(k,is), k=1,maxlayer), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_v_tser_vtk,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', vbar(is), (',', v_cell(k,is), k=1,maxlayer), is=1,tser_station_num)
			end if			
		else if(tser_hloc == 2) then ! at node
			! vbar_node is already calculated in update_variables.f90
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_v_tser_tec,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', vbar_node(tser_station_node(is)), &
			   &	(',', v_node(k,tser_station_node(is)), k=1,maxlayer), is=1,tser_station_num)				
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_v_tser_vtk,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', vbar_node(tser_station_node(is)), &
			   &	(',', v_node(k,tser_station_node(is)), k=1,maxlayer), is=1,tser_station_num)				
			end if			
		end if
 	end if

	if(tser_salt == 1) then ! salt ============================================!
		if(tser_hloc == 1) then ! at cell
			! calculate salinity at cell center (prism) first
			! The global variable, salt_cell, is already in prism center, so, I don't need to calculate it again here.
			do is=1,tser_station_num
				ie = tser_station_cell(is)

				! note: both salt_cell and salt_node are global variable
				sbar(is) = 0.0_dp

				! If it is a dry element, skip rest part
		      if(top_layer_at_element(ie) == 0) then
		      	cycle ! do nothing and go to the next element
		      end if

				num_vertical_layer(is) = top_layer_at_element(ie) - bottom_layer_at_element(ie) + 1
								
				! calculate vertically averaged value
   			do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
   				sbar(is) = sbar(is) + salt_cell(k,ie)
   			end do
   			sbar(is) = sbar(is)/num_vertical_layer(is)
			end do
			
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_salt_tser_tec,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', sbar(is), (',', salt_cell(k,tser_station_cell(is)), k=1,maxlayer), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_salt_tser_vtk,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', sbar(is), (',', salt_cell(k,tser_station_cell(is)), k=1,maxlayer), is=1,tser_station_num)
			end if			
		else if(tser_hloc == 2) then ! at node
			! sbar_node is already calculated in update_variables.f90
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_salt_tser_tec,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', sbar_node(tser_station_node(is)), &
			   &	(',', salt_node(k,tser_station_node(is)), k=1,maxlayer),is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_salt_tser_vtk,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', sbar_node(tser_station_node(is)), &
			   &	(',', salt_node(k,tser_station_node(is)), k=1,maxlayer),is=1,tser_station_num)
			end if			
		end if
	end if
	
	if(tser_temp == 1) then ! temperature =====================================!
		if(tser_hloc == 1) then ! at cell
			! calculate temperature at cell center (prism) first
			! The global variable, temp_cell, is already in prism center, so, I don't need to calculate it again here.
			do is=1,tser_station_num
				ie = tser_station_cell(is)

				! note: both temp_cell and temp_node are global variable
				tbar(is) = 0.0_dp

				! If it is a dry element, skip rest part
		      if(top_layer_at_element(ie) == 0) then
		      	cycle ! do nothing and go to the next element
		      end if

				num_vertical_layer(is) = top_layer_at_element(ie) - bottom_layer_at_element(ie) + 1
								
				! calculate vertically averaged value
   			do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
   				tbar(is) = tbar(is) + temp_cell(k,ie)
   			end do
   			tbar(is) = tbar(is)/num_vertical_layer(is)
			end do
			
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_temp_tser_tec,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', tbar(is), (',', temp_cell(k,tser_station_cell(is)), k=1,maxlayer), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_temp_tser_vtk,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', tbar(is), (',', temp_cell(k,tser_station_cell(is)), k=1,maxlayer), is=1,tser_station_num)
			end if			
		else if(tser_hloc == 2) then ! at node
			! sbar_node is already calculated in update_variables.f90
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_temp_tser_tec,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', tbar_node(tser_station_node(is)), &
			   &	(',', temp_node(k,tser_station_node(is)), k=1,maxlayer),is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_temp_tser_vtk,format2) &
			   &	julian_day*tser_time_conv, &
			   &	(',', tbar_node(tser_station_node(is)), &
			   &	(',', temp_node(k,tser_station_node(is)), k=1,maxlayer),is=1,tser_station_num)
			end if			
		end if		
	end if

	if(tser_airp == 1) then ! air pressure ====================================!
		if(tser_hloc == 1) then
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_airp_tser_tec,format1) julian_day*tser_time_conv, &
			   &		(',', air_p1(tser_station_cell(is)), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_airp_tser_vtk,format1) julian_day*tser_time_conv, &
			   &		(',', air_p1(tser_station_cell(is)), is=1,tser_station_num)
			end if
		else if(tser_hloc == 2) then
			if(tser_format == 1 .or. tser_format == 3) then
			   write(pw_airp_tser_tec,format1) julian_day*tser_time_conv, &
			   &		(',', airp_at_node(tser_station_node(is)), is=1,tser_station_num)
			end if
			if(tser_format == 2 .or. tser_format == 3) then
			   write(pw_airp_tser_vtk,format1) julian_day*tser_time_conv, &
			   &		(',', airp_at_node(tser_station_node(is)), is=1,tser_station_num)
			end if			
		end if
	end if
end subroutine write_tser_out_files

