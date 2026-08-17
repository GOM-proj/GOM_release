!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Tecplot & VTK files require different file formats at the header lines.
!! Thus, this subroutine will write different header lines for each file format.
!! 
subroutine open_tser_outfiles_sub(pw_tser,id_tser,tser_list)
	use mod_global_variables
	use mod_file_definition
	implicit none
	
   integer, intent(in) :: pw_tser, tser_list
   character(len=200), intent(in) :: id_tser	

	integer :: i,k,l
	integer :: is,ie,nd
   character(len=200),dimension(0:maxlayer,tser_station_num) :: variable_name
   real(dp),dimension(tser_station_num,maxlayer) :: u_cell, v_cell
   real(dp),dimension(tser_station_num) :: ubar, vbar, sbar, tbar
   integer,dimension(tser_station_num) ::  num_vertical_layer
	character(len=40) :: format_string
	character(len= 2) :: layer_num_buff   
	character(len=40) :: format1, format2, format3
	! jw
	
	! jw
	ubar = 0.0_dp
	vbar = 0.0_dp
	sbar = 0.0_dp
	tbar = 0.0_dp
	u_cell = 0.0_dp
	v_cell = 0.0_dp
	num_vertical_layer = 1

	write(format1,'(A,I0,A)') '(E15.5E4, ',tser_station_num,'(A,E15.5E4))' 					! jw
	write(format2,'(A,I0,A)') '(E15.5E4, ',(maxlayer+1)*tser_station_num,'(A,E15.5E4))' ! jw
	write(format3,'(A,I0,A)') '(',(maxlayer+1)*tser_station_num,'A)' 							! jw
	
	! jw
	if(tser_list <= 2) then ! jw
		do i=1,tser_station_num
			write(pw_tser,'(A3,A,A1)', advance = 'no') ', "', trim(tser_station_name(i)), '"'
	 	end do
		write(pw_tser,*) ! jw
   else ! jw
	   format_string = '(I2.2)'	! jw
	   do i=1,tser_station_num
	   	do k=0,maxlayer
	   		! jw
				write(layer_num_buff,format_string) k
	   		variable_name(k,i) = ', "'//trim(tser_station_name(i))//'_('//layer_num_buff//')"' ! jw
	   	end do
	   end do
	   ! jw
	   write(pw_tser,format3) ((trim(variable_name(k,i)), k=0,maxlayer), i=1,tser_station_num)
  	end if
	! jw
	
   ! jw
   ! jw
   if(tser_list == 1) then ! jw
   	if(tser_hloc == 1) then ! jw
		   write(pw_tser,format1) &
		   &		julian_day*tser_time_conv, (',', eta_cell(tser_station_cell(is)), is=1,tser_station_num)
		else if(tser_hloc == 2) then ! jw
		   write(pw_tser,format1) &
		   &		julian_day*tser_time_conv, (',', eta_node(tser_station_node(is)), is=1,tser_station_num)
		end if
	else if(tser_list == 2) then ! jw
   	if(tser_hloc == 1) then ! jw
		   write(pw_tser,format1) &
		   &		julian_day*tser_time_conv, &
		   &		(',', h_cell(tser_station_cell(is)) + eta_cell(tser_station_cell(is)), is=1,tser_station_num)
		else if(tser_hloc == 2) then ! jw
		   write(pw_tser,format1) &
		   &		julian_day*tser_time_conv, &
		   &		(',', h_node(tser_station_cell(is)) + eta_node(tser_station_node(is)), is=1,tser_station_num)
		end if
	else if(tser_list == 3) then ! jw
		if(tser_hloc == 1) then ! jw
			! jw
			do is=1,tser_station_num
				ie = tser_station_cell(is)				

				! jw
				ubar(is) = 0.0_dp
				u_cell(is,:) = 0.0_dp 
				
				! jw
		      if(top_layer_at_element(ie) == 0) then
		      	cycle ! jw
		      end if
		      
  				num_vertical_layer(is) = top_layer_at_element(ie) - bottom_layer_at_element(ie) + 1

				! jw
				do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
					do l=1,tri_or_quad(ie)
						nd = nodenum_at_cell(l,ie)
   					u_cell(is,k) = u_cell(is,k) + u_node(k-1,nd) + u_node(k,nd)
   				end do
   				u_cell(is,k) = u_cell(is,k)/(tri_or_quad(ie)*2) ! jw
				end do
				
				! jw
   			do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
   				ubar(is) = ubar(is) + u_cell(is,k)
   			end do
   			ubar(is) = ubar(is)/num_vertical_layer(is)
			end do

		   write(pw_tser,format2) &
		   &	julian_day*tser_time_conv, &
		   &	(',', ubar(is), (',', u_cell(is,k), k=1,maxlayer), is=1,tser_station_num)
   		
		else if(tser_hloc == 2) then ! jw
			! jw
			! jw
			do is=1,tser_station_num
				nd = tser_station_node(is)

				! jw
				ubar(is) = 0.0_dp
					
				if(top_layer_at_node(nd) == 0) then
					cycle ! jw
				end if
					
				num_vertical_layer(is) = top_layer_at_node(nd) - bottom_layer_at_node(nd) + 1
					
				! jw
   			do k=bottom_layer_at_node(nd)-1,top_layer_at_node(nd)
   				ubar(is) = ubar(is) + u_node(k,nd)
   			end do
   			ubar(is) = ubar(is)/(num_vertical_layer(is)+1) ! jw
			end do
				
		   write(pw_tser,format2) &
		   &	julian_day*tser_time_conv, &
		   &	(',', ubar(is), (',', u_node(k,tser_station_node(is)), k=1,maxlayer), is=1,tser_station_num)
		end if
	else if(tser_list == 4) then ! jw
		if(tser_hloc == 1) then ! jw
			! jw
			do is=1,tser_station_num
				ie = tser_station_cell(is)
				
				! jw
				vbar(is) = 0.0_dp
				v_cell(is,:) = 0.0_dp 

				! jw
		      if(top_layer_at_element(ie) == 0) then
		      	cycle ! jw
		      end if

				num_vertical_layer(is) = top_layer_at_element(ie) - bottom_layer_at_element(ie) + 1
				
				! jw
				do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
					do l=1,tri_or_quad(ie)
						nd = nodenum_at_cell(l,ie)
   					v_cell(is,k) = v_cell(is,k) + v_node(k-1,nd) + v_node(k,nd)
   				end do
   				v_cell(is,k) = v_cell(is,k)/(tri_or_quad(ie)*2) ! jw
				end do
				
				! jw
   			do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
   				vbar(is) = vbar(is) + v_cell(is,k)
   			end do
   			vbar(is) = vbar(is)/num_vertical_layer(is)
			end do

		   write(pw_tser,format2) &
		   &	julian_day*tser_time_conv, &
		   &	(',', vbar(is), (',', v_cell(is,k), k=1,maxlayer), is=1,tser_station_num)		   
		else if(tser_hloc == 2) then ! jw
			do is=1,tser_station_num
				nd = tser_station_node(is)

				! jw
				vbar(is) = 0.0_dp

				if(top_layer_at_node(nd) == 0) then
					cycle ! jw
				end if
					
				num_vertical_layer(is) = top_layer_at_node(nd) - bottom_layer_at_node(nd) + 1
					
				! jw
   			do k=bottom_layer_at_node(nd)-1,top_layer_at_node(nd)
   				vbar(is) = vbar(is) + v_node(k,nd)
   			end do
   			vbar(is) = vbar(is)/(num_vertical_layer(is)+1) ! jw
			end do
				
		   write(pw_tser,format2) &
		   &	julian_day*tser_time_conv, &
		   &	(',', vbar(is), (',', v_node(k,tser_station_node(is)), k=1,maxlayer), is=1,tser_station_num)
		end if
	else if(tser_list == 5) then ! jw
		if(tser_hloc == 1) then ! jw
			! jw
			! jw
			do is=1,tser_station_num
				ie = tser_station_cell(is)
				
				sbar(is) = 0.0_dp

				! jw
		      if(top_layer_at_element(ie) == 0) then
		      	cycle ! jw
		      end if

				num_vertical_layer(is) = top_layer_at_element(ie) - bottom_layer_at_element(ie) + 1
								
				! jw
   			do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
   				sbar(is) = sbar(is) + salt_cell(k,ie)
   			end do
   			sbar(is) = sbar(is)/num_vertical_layer(is)
			end do

		   write(pw_tser,format2) &
		   &	julian_day*tser_time_conv, &
		   &	(',', sbar(is), (',', salt_cell(k,tser_station_cell(is)), k=1,maxlayer), is=1,tser_station_num)
		else if(tser_hloc == 2) then ! jw
			do is=1,tser_station_num
				nd = tser_station_node(is)

				sbar(is) = 0.0_dp
				
				if(top_layer_at_node(nd) == 0) then
					cycle ! jw
				end if

				num_vertical_layer(is) = top_layer_at_node(nd) - bottom_layer_at_node(nd) + 1
					
				! jw
   			do k=bottom_layer_at_node(nd)-1,top_layer_at_node(nd)
   				sbar(is) = sbar(is) + salt_node(k,nd)
   			end do
   			sbar(is) = sbar(is)/num_vertical_layer(is) ! jw
			end do
				
		   write(pw_tser,format2) &
		   &	julian_day*tser_time_conv, &
		   &	(',', sbar(is), (',', salt_node(k,tser_station_node(is)), k=1,maxlayer), is=1,tser_station_num)
		end if
	else if(tser_list == 6) then ! jw
		if(tser_hloc == 1) then ! jw
			! jw
			! jw
			do is=1,tser_station_num
				ie = tser_station_cell(is)
				
				tbar(is) = 0.0_dp

				! jw
		      if(top_layer_at_element(ie) == 0) then
		      	cycle ! jw
		      end if

				num_vertical_layer(is) = top_layer_at_element(ie) - bottom_layer_at_element(ie) + 1
								
				! jw
   			do k=bottom_layer_at_element(ie),top_layer_at_element(ie)
   				tbar(is) = tbar(is) + temp_cell(k,ie)
   			end do
   			tbar(is) = tbar(is)/num_vertical_layer(is)
			end do

		   write(pw_tser,format2) &
		   &	julian_day*tser_time_conv, &
		   &	(',', tbar(is), (',', temp_cell(k,tser_station_cell(is)), k=1,maxlayer), is=1,tser_station_num)
		else if(tser_hloc == 2) then ! jw
			do is=1,tser_station_num
				nd = tser_station_node(is)

				tbar(is) = 0.0_dp
				
				if(top_layer_at_node(nd) == 0) then
					cycle ! jw
				end if

				num_vertical_layer(is) = top_layer_at_node(nd) - bottom_layer_at_node(nd) + 1
					
				! jw
   			do k=bottom_layer_at_node(nd)-1,top_layer_at_node(nd)
   				tbar(is) = tbar(is) + temp_node(k,nd)
   			end do
   			tbar(is) = tbar(is)/num_vertical_layer(is) ! jw
			end do
				
		   write(pw_tser,format2) &
		   &	julian_day*tser_time_conv, &
		   &	(',', tbar(is), (',', temp_node(k,tser_station_node(is)), k=1,maxlayer), is=1,tser_station_num)
		end if
	else if(tser_list == 7) then ! jw
	end if
	
	! jw
end subroutine open_tser_outfiles_sub