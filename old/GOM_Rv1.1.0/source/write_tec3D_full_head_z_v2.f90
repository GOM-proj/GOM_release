!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! This subroutine will write a 3D tecplot format in z-grid system.
!!  
subroutine write_tec3D_full_head_z_v2
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, j, k, l
	integer :: count1, tot_node, tot_cell, b_layer, t_layer
	real(dp):: sum1, u_cell, v_cell, w_cell
	integer :: tec_node
	integer :: tec_node_bottom(4), tec_node_up(4)
	character(len=15) :: zonetype = 'FEBRICK'
	character(len= 4) :: time_char
	! jw
	
	! jw
	! jw
	
	! jw
	write(pw_tec3D_full,'(A)') 'Title = "3D contour plot"'

	! jw
	! jw
	! jw
	if(IS3D_unit_conv > 0.5_dp) then ! jw
		write(pw_tec3D_full,'(A)',advance = 'no') 'Variables = "X [m]", "Y [m]", "Z [m]"'
	else 										! jw
		write(pw_tec3D_full,'(A)',advance = 'no') 'Variables = "X [km]", "Y [km]", "Z [m]"'
	end if

	do i=1,6
		if(IS3D_variable(i) == 1 .and. i == 1) then
			write(pw_tec3D_full,'(A)',advance = 'no') ', "u [m/s]"'
		else if(IS3D_variable(i) == 1 .and. i == 2) then
			write(pw_tec3D_full,'(A)',advance = 'no') ', "v [m/s]"'
		else if(IS3D_variable(i) == 1 .and. i == 3) then
			write(pw_tec3D_full,'(A)',advance = 'no') ', "w [m/s]"'
		else if(IS3D_variable(i) == 1 .and. i == 4) then
			write(pw_tec3D_full,'(A)',advance = 'no') ', "Salt [psu]"'
		else if(IS3D_variable(i) == 1 .and. i == 5) then
			write(pw_tec3D_full,'(A)',advance = 'no') ', "Temp [C]"'
		else if(IS3D_variable(i) == 1 .and. i == 6) then
			write(pw_tec3D_full,'(A)',advance = 'no') ', "Rho [kg/m3]"'
		end if
	end do
	write(pw_tec3D_full,*) ! jw
	
	! jw
	! jw
	
	tot_node = 0
	tot_cell = 0
	do i=1,maxele
		b_layer = bottom_layer_at_element(i)
		t_layer = top_layer_at_element(i)
		
		! jw
		if(t_layer == 0) then
			t_layer = b_layer
		end if

		tot_node = tot_node + ((t_layer+1 - b_layer) + 1) * tri_or_quad(i) ! jw
		tot_cell = tot_cell + (t_layer+1 - b_layer) ! jw
		
		if(t_layer /= 0 .and. t_layer /= maxlayer) then
			tot_node = tot_node + ((maxlayer - t_layer) + 1) * tri_or_quad(i) ! jw
			tot_cell = tot_cell + (maxlayer - t_layer) ! jw
		end if
	end do

	write(pw_tec3D_full,'(A,F15.5, A,I10, A,I10, A)', advance = 'no') &
	&	'ZONE T = "', julian_day * IS3D_time_conv,	&
	&	'", N =', tot_node, &
	&	', E = ', tot_cell, &
	&	', DATAPACKING=BLOCK'

	! jw
	! jw
	if(sum(IS3D_variable) > 0) then ! jw
		write(pw_tec3D_full,'(A)',advance = 'no') &
		&	', VARLOCATION=(['
		do i=1,6 ! jw
			if(IS3D_variable(i) == 1) then
				write(pw_tec3D_full,'(I1,A)',advance = 'no') 3+i, ',' ! jw
			end if
			! jw
		end do
		write(pw_tec3D_full,'(A)',advance = 'no') &
		&	']=CELLCENTERED)'
	end if
	
	! jw
	! jw
	! jw
	write(pw_tec3D_full,'(A,A, A,F15.5, A)') &
	&	', ZONETYPE=',zonetype, &
	&	', SOLUTIONTIME=', julian_day * IS3D_time_conv, &
	&	', STRANDID=1' 

	! jw
	! jw
	! jw
	! jw
	! jw
	do i=1,maxele
		b_layer = bottom_layer_at_element(i)
		t_layer = top_layer_at_element(i)
		
		if(t_layer == 0) then
			t_layer = b_layer
		end if
		do k=b_layer-1,t_layer
			do l=1,tri_or_quad(i)
				! jw
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') (x_node(nodenum_at_cell(l,i))-xn_min) * IS3D_unit_conv ! jw
			end do
		end do
		! jw
		if(t_layer /= 0 .and. t_layer /= maxlayer) then
			do k=t_layer-1,maxlayer
				do l=1,tri_or_quad(i)
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') (x_node(nodenum_at_cell(l,i))-xn_min) * IS3D_unit_conv ! jw
				end do
			end do
		end if
		write(pw_tec3D_full,*) ! jw
	end do
	write(pw_tec3D_full,*) ! jw
	
	
	! jw
	do i=1,maxele
		b_layer = bottom_layer_at_element(i)
		t_layer = top_layer_at_element(i)
		
		if(t_layer == 0) then
			t_layer = b_layer
		end if

		do k=b_layer-1,t_layer
			do l=1,tri_or_quad(i)
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') (y_node(nodenum_at_cell(l,i))-yn_min) * IS3D_unit_conv
			end do
		end do
		! jw
		if(t_layer /= 0 .and. t_layer /= maxlayer) then
			do k=t_layer-1,maxlayer
				do l=1,tri_or_quad(i)
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') (y_node(nodenum_at_cell(l,i))-yn_min) * IS3D_unit_conv ! jw
				end do
			end do
		end if	
		write(pw_tec3D_full,*) ! jw
	end do
	write(pw_tec3D_full,*) ! jw
		
	! jw
	do i=1,maxele
		b_layer = bottom_layer_at_element(i)
		t_layer = top_layer_at_element(i)
		
		
		! jw
		do l=1,tri_or_quad(i)
			write(pw_tec3D_full,'(F0.5,1x)', advance='no') MSL-h_cell(i)
		end do
			
			
		! jw
		if(t_layer == 0) then
			! jw
			do l=1,tri_or_quad(i)
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') MSL-h_cell(i)
			end do
		else
			do k=b_layer,t_layer
				if(z_level(k) <= eta_cell(i)) then
					! jw
					do l=1,tri_or_quad(i)
						write(pw_tec3D_full,'(F0.5,1x)', advance='no') z_level(k)
					end do
				else
					! jw
					do l=1,tri_or_quad(i)
						write(pw_tec3D_full,'(F0.5,1x)', advance='no') MSL + eta_node(nodenum_at_cell(l,i))
					end do
				end if
			end do
		end if
		
		! jw
		if(t_layer /= 0 .and. t_layer /= maxlayer) then
			do k=t_layer-1,maxlayer
				do l=1,tri_or_quad(i)
					! jw
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') MSL-h_cell(i)
					! jw
				end do
			end do
		end if
		write(pw_tec3D_full,*) ! jw
	end do
	write(pw_tec3D_full,*) ! jw
	
	! jw
	! jw
	! jw
	! jw
	if(IS3D_variable(1) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! jw
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp			
			else
				do k=b_layer,t_layer
					sum1 = 0.0_dp
					do l=1,tri_or_quad(i)
						sum1 = sum1 + u_node(k,nodenum_at_cell(l,i)) ! jw
					end do
					u_cell = sum1/l
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') u_cell 
				end do
			end if
			
			! jw
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! jw
		end do
		write(pw_tec3D_full,*)	! jw
	end if

	! jw
	if(IS3D_variable(2) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! jw
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp			
			else
				do k=b_layer,t_layer
					sum1 = 0.0_dp
					do l=1,tri_or_quad(i)
						sum1 = sum1 + v_node(k,nodenum_at_cell(l,i)) ! jw
					end do
					v_cell = sum1/l
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') v_cell 
				end do
			end if

			! jw
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! jw
		end do
		write(pw_tec3D_full,*)	! jw
	end if
	
	! jw
	if(IS3D_variable(3) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! jw
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp			
			else
				do k=b_layer,t_layer
					sum1 = 0.0_dp
					do l=1,tri_or_quad(i)
						sum1 = sum1 + w_node(k,nodenum_at_cell(l,i)) ! jw
					end do
					w_cell = sum1/l
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') w_cell 
				end do
			end if

			! jw
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! jw
		end do
		write(pw_tec3D_full,*)	! jw
	end if
	
	! jw
	if(IS3D_variable(4) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! jw
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp							
			else
				do k=b_layer,t_layer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') salt_cell(k,i)
				end do
			end if

			! jw
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! jw
		end do
		write(pw_tec3D_full,*)	! jw
	end if
	
	! jw
	if(IS3D_variable(5) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! jw
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp							
			else
				do k=b_layer,t_layer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') temp_cell(k,i)
				end do
			end if

			! jw
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! jw
		end do
		write(pw_tec3D_full,*)	! jw
	end if
	
	! jw
	if(IS3D_variable(6) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! jw
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
			else				
				do k=b_layer,t_layer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') rho_cell(k,i)
				end do
			end if

			! jw
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! jw
		end do
		write(pw_tec3D_full,*)	! jw
	end if

	! jw
	! jw
	tec_node = 0
	do i=1,maxele
		b_layer = bottom_layer_at_element(i)
		t_layer = top_layer_at_element(i)
		
		if(t_layer == 0) then
			t_layer = b_layer
		end if
		
		if(tri_or_quad(i) == 3) then
			! jw
			do k=b_layer,maxlayer
				tec_node_bottom(1) = tec_node + 1
				tec_node_bottom(2) = tec_node + 2
				tec_node_bottom(3) = tec_node + 3
				tec_node_bottom(4) = tec_node + 3
				tec_node_up(1) = tec_node + 4
				tec_node_up(2) = tec_node + 5
				tec_node_up(3) = tec_node + 6
				tec_node_up(4) = tec_node + 6
				
				! jw
				tec_node = tec_node_bottom(4)
				write(pw_tec3D_full,'(8(I0,1x))') &
				&	tec_node_bottom(1), tec_node_bottom(2), tec_node_bottom(3), tec_node_bottom(4), &
				&	tec_node_up(1), tec_node_up(2), tec_node_up(3), tec_node_up(4)
			end do
			tec_node = tec_node_up(4)
		else if(tri_or_quad(i) == 4) then
			! jw
			do k=b_layer,maxlayer
				tec_node_bottom(1) = tec_node + 1
				tec_node_bottom(2) = tec_node + 2
				tec_node_bottom(3) = tec_node + 3
				tec_node_bottom(4) = tec_node + 4
				tec_node_up(1) = tec_node + 5
				tec_node_up(2) = tec_node + 6
				tec_node_up(3) = tec_node + 7
				tec_node_up(4) = tec_node + 8
				
				! jw
				tec_node = tec_node_bottom(4)
				write(pw_tec3D_full,'(8(I0,1x))') &
				&	tec_node_bottom(1), tec_node_bottom(2), tec_node_bottom(3), tec_node_bottom(4), &
				&	tec_node_up(1), tec_node_up(2), tec_node_up(3), tec_node_up(4)				
			end do
			tec_node = tec_node_up(4)
		end if
	end do
	write(pw_tec3D_full,*)	! jw

	! jw
	if(IS3D_time == 1) then	! jw
		time_char = ' sec'
	else if(IS3D_time == 2) then ! jw
		time_char = ' min'
	else if(IS3D_time == 3) then ! jw
		time_char = '  hr'
	else if(IS3D_time == 4) then ! jw
		time_char = ' day'
	end if
	write(pw_tec3D_full,'(A,F15.5,A,A,I5)') &
	&	'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, T="TIME = ', &
	&	julian_day * IS3D_time_conv, time_char, '", ZN =', zone_num_3D_full

	write(pw_tec3D_full,*)	! jw
	
	if(IS3D_binary == 1) then
		! jw
	end if
end subroutine write_tec3D_full_head_z_v2