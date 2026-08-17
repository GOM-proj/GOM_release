!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! This subroutine is almost identical to write_tec3D_full_head_z.f90
!!  
subroutine write_tec3D_full_body_z_v2
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, j, k, l
	integer :: count1, tot_node, tot_cell, b_layer, t_layer
	real(dp):: sum1, u_cell, v_cell, w_cell
	! jw
	! jw
	character(len=15) :: zonetype = 'FEBRICK'
	character(len= 4) :: time_char
	! jw

	! jw
	! jw

	! jw
	! jw

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
	write(pw_tec3D_full,'(A,A, A,F15.5, A,F15.5, A)') &
	&	', ZONETYPE=',zonetype, &
	&	', SOLUTIONTIME=', julian_day * IS3D_time_conv, &
	&	', VARSHARELIST=([1-2]=1), CONNECTIVITYSHAREZONE=1, SOLUTIONTIME=', julian_day * IS3D_time_conv, & ! jw
	&	', STRANDID=1' 

	! jw
	! jw
	! jw
	! jw
	! jw
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
end subroutine write_tec3D_full_body_z_v2