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
	! integer :: tec_node
	! integer :: tec_node_bottom(4), tec_node_up(4)
	character(len=15) :: zonetype = 'FEBRICK'
	character(len= 4) :: time_char
	! End of local variables ==================================================!

	! I will use "F" format for x,y,z coordinate and "E" format for other variables.
	! This will allow to show extream values even though we have extreme (wrong) values. Let's re-think about this...

	! Write header lines ======================================================!
	! Note: the starting line is different from write_tec3D_full_head_zto_sigma.f90

	! find the maxnod & maxele for 3D writing:
	! note: this finding should be done every timestep because of the wetting & drying
	
	tot_node = 0
	tot_cell = 0
	do i=1,maxele
		b_layer = bottom_layer_at_element(i)
		t_layer = top_layer_at_element(i)
		
		! even for the dry element, I will show the element
		if(t_layer == 0) then
			t_layer = b_layer
		end if

		tot_node = tot_node + ((t_layer+1 - b_layer) + 1) * tri_or_quad(i) ! i.e., total level numer * tri_or_quad
		tot_cell = tot_cell + (t_layer+1 - b_layer) ! i.e., total layer number
		
		if(t_layer /= 0 .and. t_layer /= maxlayer) then
			tot_node = tot_node + ((maxlayer - t_layer) + 1) * tri_or_quad(i) ! i.e., total level numer * tri_or_quad
			tot_cell = tot_cell + (maxlayer - t_layer) ! i.e., total layer number
		end if
	end do

	write(pw_tec3D_full,'(A,F15.5, A,I10, A,I10, A)', advance = 'no') &
	&	'ZONE T = "', julian_day * IS3D_time_conv,	&
	&	'", N =', tot_node, &
	&	', E = ', tot_cell, &
	&	', DATAPACKING=BLOCK'

	! I will add all variables at cell center in the z-grid system output.
	! if there is no cell centered values selected, 'VARLOCATION=' term should not be included.
	if(sum(IS3D_variable) > 0) then ! if at least one of variable is activated (of couse it should be, but anyway to make sure)
		write(pw_tec3D_full,'(A)',advance = 'no') &
		&	', VARLOCATION=(['
		do i=1,6 ! u,v,w,s,t,rho
			if(IS3D_variable(i) == 1) then
				write(pw_tec3D_full,'(I1,A)',advance = 'no') 3+i, ',' ! here 3 is for (x,y,z)
			end if
			! note: the last comma will not be a problem, so ignore the last comma here
		end do
		write(pw_tec3D_full,'(A)',advance = 'no') &
		&	']=CELLCENTERED)'
	end if
	
	! Note: this is different from write_tec3D_full_head_z.f90
	!   	  there is no "VARSHARELIST" in z-grid system because of the wetting & drying process
	!       so, I cannot share x, y coordinate...
	write(pw_tec3D_full,'(A,A, A,F15.5, A,F15.5, A)') &
	&	', ZONETYPE=',zonetype, &
	&	', SOLUTIONTIME=', julian_day * IS3D_time_conv, &
	&	', VARSHARELIST=([1-2]=1), CONNECTIVITYSHAREZONE=1, SOLUTIONTIME=', julian_day * IS3D_time_conv, & ! note: this is different from write_tec3D_full_head_zto_sigma.f90
	&	', STRANDID=1' 

	! Write selected variables ================================================!	
	! skip writing x & y coordinates but z coordinate should be updated =======!
	! write x,y,z coordiantes at node at vertical levels ----------------------!
	! I will write each element's nodes in vertical (from bottom to top) direction
	! I need each value at (layer# + 1) levels
	! write z-coordinate:
	do i=1,maxele
		b_layer = bottom_layer_at_element(i)
		t_layer = top_layer_at_element(i)
		
		
		! write bottom elevation first
		do l=1,tri_or_quad(i)
			write(pw_tec3D_full,'(F0.5,1x)', advance='no') MSL-h_cell(i)
		end do
			
			
		! then write to the surface
		if(t_layer == 0) then
			! write bottom elevation as the top layer also
			do l=1,tri_or_quad(i)
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') MSL-h_cell(i)
			end do
		else
			do k=b_layer,t_layer
				if(z_level(k) <= eta_cell(i)) then
					! at the mid-layers, use fixed vertical level info:
					do l=1,tri_or_quad(i)
						write(pw_tec3D_full,'(F0.5,1x)', advance='no') z_level(k)
					end do
				else
					! at the surface, use nodal elevation:
					do l=1,tri_or_quad(i)
						write(pw_tec3D_full,'(F0.5,1x)', advance='no') MSL + eta_node(nodenum_at_cell(l,i))
					end do
				end if
			end do
		end if
		
		! insert ghost layers to the bottom
		if(t_layer /= 0 .and. t_layer /= maxlayer) then
			do k=t_layer-1,maxlayer
				do l=1,tri_or_quad(i)
					! write to the all bottom
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') MSL-h_cell(i)
					! write(pw_tec3D_full,'(F0.5,1x)', advance='no') MSL + eta_node(nodenum_at_cell(l,i))
				end do
			end do
		end if
		write(pw_tec3D_full,*) ! line change; this is to avoid "Line Too Long" error in tecplot
	end do
	write(pw_tec3D_full,*) ! line change
	
	! write all variables at cell center --------------------------------------!		
	! Note: 
	! 		all variables are alreay updated with the new values, so old_value = new_value in this time
	! write u at cell ---------------------------------------------------------!
	if(IS3D_variable(1) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! dry cell
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp			
			else
				do k=b_layer,t_layer
					sum1 = 0.0_dp
					do l=1,tri_or_quad(i)
						sum1 = sum1 + u_node(k,nodenum_at_cell(l,i)) ! let's just use top-level velocity for simplicity
					end do
					u_cell = sum1/l
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') u_cell 
				end do
			end if
			
			! insert ghost layers to the bottom
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! line change
		end do
		write(pw_tec3D_full,*)	! put a line
	end if

	! write v at cell ---------------------------------------------------------!
	if(IS3D_variable(2) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! dry cell
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp			
			else
				do k=b_layer,t_layer
					sum1 = 0.0_dp
					do l=1,tri_or_quad(i)
						sum1 = sum1 + v_node(k,nodenum_at_cell(l,i)) ! just use top-level velocity
					end do
					v_cell = sum1/l
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') v_cell 
				end do
			end if

			! insert ghost layers to the bottom
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! line change; this is to avoid "Line Too Long" error in tecplot
		end do
		write(pw_tec3D_full,*)	! put a line
	end if
	
	! write w at cell ---------------------------------------------------------!
	if(IS3D_variable(3) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! dry cell
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp			
			else
				do k=b_layer,t_layer
					sum1 = 0.0_dp
					do l=1,tri_or_quad(i)
						sum1 = sum1 + w_node(k,nodenum_at_cell(l,i)) ! just use top-level velocity
					end do
					w_cell = sum1/l
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') w_cell 
				end do
			end if

			! insert ghost layers to the bottom
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! line change; this is to avoid "Line Too Long" error in tecplot
		end do
		write(pw_tec3D_full,*)	! put a line
	end if
	
	! write salinity at cell center -------------------------------------------!
	if(IS3D_variable(4) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! dry cell
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp							
			else
				do k=b_layer,t_layer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') salt_cell(k,i)
				end do
			end if

			! insert ghost layers to the bottom
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! line change; this is to avoid "Line Too Long" error in tecplot
		end do
		write(pw_tec3D_full,*)	! put a line
	end if
	
	! write temperature at cell center ----------------------------------------!
	if(IS3D_variable(5) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! dry cell
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp							
			else
				do k=b_layer,t_layer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') temp_cell(k,i)
				end do
			end if

			! insert ghost layers to the bottom
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! line change; this is to avoid "Line Too Long" error in tecplot
		end do
		write(pw_tec3D_full,*)	! put a line
	end if
	
	! write density (rho_cell) at cell center ---------------------------------!
	if(IS3D_variable(6) == 1) then
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			if(t_layer == 0) then
				! dry cell
				write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
			else				
				do k=b_layer,t_layer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') rho_cell(k,i)
				end do
			end if

			! insert ghost layers to the bottom
			if(t_layer /= 0 .and. t_layer /= maxlayer) then
				do k=t_layer+1,maxlayer
					write(pw_tec3D_full,'(F0.5,1x)', advance='no') 0.0_dp
				end do
			end if
			write(pw_tec3D_full,*) ! line change; this is to avoid "Line Too Long" error in tecplot
		end do
		write(pw_tec3D_full,*)	! put a line
	end if

	! Write cell connectivity =================================================!
	! skip writing cell connectivity ==========================================!

	! Write text information at the bottom of each zone =======================!	
	if(IS3D_time == 1) then	! second
		time_char = ' sec'
	else if(IS3D_time == 2) then ! minute
		time_char = ' min'
	else if(IS3D_time == 3) then ! hour
		time_char = '  hr'
	else if(IS3D_time == 4) then ! day
		time_char = ' day'
	end if
	write(pw_tec3D_full,'(A,F15.5,A,A,I5)') &
	&	'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, T="TIME = ', &
	&	julian_day * IS3D_time_conv, time_char, '", ZN =', zone_num_3D_full

	write(pw_tec3D_full,*)	! put a line
	
	if(IS3D_binary == 1) then
		! not yet activated
	end if
end subroutine write_tec3D_full_body_z_v2