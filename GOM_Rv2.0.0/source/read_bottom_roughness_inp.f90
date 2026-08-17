!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Read spatially varying bottom friction at each node
!! 
subroutine read_bottom_roughness_inp
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, j
	integer :: itemp, serial_num, n1, n2
	real(dp),dimension(maxnod) :: rtemp1, rtemp2
	! jw
	
	rtemp1 = 0.0_dp
	rtemp2 = 0.0_dp
	
	open(pw_bottom_roughness_inp,file=id_bottom_roughness_inp,form='formatted',status='old')

	! jw
	call skip_header_lines(pw_bottom_roughness_inp, id_bottom_roughness_inp)
	
	! jw
	read(pw_bottom_roughness_inp,*) itemp ! jw
	if(itemp /= maxnod) then
		write(*,*) 'Error: maximum node does not match in bottom_roughness.inp'
		write(pw_run_log,*) 'Error: maximum node does not match in bottom_roughness.inp'
		stop
	end if
	do i=1,maxnod
		read(pw_bottom_roughness_inp,*) serial_num, rtemp1(i), rtemp2(i)
	end do	
	close(pw_bottom_roughness_inp)	
	! jw
	
	! jw
	do j=1,maxface
		n1 = nodenum_at_face(1,j)
		n2 = nodenum_at_face(2,j)
		
		! jw
		manning(j) = (rtemp1(n1) + rtemp1(n2)) * 0.5_dp
		bf_height(j) = (rtemp2(n1) + rtemp2(n2)) * 0.5_dp
	end do
end subroutine read_bottom_roughness_inp