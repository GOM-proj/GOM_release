!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! This subroutine will write a 3D tecplot format in z-grid system.
!!  
subroutine write_tec3D_full_head_z
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i
	! End of local variables ==================================================!
	
	! I will use "F" format for x,y,z coordinate and "E" format for other variables.
	! This will allow to show extream values even though we have extreme (wrong) values. Let's re-think about this...
	
	! Write header lines ======================================================!
	write(pw_tec3D_full,'(A)') 'Title = "3D contour plot"'

	! note: do not use:
	!   	if(check_grid_unit_conv == 0.1)
	! 		since, real number annnot be equal to other real number
	if(IS3D_unit_conv > 0.5_dp) then ! 1.0
		write(pw_tec3D_full,'(A)',advance = 'no') 'Variables = "X [m]", "Y [m]", "Z [m]"'
	else 										! 0.001
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
	write(pw_tec3D_full,*) ! change the line		
	
	! write the main body of 3D data
	call write_tec3D_full_body_z

	if(IS3D_binary == 1) then
		! not yet activated
	end if
end subroutine write_tec3D_full_head_z