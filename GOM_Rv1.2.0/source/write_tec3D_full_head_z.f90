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
	call write_tec3D_full_body_z

	if(IS3D_binary == 1) then
		! jw
	end if
end subroutine write_tec3D_full_head_z