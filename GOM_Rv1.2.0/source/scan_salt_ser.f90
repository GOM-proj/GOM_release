!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Find "max_salt_data_num" & "max_salt_ser"
!! 
subroutine scan_salt_ser
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, t, i1, i2
	! jw

	! jw
	if(salt_ser_shape == 1) then
		open(pw_salt_ser, file=id_salt_ser1, form='formatted', status = 'old')

		! jw
		call skip_header_lines(pw_salt_ser,id_salt_ser1)
		
		max_salt_data_num = 1 ! jw
		max_salt_ser = num_salt_ser
		
		do i = 1, num_salt_ser
			read(pw_salt_ser,*) i1, i2
			max_salt_data_num = MAX(max_salt_data_num,i2) ! jw
			do t=1,i2
				read(pw_salt_ser,*)
			end do
		end do
	else if(salt_ser_shape == 2) then
		open(pw_salt_ser, file=id_salt_ser2, form='formatted', status = 'old')

		! jw
		call skip_header_lines(pw_salt_ser,id_salt_ser2)
		
		max_salt_data_num = 1 ! jw
		max_salt_ser = num_salt_ser
		
		read(pw_salt_ser,*) i1, i2
		max_salt_data_num = MAX(max_salt_data_num,i2) ! jw
	end if

	close(pw_salt_ser)
end subroutine scan_salt_ser