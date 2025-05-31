!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Find "max_temp_data_num" & "max_temp_ser"
!! 
subroutine scan_temp_ser
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, t, i1, i2
	! jw

	! jw
	if(temp_ser_shape == 1) then
		open(pw_temp_ser, file=id_temp_ser1, form='formatted', status = 'old')

		! jw
		call skip_header_lines(pw_temp_ser,id_temp_ser1)
		
		max_temp_data_num = 1 ! jw
		max_temp_ser = num_temp_ser
		
		do i = 1, num_temp_ser
			read(pw_temp_ser,*) i1, i2
			max_temp_data_num = MAX(max_temp_data_num,i2) ! jw
			do t=1,i2
				read(pw_temp_ser,*)
			end do
		end do
	else if(temp_ser_shape == 2) then
		open(pw_temp_ser, file=id_temp_ser2, form='formatted', status = 'old')

		! jw
		call skip_header_lines(pw_temp_ser,id_temp_ser2)
		
		max_temp_data_num = 1 ! jw
		max_temp_ser = num_temp_ser
		
		read(pw_temp_ser,*) i1, i2
		max_temp_data_num = MAX(max_temp_data_num,i2) ! jw
	end if

	close(pw_temp_ser)
end subroutine scan_temp_ser