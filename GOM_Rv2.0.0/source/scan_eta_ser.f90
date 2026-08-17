!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Find "max_eta_data_num" and "max_eta_ser"
!! 
subroutine scan_eta_ser
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, t, i1, i2
	! jw

	! jw
	if(eta_ser_shape == 1) then
		open(pw_eta_ser, file=id_eta_ser1, form='formatted', status = 'old')

		! jw
		call skip_header_lines(pw_eta_ser,id_eta_ser1)
		
		max_eta_data_num = 1 ! jw
		max_eta_ser = num_eta_ser
		
		do i = 1, num_eta_ser
			read(pw_eta_ser,*) i1, i2
			max_eta_data_num = MAX(max_eta_data_num,i2) ! jw
			do t=1,i2
				read(pw_eta_ser,*)
			end do
		end do	
	else if(eta_ser_shape == 2) then
		open(pw_eta_ser, file=id_eta_ser2, form='formatted', status = 'old')

		! jw
		call skip_header_lines(pw_eta_ser,id_eta_ser2)
		
		max_eta_data_num = 1 ! jw
		max_eta_ser = num_eta_ser
		
		read(pw_eta_ser,*) i1, i2
		max_eta_data_num = MAX(max_eta_data_num,i2) ! jw
	end if
	
	close(pw_eta_ser)		
end subroutine scan_eta_ser