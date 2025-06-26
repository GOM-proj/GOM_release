!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Find "max_air_data_num"
!! 
subroutine scan_air_ser
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i
	! jw

	! jw
	open(pw_air_ser, file=id_air_ser, form='formatted', status = 'old')

	! jw
	call skip_header_lines(pw_air_ser,id_air_ser)
	
	max_air_data_num = 1 ! jw
		
	read(pw_air_ser,*) i
	max_air_data_num = MAX(max_air_data_num,i) ! jw

	close(pw_air_ser)
end subroutine scan_air_ser