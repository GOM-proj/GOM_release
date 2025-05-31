!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Find: max_windp_data_num
!! 
subroutine scan_windp_ser
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, t
	integer :: i1, i2
	! jw
	
	! jw
	open(pw_windp_ser,file=id_windp_ser,form='formatted',status='old')
	
	! jw
	call skip_header_lines(pw_windp_ser, id_windp_ser)
	
	max_windp_data_num = 1 ! jw
	
	! jw
	do i=1,num_windp_ser
		read(pw_windp_ser,*) i1, i2
		max_windp_data_num = MAX(max_windp_data_num,i2)
		do t=1,i2
			read(pw_windp_ser,*) 
		end do
	end do
	close(pw_windp_ser)
end subroutine scan_windp_ser