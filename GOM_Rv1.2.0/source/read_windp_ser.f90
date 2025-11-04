!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will read windp_ser.inp one time before simulation start,
!! and the data will be stored as:
!!		windp_ser_time(t,i)
!!		windp_u(t,i)
!!		windp_v(t,i)
!!		windp_p(t,i)
!! 
subroutine read_windp_ser
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, t
	real(dp):: windp_time_conv, windp_time_adjust, windp_w_unit_conv, windp_p_unit_conv
	real(dp):: temp1, temp2, temp3, temp4
	integer :: windp_ser_id_dummy
	! jw
	
	! jw
	open(pw_windp_ser,file=id_windp_ser,form='formatted',status='old')
	
	! jw
	call skip_header_lines(pw_windp_ser, id_windp_ser)

	! jw
	do i=1,num_windp_ser
		read(pw_windp_ser,*) windp_ser_id_dummy, windp_ser_data_num(i), windp_time_conv, windp_time_adjust, &
		&	windp_w_unit_conv, windp_p_unit_conv, windp_station_node(i)
		do t=1,windp_ser_data_num(i)
			read(pw_windp_ser,*) temp1, temp2, temp3, temp4 ! jw
			windp_ser_time(t,i) = (temp1 + windp_time_adjust) * windp_time_conv ! jw
			windp_u(t,i) = temp2*windp_w_unit_conv
			windp_v(t,i) = temp3*windp_w_unit_conv
			windp_p(t,i) = temp4*windp_p_unit_conv
		end do
	end do
	close(pw_windp_ser)
end subroutine read_windp_ser