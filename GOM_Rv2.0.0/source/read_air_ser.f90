!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine read_air_ser
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, t
	integer :: air_data_num
	real(dp):: rtemp1, air_time_conv, air_time_adjust
	! jw
	
	! jw
	open(pw_air_ser,file=id_air_ser,form='formatted',status='old')
	
	! jw
	call skip_header_lines(pw_air_ser,id_air_ser)
	
	! jw
	read(pw_air_ser,*) air_data_num, air_time_conv, air_time_adjust
	
	do t=1,air_data_num
		read(pw_air_ser,*) rtemp1, T_air(t), T_dew(t), air_wind_speed(t), solar_radiation(t), cloud(t), rain(t), evaporation(t)
		air_ser_time(t) = (rtemp1 + air_time_adjust) * air_time_conv ! jw
	end do
	
	close(pw_air_ser)
end subroutine read_air_ser