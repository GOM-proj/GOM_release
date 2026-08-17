!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine reads only first two data sets from hurricane_ser.inp, and
!! the next data will be read from:
!! 	read_hurricane_ser_1.f90 for the linear interpolation version
!! 	read_hurricane_ser_2.f90 for the time & spatial interpolation version
!! 
subroutine read_hurricane_ser_0
	use mod_global_variables
	use mod_file_definition
	implicit none

	integer :: i
	integer :: serial_num
	! jw
	   
   open(pw_hurricane_ser, file = id_hurricane_ser, form = 'formatted', status = 'old')

	! jw
	call skip_header_lines(pw_hurricane_ser,id_hurricane_ser)
	
	! jw

	! jw
   hurricane_read_interval = int(hurricane_dt / dt) ! jw
	
	
	! jw
   if(hurricane_data_type == 1) then	! jw
      ! jw
      do i = 1, maxnod
         read(pw_hurricane_ser,*) serial_num, wind_u0(i), wind_v0(i), air_p0(i)
         air_p0(i) = air_p0(i) * 100.0_dp	! jw
      end do
      
      ! jw
      do i = 1, maxnod
         read(pw_hurricane_ser,*) serial_num, wind_u2(i), wind_v2(i), air_p2(i)
         air_p2(i) = air_p2(i) * 100.0_dp	! jw
      end do
   elseif(hurricane_data_type == 2) then ! jw
      ! jw
      do i = 1, maxnod
         read(pw_hurricane_ser) serial_num, wind_u0(i), wind_v0(i), air_p0(i)
         air_p0(i) = air_p0(i) * 100.0_dp	! jw
      end do
      
      ! jw
      do i = 1, maxnod
         read(pw_hurricane_ser) serial_num, wind_u2(i), wind_v2(i), air_p2(i)
         air_p2(i) = air_p2(i) * 100.0_dp	! jw
      end do
   end if

	! jw
   if(hurricane_interp_method == 2) then
   	! jw
   	! jw
   	open(pw_hurricane_center,file=id_hurricane_center,form='formatted',status='old')
      read(pw_hurricane_center,*) ! jw
      read(pw_hurricane_center,*) hurric_year_1, hurric_month_1, hurric_day_1 , hurric_hour_1 , hurric_minute_1,     &
      &									 hurric_x_1   , hurric_y_1    , hurric_latitude_1, hurric_delta_pressure_1 , hurric_mwr_1
      
      call julian(hurric_year_1, hurric_month_1, hurric_day_1, hurric_hour_1, hurric_minute_1, 0, &
      &				hurric_julian_day_1, hurric_julian_day_1900_1)
			
      hurricane_time_1 = hurric_julian_day_1
		
		! jw
      hurricane_center_x  = hurric_x_1
      hurricane_center_y  = hurric_y_1

      ! jw
      read(pw_hurricane_center,*) hurric_year_2, hurric_month_2, hurric_day_2 , hurric_hour_2 , hurric_minute_2,     &
      &									 hurric_x_2   , hurric_y_2    , hurric_latitude_2, hurric_delta_pressure_2 , hurric_mwr_2
     
      call julian(hurric_year_2, hurric_month_2, hurric_day_2, hurric_hour_2, hurric_minute_2, 0, &
      &				hurric_julian_day_2, hurric_julian_day_1900_2)
	
      hurricane_time_2 = hurric_julian_day_2
   end if
end subroutine read_hurricane_ser_0