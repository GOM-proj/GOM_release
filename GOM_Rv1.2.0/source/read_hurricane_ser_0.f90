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
	! End of local variables ==================================================!
	   
   open(pw_hurricane_ser, file = id_hurricane_ser, form = 'formatted', status = 'old')

	! skip header lines 
	call skip_header_lines(pw_hurricane_ser,id_hurricane_ser)
	
	! read main body of the hurricane_ser.inp =================================!

	! calculate wind and pressure data reading interval
   hurricane_read_interval = int(hurricane_dt / dt) ! 36 (=3hr/5min)
	
	
	! First, we have to read two data sets since we have to use interpolated data at each elapsed_time.
   if(hurricane_data_type == 1) then	! hurricane_ser.inp is ascii
      ! read first set of wind and pressure data set from hurricane_ser.inp
      do i = 1, maxnod
         read(pw_hurricane_ser,*) serial_num, wind_u0(i), wind_v0(i), air_p0(i)
         air_p0(i) = air_p0(i) * 100.0_dp	! Millibar to Pascal (N/m2)
      end do
      
      ! read second set of wind and pressure data set from hurricane_ser.inp
      do i = 1, maxnod
         read(pw_hurricane_ser,*) serial_num, wind_u2(i), wind_v2(i), air_p2(i)
         air_p2(i) = air_p2(i) * 100.0_dp	! Millibar to Pascal (N/m2)
      end do
   elseif(hurricane_data_type == 2) then ! hurricane_ser.inp is binary
      ! read first set of wind and pressure data set from hurricane_ser.inp
      do i = 1, maxnod
         read(pw_hurricane_ser) serial_num, wind_u0(i), wind_v0(i), air_p0(i)
         air_p0(i) = air_p0(i) * 100.0_dp	! Millibar to Pascal (N/m2)
      end do
      
      ! read second set of wind and pressure data set from hurricane_ser.inp
      do i = 1, maxnod
         read(pw_hurricane_ser) serial_num, wind_u2(i), wind_v2(i), air_p2(i)
         air_p2(i) = air_p2(i) * 100.0_dp	! Millibar to Pascal (N/m2)
      end do
   end if

	! check wind interpolation method and open hurricane center data file, hurricane_center.inp
   if(hurricane_interp_method == 2) then
   	! Additional hurricane center location data is required: "hurricane_center.inp"
   	! reading the 1st hurricane center (eye) data
   	open(pw_hurricane_center,file=id_hurricane_center,form='formatted',status='old')
      read(pw_hurricane_center,*) ! header line
      read(pw_hurricane_center,*) hurric_year_1, hurric_month_1, hurric_day_1 , hurric_hour_1 , hurric_minute_1,     &
      &									 hurric_x_1   , hurric_y_1    , hurric_latitude_1, hurric_delta_pressure_1 , hurric_mwr_1
      
      call julian(hurric_year_1, hurric_month_1, hurric_day_1, hurric_hour_1, hurric_minute_1, 0, &
      &				hurric_julian_day_1, hurric_julian_day_1900_1)
			
      hurricane_time_1 = hurric_julian_day_1
		
		! current hurricane center
      hurricane_center_x  = hurric_x_1
      hurricane_center_y  = hurric_y_1

      ! reading the 2nd hurricane center (eye) data
      read(pw_hurricane_center,*) hurric_year_2, hurric_month_2, hurric_day_2 , hurric_hour_2 , hurric_minute_2,     &
      &									 hurric_x_2   , hurric_y_2    , hurric_latitude_2, hurric_delta_pressure_2 , hurric_mwr_2
     
      call julian(hurric_year_2, hurric_month_2, hurric_day_2, hurric_hour_2, hurric_minute_2, 0, &
      &				hurric_julian_day_2, hurric_julian_day_1900_2)
	
      hurricane_time_2 = hurric_julian_day_2
   end if
end subroutine read_hurricane_ser_0