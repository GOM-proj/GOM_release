!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Read hurricane_ser.inp,
!! for the space and time interpolation version
!! 
subroutine read_hurricane_ser_2
   use mod_global_variables
   use mod_file_definition
   implicit none
   
   integer :: i, j, n1, n2 
   real(dp) :: read_count
   real(dp) :: time_ratio
   real(dp) :: dx_hurricane_center, dy_hurricane_center ! , hurricane_center_x, hurricane_center_y
   ! end of local variables ==================================================!

   ! hurricane_read_interval is defined in read_hurricane.f90: 
   ! 		hurricane_read_interval = int(hurricane_dt / dt)
   read_count = mod(it,hurricane_read_interval)

	! wind_u0: left boundary value
	! wind_u1: middle (interpolated) value
	! wind_u2: right boundary value
   if(read_count == 0)then ! use exact given values
   	!$omp parallel do private(i)
      do i = 1, maxnod
      	! update the left boundary & interpolated values as the previous right boundary value
         wind_u0(i) = wind_u2(i)
         wind_v0(i) = wind_v2(i) 
         air_p0(i) = air_p2(i)
         
         wind_u1(i) = wind_u2(i)
         wind_v1(i) = wind_v2(i) 
         air_p1(i) = air_p2(i)
      end do
		!$omp end parallel do
		! update the current hurricane to the previous second one
      hurric_x_1 = hurric_x_2
      hurric_y_1 = hurric_y_2

		! update the current hurricane eye
      hurricane_center_x = hurric_x_1
      hurricane_center_y = hurric_y_1

      ! read the next hurricane center (eye) data ============================!
      read(pw_hurricane_center,*) hurric_year_2, hurric_month_2, hurric_day_2 , hurric_hour_2 , hurric_minute_2,     &
      &									 hurric_x_2   , hurric_y_2    , hurric_latitude_2, hurric_delta_pressure_2 , hurric_mwr_2
		
		! convert to julian day
      call julian(hurric_year_2, hurric_month_2, hurric_day_2, hurric_hour_2, hurric_minute_2, 0, &
      &				hurric_julian_day_2, hurric_julian_day_1900_2)

      hurricane_time_2 = hurric_julian_day_2
      
      ! read the next hurricane data =========================================!
      ! new right boundary values, from hurrican_ser.inp
      if(hurricane_data_type == 1) then	! hurricane_ser.inp is ascii
         do i = 1, maxnod
            read(pw_hurricane_ser,*) j, wind_u2(i), wind_v2(i), air_p2(i)
            air_p2(i) =  air_p2(i) * 100.0_dp
         end do
      else if(hurricane_data_type == 2) then	! hurricane_ser.inp is binary
      	do i=1,maxnod
	         read(pw_hurricane_ser) j, wind_u2(i),wind_v2(i),air_p2(i)
	         air_p2(i) =  air_p2(i) * 100.0_dp
	      end do
      end if

   	! use exact value from given data
   	!$omp parallel do private(j,n1,n2)
      do j = 1, maxface
         n1 = nodenum_at_face(1,j)
         n2 = nodenum_at_face(2,j)
         wind_u_at_face(j) = (wind_u1(n1) + wind_u1(n2)) * 0.5_dp
         wind_v_at_face(j) = (wind_v1(n1) + wind_v1(n2)) * 0.5_dp
      end do
      !$omp end parallel do
   else	! if(read_count /= 0)then
   	! time & space interpolation ===========================================!
   	! time_ration is the theta in my note (equation).
      time_ratio = (julian_day - hurricane_time_1) / (hurricane_time_2 - hurricane_time_1)
      
      ! note: dx & dy has directions (+ or -), thus they show the direction the hurricane is moving 
      dx_hurricane_center = hurric_x_2 - hurric_x_1
      dy_hurricane_center = hurric_y_2 - hurric_y_1
      
      ! Find expected hurricane center at the current time (linear interpolation)
      hurricane_center_x = hurric_x_1 + dx_hurricane_center * time_ratio
      hurricane_center_y = hurric_y_1 + dy_hurricane_center * time_ratio
      
      ! update the current hurricane location
      hurric_x_1 = hurricane_center_x
      hurric_y_1 = hurricane_center_y
      hurricane_time_1 = hurricane_time_1 + dt/86400.0	! [day]
		
		! write(*,*) hurricane_time_1, hurricane_time_2, julian_day, time_ratio
		
		! move and interpolate with the 1st snapshot ===========================!
		!$omp parallel do private(i)
      do i = 1, maxnod
      	! the moved grid points are at:
      	! 		from the first snapshot, it should moved forward as dx*(theta)
      	! this is the original method
         ! shiftx(i) = x_node(i) + (dx_hurricane_center*time_ratio) ! theta = (x0 - x1)/dx, and x0 = x1 + time_ratio*dx
         ! shifty(i) = y_node(i) + (dy_hurricane_center*time_ratio) !
         
         ! new backtracking approach, jw, check this again
         ! this is my method
         shiftx(i) = x_node(i) - (dx_hurricane_center*time_ratio) ! theta = (x0 - x1)/dx, and x0 = x1 + time_ratio*dx
         shifty(i) = y_node(i) - (dy_hurricane_center*time_ratio) !
         
			! write(*,*) x_node(1), dx_hurricane_center, shiftx(1)
			! write(*,*) y_node(1), dy_hurricane_center, shifty(1)
         ! stop 'aimmmm'
      end do
      !$omp end parallel do
      	
		! interpolate the wind field with distance at each node
      ! call IDW_3(maxnod, shiftx, shifty, wind_u1, wind_v1, wind_u1_new, wind_v1_new)
      ! call IDW_4(shiftx, shifty, wind_u1, wind_v1, air_p1, wind_u1_new, wind_v1_new, air_p1_new)
      call IDW_5(shiftx, shifty, wind_u1, wind_v1, air_p1, wind_u1_new, wind_v1_new, air_p1_new)
		
		! move and interpolate with the 2nd snapshot ===========================!
		!$omp parallel do private(i)
      do i = 1, maxnod
      	! the moved grid points are at:
      	! 		from the second snapshot, it should moved backward as dx*(1-theta)
      	! this is the original method      	
         shiftx(i) = x_node(i) - (dx_hurricane_center*(1.0_dp - time_ratio)) 
         shifty(i) = y_node(i) - (dy_hurricane_center*(1.0_dp - time_ratio))
			
         ! new backtracking approach, jw, check this again
         ! this is my method
         shiftx(i) = x_node(i) + (dx_hurricane_center*(1.0_dp - time_ratio)) 
         shifty(i) = y_node(i) + (dy_hurricane_center*(1.0_dp - time_ratio))
			! write(*,*) x_node(1), dx_hurricane_center, shiftx(1)
			! write(*,*) y_node(1), dy_hurricane_center, shifty(1)
         ! stop 'aimmmm'
      end do
      !$omp end parallel do

		! interpolate the wind field with distance at each node
      ! call IDW_3(maxnod, shiftx, shifty, wind_u2, wind_v2, wind_u2_new, wind_v2_new)
      ! call IDW_4(shiftx, shifty, wind_u2, wind_v2, air_p2, wind_u2_new, wind_v2_new, air_p2_new)
      call IDW_5(shiftx, shifty, wind_u2, wind_v2, air_p2, wind_u2_new, wind_v2_new, air_p2_new)
		
		! get new values at the current timestep ===============================!
		! linear time interpolation for winds
		!$omp parallel
		!$omp do private(i)
      do i = 1, maxnod
         wind_u1(i) = ((1.0 - time_ratio) * wind_u1_new(i)) + (time_ratio  * wind_u2_new(i))
         wind_v1(i) = ((1.0 - time_ratio) * wind_v1_new(i)) + (time_ratio  * wind_v2_new(i))
         air_p1(i)  = ((1.0 - time_ratio) * air_p1_new(i)) + (time_ratio  * air_p2_new(i))
      end do
		!$omp end do
		
		! calculate wind speed at faces
		!$omp do private(j,n1,n2)
      do j = 1, maxface
         n1 = nodenum_at_face(1,j)
         n2 = nodenum_at_face(2,j)
         wind_u_at_face(j) = (wind_u1(n1) + wind_u1(n2)) * 0.5
         wind_v_at_face(j) = (wind_v1(n1) + wind_v1(n2)) * 0.5
      end do
   	!$omp end do
   	!$omp end parallel
   end if   ! if(read_count /= 0.0 )then
end subroutine read_hurricane_ser_2
