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
   real(dp) :: dx_hurricane_center, dy_hurricane_center ! jw
   ! jw

   ! jw
   ! jw
   read_count = mod(it,hurricane_read_interval)

	! jw
	! jw
	! jw
   if(read_count == 0)then ! jw
   	!$omp parallel do private(i)
      do i = 1, maxnod
      	! jw
         wind_u0(i) = wind_u2(i)
         wind_v0(i) = wind_v2(i) 
         air_p0(i) = air_p2(i)
         
         wind_u1(i) = wind_u2(i)
         wind_v1(i) = wind_v2(i) 
         air_p1(i) = air_p2(i)
      end do
		!$omp end parallel do
		! jw
      hurric_x_1 = hurric_x_2
      hurric_y_1 = hurric_y_2

		! jw
      hurricane_center_x = hurric_x_1
      hurricane_center_y = hurric_y_1

      ! jw
      read(pw_hurricane_center,*) hurric_year_2, hurric_month_2, hurric_day_2 , hurric_hour_2 , hurric_minute_2,     &
      &									 hurric_x_2   , hurric_y_2    , hurric_latitude_2, hurric_delta_pressure_2 , hurric_mwr_2
		
		! jw
      call julian(hurric_year_2, hurric_month_2, hurric_day_2, hurric_hour_2, hurric_minute_2, 0, &
      &				hurric_julian_day_2, hurric_julian_day_1900_2)

      hurricane_time_2 = hurric_julian_day_2
      
      ! jw
      ! jw
      if(hurricane_data_type == 1) then	! jw
         do i = 1, maxnod
            read(pw_hurricane_ser,*) j, wind_u2(i), wind_v2(i), air_p2(i)
            air_p2(i) =  air_p2(i) * 100.0_dp
         end do
      else if(hurricane_data_type == 2) then	! jw
      	do i=1,maxnod
	         read(pw_hurricane_ser) j, wind_u2(i),wind_v2(i),air_p2(i)
	         air_p2(i) =  air_p2(i) * 100.0_dp
	      end do
      end if

   	! jw
   	!$omp parallel do private(j,n1,n2)
      do j = 1, maxface
         n1 = nodenum_at_face(1,j)
         n2 = nodenum_at_face(2,j)
         wind_u_at_face(j) = (wind_u1(n1) + wind_u1(n2)) * 0.5_dp
         wind_v_at_face(j) = (wind_v1(n1) + wind_v1(n2)) * 0.5_dp
      end do
      !$omp end parallel do
   else	! jw
   	! jw
   	! jw
      time_ratio = (julian_day - hurricane_time_1) / (hurricane_time_2 - hurricane_time_1)
      
      ! jw
      dx_hurricane_center = hurric_x_2 - hurric_x_1
      dy_hurricane_center = hurric_y_2 - hurric_y_1
      
      ! jw
      hurricane_center_x = hurric_x_1 + dx_hurricane_center * time_ratio
      hurricane_center_y = hurric_y_1 + dy_hurricane_center * time_ratio
      
      ! jw
      hurric_x_1 = hurricane_center_x
      hurric_y_1 = hurricane_center_y
      hurricane_time_1 = hurricane_time_1 + dt/86400.0	! jw
		
		! jw
		
		! jw
		!$omp parallel do private(i)
      do i = 1, maxnod
      	! jw
      	! jw
      	! jw
         ! jw
         ! jw
         
         ! jw
         ! jw
         shiftx(i) = x_node(i) - (dx_hurricane_center*time_ratio) ! jw
         shifty(i) = y_node(i) - (dy_hurricane_center*time_ratio) !
         
			! jw
			! jw
         ! jw
      end do
      !$omp end parallel do
      	
		! jw
      ! jw
      ! jw
      call IDW_5(shiftx, shifty, wind_u1, wind_v1, air_p1, wind_u1_new, wind_v1_new, air_p1_new)
		
		! jw
		!$omp parallel do private(i)
      do i = 1, maxnod
      	! jw
      	! jw
      	! jw
         shiftx(i) = x_node(i) - (dx_hurricane_center*(1.0_dp - time_ratio)) 
         shifty(i) = y_node(i) - (dy_hurricane_center*(1.0_dp - time_ratio))
			
         ! jw
         ! jw
         shiftx(i) = x_node(i) + (dx_hurricane_center*(1.0_dp - time_ratio)) 
         shifty(i) = y_node(i) + (dy_hurricane_center*(1.0_dp - time_ratio))
			! jw
			! jw
         ! jw
      end do
      !$omp end parallel do

		! jw
      ! jw
      ! jw
      call IDW_5(shiftx, shifty, wind_u2, wind_v2, air_p2, wind_u2_new, wind_v2_new, air_p2_new)
		
		! jw
		! jw
		!$omp parallel
		!$omp do private(i)
      do i = 1, maxnod
         wind_u1(i) = ((1.0 - time_ratio) * wind_u1_new(i)) + (time_ratio  * wind_u2_new(i))
         wind_v1(i) = ((1.0 - time_ratio) * wind_v1_new(i)) + (time_ratio  * wind_v2_new(i))
         air_p1(i)  = ((1.0 - time_ratio) * air_p1_new(i)) + (time_ratio  * air_p2_new(i))
      end do
		!$omp end do
		
		! jw
		!$omp do private(j,n1,n2)
      do j = 1, maxface
         n1 = nodenum_at_face(1,j)
         n2 = nodenum_at_face(2,j)
         wind_u_at_face(j) = (wind_u1(n1) + wind_u1(n2)) * 0.5
         wind_v_at_face(j) = (wind_v1(n1) + wind_v1(n2)) * 0.5
      end do
   	!$omp end do
   	!$omp end parallel
   end if   ! jw
end subroutine read_hurricane_ser_2
