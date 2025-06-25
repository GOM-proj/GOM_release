!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Read hurricane_ser.inp,
!! for the linear interpolation version
!! 
subroutine read_hurricane_ser_1
   use mod_global_variables
   use mod_file_definition

   implicit none
   integer :: i, j, n1, n2
   integer :: read_count
   real(dp):: windp_increment
   ! end of local variables ==================================================!
   
   ! hurricane_read_interval is defined in read_hurricane.f90: 
   ! 		hurricane_read_interval = int(hurricane_dt / dt)
   ! read_count = mod(real(it),real(hurricane_read_interval))
	read_count = mod(it,hurricane_read_interval)
	
	! wind_u0: left boundary value
	! wind_u1: middle (interpolated) value
	! wind_u2: right boundary value
   if(read_count == 0) then
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
      	
      ! read next wind and pressure data, new right boundary values, from hurrican_ser.inp
      if(hurricane_data_type == 1) then	! hurricane_ser.inp is ascii
         do i = 1, maxnod
            read(pw_hurricane_ser,*) j, wind_u2(i), wind_v2(i), air_p2(i)
            air_p2(i) =  air_p2(i) * 100.0_dp	! Millibar to Pascal (N/m2)
         end do
      else if(hurricane_data_type == 2) then	! hurricane_ser.inp is binary
         do i = 1, maxnod
	         read(pw_hurricane_ser) j, wind_u2(i), wind_v2(i), air_p2(i)
            air_p2(i) =  air_p2(i) * 100.0_dp	! Millibar to Pascal (N/m2)
         end do
      end if
      
   	! use exact value from given data
   	!$omp parallel do private(j,n1,n2)
      do j = 1, maxface
         n1 = nodenum_at_face(1,j)
         n2 = nodenum_at_face(2,j)
         wind_u_at_face(j) = (wind_u1(n1) + wind_u1(n2)) * 0.5
         wind_v_at_face(j) = (wind_v1(n1) + wind_v1(n2)) * 0.5
      end do
      !$omp end parallel do
   else   ! if(read_count /= 0 )then
   	! time-only linear interpolation =======================================!
   	!$omp parallel 
   	!$omp do private(i,windp_increment)
      do i = 1, maxnod
      	windp_increment = (wind_u2(i) - wind_u0(i)) / real(hurricane_read_interval)
         wind_u1(i) = wind_u1(i) + windp_increment ! time interpolated wind

         windp_increment = (wind_v2(i) - wind_v0(i)) / real(hurricane_read_interval)
         wind_v1(i) = wind_v1(i) + windp_increment	! time interpolated wind

         windp_increment = (air_p2(i) - air_p0(i)) / real(hurricane_read_interval)
         air_p1(i) = air_p1(i) + windp_increment	! time interpolated airp
      end do
   	!$omp end do
   
      ! use interpolated data
      !$omp do private(j,n1,n2)
      do j = 1, maxface
         n1 = nodenum_at_face(1,j)
         n2 = nodenum_at_face(2,j)
         wind_u_at_face(j) = (wind_u1(n1) + wind_u1(n2)) * 0.5
         wind_v_at_face(j) = (wind_v1(n1) + wind_v1(n2)) * 0.5
      end do
    	!$omp end do
    	!$omp end parallel
   end if
end subroutine read_hurricane_ser_1
