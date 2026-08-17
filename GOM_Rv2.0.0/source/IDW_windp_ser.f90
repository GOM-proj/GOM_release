!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Calculate wind speed and air pressure at each node from wind station data 
!! using the Inverse Distance Weighting (IDW) interpolation
!! Note that maximum 3 windp stations are allowed.
!! If you want to use more than 3 windp stations, consider to use wind model data.
!! 
subroutine IDW_windp_ser
   use mod_global_variables
   implicit none
	
	integer :: i, k, t
	integer :: power
   real(dp):: u1, u2, u3, v11, v12, v13, v21, v22, v23, v31, v32, v33
   real(dp):: x, y, x1, y1, r1, r2, r3
	real(dp):: dist(3), interp_wind_u(3), interp_wind_v(3), interp_p(3) ! jw
	! jw
   
   ! jw
   power = 1
	
	! jw
	dist = 0.0_dp
	interp_wind_u = 0.0_dp
	interp_wind_v = 0.0_dp
	interp_p = 0.0_dp
	
	! jw
	do k=1,num_windp_ser ! jw
		! jw
		u2 = julian_day*86400.0 ! jw
		
		! jw
		v12 = 0.0_dp
		v22 = 0.0_dp
		v32 = 0.0_dp

		do t=2,windp_ser_data_num(k)
			! jw
			! jw
			! jw
			! jw
			u1 = windp_ser_time(t-1,k) - reference_diff_days * 86400.0	! jw
			u3 = windp_ser_time(t  ,k) - reference_diff_days * 86400.0	! jw

			if(u2 >= u1 .and. u2 <= u3) then			
				! jw
				v11 = windp_u(t-1,k)	! jw
				v13 = windp_u(t,k)	! jw
				
				! jw
				v21 = windp_v(t-1,k)	! jw
				v23 = windp_v(t,k)	! jw
				
				! jw
				v31 = windp_p(t-1,k)	! jw
				v33 = windp_p(t,k)	! jw
					
				! jw
				v12 = (v13 - v11)/(u3 - u1) * (u2 - u1) + v11	! jw
				v22 = (v23 - v21)/(u3 - u1) * (u2 - u1) + v21	! jw
				v32 = (v33 - v31)/(u3 - u1) * (u2 - u1) + v31	! jw
				
				exit
			end if
		end do
		interp_wind_u(k) = v12 ! jw
		interp_wind_v(k) = v22 ! jw
		interp_p(k) = v32 ! jw
	end do

	! jw
	! jw
	! jw
	! jw
	if(num_windp_ser == 2) then
		!$omp parallel do private(i,k,x1,y1,dist,r1,r2)
		do i=1,maxnod
	      x = x_node(i)
	      y = y_node(i)
	 		
	 		! jw
	     	! jw
	      do k=1,num_windp_ser
	      	x1 = x_node(windp_station_node(k))
	      	y1 = y_node(windp_station_node(k))
	      	dist(k) = sqrt((x - x1)**2 + (y - y1)**2)
	      end do

         r1 = (dist(2))**power
         r2 = (dist(1))**power
         
         wind_u_at_node(i) = (r1*interp_wind_u(1) + r2*interp_wind_u(2)) / (r1 + r2)
         wind_v_at_node(i) = (r1*interp_wind_v(1) + r2*interp_wind_v(2)) / (r1 + r2)
         airp_at_node(i)   = (r1*interp_p(1) + r2*interp_p(2)) / (r1 + r2)
		end do
		!$omp end parallel do
	else if(num_windp_ser == 3) then ! jw
		!$omp parallel do private(i,k,x,y,dist,r1,r2,r3)
		do i=1,maxnod
	      x = x_node(i)
	      y = y_node(i)
	 		
	 		! jw
	     	! jw
	      do k=1,num_windp_ser
	      	x1 = x_node(windp_station_node(k))
	      	y1 = y_node(windp_station_node(k))
	      	dist(k) = sqrt((x - x1)**2 + (y - y1)**2)
	      end do

         r1 = (dist(2)*dist(3))**power
         r2 = (dist(1)*dist(3))**power
         r3 = (dist(1)*dist(2))**power

         wind_u_at_node(i) = (r1*interp_wind_u(1) + r2*interp_wind_u(2) + r3*interp_wind_u(3)) / (r1 + r2 + r3)
         wind_v_at_node(i) = (r1*interp_wind_v(1) + r2*interp_wind_v(2) + r3*interp_wind_v(3)) / (r1 + r2 + r3)
         airp_at_node(i)   = (r1*interp_p(1) + r2*interp_p(2) + r3*interp_p(3)) / (r1 + r2 + r3)
		end do
		!$omp end parallel do
	end if
end subroutine IDW_windp_ser