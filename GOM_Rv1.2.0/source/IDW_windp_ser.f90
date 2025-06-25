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
	real(dp):: dist(3), interp_wind_u(3), interp_wind_v(3), interp_p(3) ! maximum three wind stations are allowed.
	! End of local variables ==================================================!
   
   ! weighting factor: normally 1 or 2
   power = 1
	
	! initialize variables
	dist = 0.0_dp
	interp_wind_u = 0.0_dp
	interp_wind_v = 0.0_dp
	interp_p = 0.0_dp
	
	! Linear interpolation in time ============================================!
	do k=1,num_windp_ser ! k'th station
		! u2 = time_begin*86400.0 + it*dt	! = elapsed time [s], current time and/or time for looking
		u2 = julian_day*86400.0 ! julian time [s], current time and/or time for looking
		
		! initialize interpolated values, otherwise it will use garbage values.
		v12 = 0.0_dp
		v22 = 0.0_dp
		v32 = 0.0_dp

		do t=2,windp_ser_data_num(k)
			! linear interpolation in time
			! Since windp_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
			! u1, and u3 should be shifted toward the simulation start year.
			! That is why I put "- reference_diff_days * 86400.0"			
			u1 = windp_ser_time(t-1,k) - reference_diff_days * 86400.0	! [s],	lower bound time
			u3 = windp_ser_time(t  ,k) - reference_diff_days * 86400.0	! [s],	upper bound time

			if(u2 >= u1 .and. u2 <= u3) then			
				! wind u component
				v11 = windp_u(t-1,k)	! [m/s], lower bound wind u
				v13 = windp_u(t,k)	! [m/s], upper bound wind u
				
				! wind v component
				v21 = windp_v(t-1,k)	! [m/s], lower bound wind v
				v23 = windp_v(t,k)	! [m/s], upper bound wind v
				
				! pressure component
				v31 = windp_p(t-1,k)	! [Millibar], lower bound air pressure
				v33 = windp_p(t,k)	! [Millibar], upper bound air pressure
					
				! linear interpolation in time 
				v12 = (v13 - v11)/(u3 - u1) * (u2 - u1) + v11	! wind u 
				v22 = (v23 - v21)/(u3 - u1) * (u2 - u1) + v21	! wind v
				v32 = (v33 - v31)/(u3 - u1) * (u2 - u1) + v31	! pressure
				
				exit
			end if
		end do
		interp_wind_u(k) = v12 ! time interpolated value at each station
		interp_wind_v(k) = v22 ! time interpolated value at each station
		interp_p(k) = v32 ! time interpolated value at each station
	end do

	! Space interpolation =====================================================!
	! Now, let's calculate wind speed at each node using calculated wind speed at each station at current time
	! Note that just three wind stations are allowed here.
	! Otherwise I will use wind model approach.	
	if(num_windp_ser == 2) then
		!$omp parallel do private(i,k,x1,y1,dist,r1,r2)
		do i=1,maxnod
	      x = x_node(i)
	      y = y_node(i)
	 		
	 		! Inverse Distance Weighting (IDW) Interpolation
	     	! calculate distance from the current node to each station     
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
	else if(num_windp_ser == 3) then ! currently, maximum three windp station is allowed.
		!$omp parallel do private(i,k,x,y,dist,r1,r2,r3)
		do i=1,maxnod
	      x = x_node(i)
	      y = y_node(i)
	 		
	 		! Inverse Distance Weighting (IDW) Interpolation
	     	! calculate distance from the current node to each station     
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