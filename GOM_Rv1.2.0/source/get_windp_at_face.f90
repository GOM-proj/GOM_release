!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Get wind speed and air pressure at each face at current time
!! 
subroutine get_windp_at_face
   use mod_global_variables
   use mod_file_definition   
   implicit none
   
   integer :: j, t
   integer :: n1, n2
   real(dp):: u1, u2, u3, v11, v12, v13, v21, v22, v23, v31, v32, v33
   ! End of local variables ==================================================!
   
	! calculate new wind speed and air pressure at each face at current time ==!
	if(num_windp_ser == 1) then 
		! if just one station data is given, 
		! linearly interpolated (in time) uniform wind and air pressure information will be used in all grids.
		! u2 = time_begin*86400.0 + it*dt	! = elapsed time [s], current time and/or time for looking
		u2 = julian_day*86400.0 ! julian time [s], current time and/or time for looking
		
		
		! initialize interpolated values, otherwise it will use garbage values.
		v12 = 0.0_dp
		v22 = 0.0_dp
		v32 = 0.0_dp
		
		do t=2,windp_ser_data_num(1)
			! linear interpolation in time
			! Since windp_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
			! u1, and u3 should be shifted toward the simulation start year.
			! That is why I put "- reference_diff_days * 86400.0"
			u1 = windp_ser_time(t-1,1) - reference_diff_days * 86400.0	! [s],	lower bound time
			u3 = windp_ser_time(t  ,1) - reference_diff_days * 86400.0	! [s],	upper bound time

			if(u2 >= u1 .and. u2 <= u3) then
				! wind u component (v1*)
				v11 = windp_u(t-1,1)	! [m/s], lower bound wind u
				v13 = windp_u(t,1)	! [m/s], upper bound wind u
				
				! wind v component (v2*)
				v21 = windp_v(t-1,1)	! [m/s], lower bound wind v
				v23 = windp_v(t,1)	! [m/s], upper bound wind v
				
				! pressure component (v3*)
				v31 = windp_p(t-1,1)	! [Millibar], lower bound air pressure
				v33 = windp_p(t,1)	! [Millibar], upper bound air pressure
				
				! linear interpolation in time 
				v12 = (v13 - v11)/(u3 - u1) * (u2 - u1) + v11	! wind u 
				v22 = (v23 - v21)/(u3 - u1) * (u2 - u1) + v21	! wind v
				v32 = (v33 - v31)/(u3 - u1) * (u2 - u1) + v31	! pressure
				
				exit
			end if
		end do
		
		! give identical wind speed across all domain
		! do j=1,maxface
		! 	wind_u_at_face(j) = v12
		! 	wind_v_at_face(j) = v22
		! 	airp_at_face(j) = v32
		! end do
		
		! this is equivalent to above do loop
		wind_u_at_face = v12
		wind_v_at_face = v22
		airp_at_face = v32		
   else if(num_windp_ser > 1) then ! multiple wind stations are given
   	! this will calculate wind velocity and pressure at each node: 
   	! 		[wind_u_at_node], [wind_v_at_node], and [airp_at_node]
   	call IDW_windp_ser
   	
   	! calculate wind velocity and air pressure at each face using calculated wind speed at each node
   	!$omp parallel do private(j,n1,n2)
      do j = 1, maxface
         n1 = nodenum_at_face(1,j)
         n2 = nodenum_at_face(2,j)
         
         wind_u_at_face(j) = (wind_u_at_node(n1) + wind_u_at_node(n2)) * 0.5
         wind_v_at_face(j) = (wind_v_at_node(n1) + wind_v_at_node(n2)) * 0.5
         airp_at_face(j) = (airp_at_node(n1) + airp_at_node(n2)) * 0.5
      end do
   end if 
end subroutine get_windp_at_face
