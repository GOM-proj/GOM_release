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
   ! jw
   
	! jw
	if(num_windp_ser == 1) then 
		! jw
		! jw
		! jw
		u2 = julian_day*86400.0 ! jw
		
		
		! jw
		v12 = 0.0_dp
		v22 = 0.0_dp
		v32 = 0.0_dp
		
		do t=2,windp_ser_data_num(1)
			! jw
			! jw
			! jw
			! jw
			u1 = windp_ser_time(t-1,1) - reference_diff_days * 86400.0	! jw
			u3 = windp_ser_time(t  ,1) - reference_diff_days * 86400.0	! jw

			if(u2 >= u1 .and. u2 <= u3) then
				! jw
				v11 = windp_u(t-1,1)	! jw
				v13 = windp_u(t,1)	! jw
				
				! jw
				v21 = windp_v(t-1,1)	! jw
				v23 = windp_v(t,1)	! jw
				
				! jw
				v31 = windp_p(t-1,1)	! jw
				v33 = windp_p(t,1)	! jw
				
				! jw
				v12 = (v13 - v11)/(u3 - u1) * (u2 - u1) + v11	! jw
				v22 = (v23 - v21)/(u3 - u1) * (u2 - u1) + v21	! jw
				v32 = (v33 - v31)/(u3 - u1) * (u2 - u1) + v31	! jw
				
				exit
			end if
		end do
		
		! jw
		! jw
		! jw
		! jw
		! jw
		! jw
		
		! jw
		wind_u_at_face = v12
		wind_v_at_face = v22
		airp_at_face = v32		
   else if(num_windp_ser > 1) then ! jw
   	! jw
   	! jw
   	call IDW_windp_ser
   	
   	! jw
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
