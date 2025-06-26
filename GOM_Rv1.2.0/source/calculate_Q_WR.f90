!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will cancluate interpolated river discharge, Q_add [m3/s],
!! at the current time step at each river boundary element
!! 
subroutine calculate_Q_WR
	use mod_global_variables
	
	implicit none
	integer :: i, t
	real(dp):: u1,v1,u2,v2,u3,v3

	integer :: order, min_index, shift_num, ii, jj
	real(dp):: lagrange_sum, lagrange_product
	real(dp), allocatable :: x(:), y(:), difference(:)
	!! end of local variables ------------------------------------------------!!
	
	! jw
	if(q_interp_method == 1) then	! jw
		! jw
		u2 = julian_day*86400.0 ! jw
		
		! jw
		v2 = 0.0_dp
		
		do i=1,num_WR_cell
			do t=2,q_data_num(WR_Q_ser_id(i))
				! jw
				! jw
				! jw
				u1 = q_ser_time(t-1,WR_Q_ser_id(i)) - reference_diff_days * 86400.0	! jw
				u3 = q_ser_time(t  ,WR_Q_ser_id(i)) - reference_diff_days * 86400.0	! jw
								
				if(u2 >= u1 .and. u2 <= u3) then
					v1 = q_ser_Q(t-1,WR_Q_ser_id(i))		! jw
					v3 = q_ser_Q(t  ,WR_Q_ser_id(i))		! jw
				
					v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
					Q_add_WR(i) = v2*WR_portion(i)
					exit
				end if
			end do
		end do
	else if(q_interp_method == 2) then	! jw
		order = 2	! jw
		allocate(x(order+1), y(order+1), difference(max_q_data_num))
		x = 0.0_dp
		y = 0.0_dp
		
		! jw
		u2 = julian_day*86400.0 ! jw
		do i=1,num_WR_cell
			! jw
			! jw
			difference = 9999999.0_dp	! jw
			do t=1,q_data_num(WR_Q_ser_id(i))
				difference(t) = q_ser_time(t,WR_Q_ser_id(i)) - u2
			end do
			min_index = minloc(difference,DIM=1,MASK=difference >= 0.0_dp)

			! jw
			if(min_index == 1) then
				min_index = min_index+1
			end if
			
			min_index = min_index - 1						
			shift_num = (min_index + order) - q_data_num(WR_Q_ser_id(i))
			if(shift_num > 0) then
				min_index = min_index - shift_num
			end if
			do ii=0,order
				x(ii+1) = q_ser_time(min_index + ii,WR_Q_ser_id(i))
				y(ii+1) = q_ser_Q(min_index + ii,WR_Q_ser_id(i))
			end do
			
			! jw
			lagrange_sum = 0.0_dp
			do ii=0,order
				lagrange_product = y(ii+1)
				do jj=0,order
					if(jj /= ii) then
						lagrange_product = lagrange_product * (u2-x(jj+1))/(x(ii+1)-x(jj+1))
					end if
				end do
				lagrange_sum = lagrange_sum + lagrange_product
			end do
			Q_add_WR(i) = lagrange_sum*WR_portion(i)
		end do
		deallocate(x,y,difference)		
	else if(q_interp_method == 3) then	! jw
		! jw
	end if
	
end subroutine calculate_Q_WR