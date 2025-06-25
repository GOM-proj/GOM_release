!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will cancluate interpolated river discharge, Q_add [m3/s],
!! at the current time step at each river boundary element
!! 
subroutine calculate_Q
	use mod_global_variables
	
	implicit none
	integer :: i, t
	real(dp):: u1,v1,u2,v2,u3,v3

	integer :: order, min_index, shift_num, ii, jj
	real(dp):: lagrange_sum, lagrange_product
	real(dp), allocatable :: x(:), y(:), difference(:)
	!! end of local variables ------------------------------------------------!!
	
	! note: q_interp_method should be an array... so, I need to update this later
	if(q_interp_method == 1) then	! Linear interpolation
		! u2 = time_begin*86400.0 + it*dt	! = elapsed time [s], current time and/or time for looking
		u2 = julian_day*86400.0 ! julian time [s], current time and/or time for looking
		
		! initialize interpolated values, otherwise it will use garbage values.
		v2 = 0.0_dp
		
		do i=1,num_Qb_cell
			do t=2,q_data_num(Q_ser_id(i))
				! Since q_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
				! u1, and u3 should be shifted toward the simulation start year.
				! That is why I put "- reference_diff_days * 86400.0"				
				u1 = q_ser_time(t-1,Q_ser_id(i)) - reference_diff_days * 86400.0	! [s],	lower bound time
				u3 = q_ser_time(t  ,Q_ser_id(i)) - reference_diff_days * 86400.0	! [s],	upper bound time
								
				if(u2 >= u1 .and. u2 <= u3) then
					v1 = q_ser_Q(t-1,Q_ser_id(i))		! [m3/s],lower bound Q
					v3 = q_ser_Q(t,Q_ser_id(i))		! [m3/s],upper bound Q
				
					v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated Q
					Q_add(i) = v2*Q_portion(i)
					exit
				end if
			end do
		end do
	else if(q_interp_method == 2) then	! Lagrange interpolation
		order = 2	! order of polynomial
		allocate(x(order+1), y(order+1), difference(max_q_data_num))
		x = 0.0_dp
		y = 0.0_dp
		
		! u2 = time_begin*86400.0 + it*dt	! = elapsed time [s], current time and/or time for looking
		u2 = julian_day*86400.0 ! julian time [s], current time and/or time for looking
		do i=1,num_Qb_cell
			! re-construct data
			! select only order+1 points letting nearest point to the 0th order point
			difference = 9999999.0_dp	! initialize difference as a big number.
			do t=1,q_data_num(q_ser_id(i))
				difference(t) = q_ser_time(t,q_ser_id(i)) - u2
			end do
			min_index = minloc(difference,DIM=1,MASK=difference >= 0.0_dp)

			! at the first data point, shift to the next point
			if(min_index == 1) then
				min_index = min_index+1
			end if
			
			min_index = min_index - 1						
			shift_num = (min_index + order) - q_data_num(q_ser_id(i))
			if(shift_num > 0) then
				min_index = min_index - shift_num
			end if
			do ii=0,order
				x(ii+1) = q_ser_time(min_index + ii,q_ser_id(i))
				y(ii+1) = q_ser_Q(min_index + ii,q_ser_id(i))
			end do
			
			! Lagrange interpolation with re-constructed data points
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
			Q_add(i) = lagrange_sum*Q_portion(i)					
		end do
		deallocate(x,y,difference)		
	else if(q_interp_method == 3) then	! Cubic Spline interpolation
		! not yet included
	end if
	
end subroutine calculate_Q