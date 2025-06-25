!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will read q_ser.inp one time before simulation start,
!! and the data will be stored as:
!!		q_ser_time(t,i)
!!		q_ser_Q(t,i)
!! 
subroutine read_q_ser
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, t
	real(dp) :: rtemp1, rtemp2
	real(dp) :: u1,u2,u3,v1,v2,v3	
	integer :: q_ser_id_dummy
	real(dp):: q_time_conv, q_time_adjust, q_unit_conv, q_adjust

	integer :: elapsed_t
	integer :: order, min_index, shift_num, ii, jj
	real(dp):: lagrange_sum, lagrange_product
	real(dp), allocatable :: x(:), y(:), difference(:)
	integer :: itemp
	real(dp),allocatable,dimension(:) :: q_buff
	! end of local variables --------------------------------------------------!

	write(pw_run_log,*) "	Read q_ser.inp"
	write(pw_run_log,*) "		Now, you are in 'read_input.f90 -> subroutine read_q_ser"

	! Read q_ser.inp ==========================================================!
	! Note:
	! 		q_interp_method should be an array, but I just use this now.
	! 		So, just use identical q_interp_method for each time series.
	if(q_ser_shape == 1) then
		! read q_ser1.inp
		open(pw_q_ser,file=id_q_ser1,form='formatted',status='old')
		
		! skip header lines
		call skip_header_lines(pw_q_ser,id_q_ser1)
		
		! read main body
		do i=1,num_q_ser
			read(pw_q_ser,*) q_ser_id_dummy, q_data_num(i), &
			&	q_time_conv, q_time_adjust, q_unit_conv, q_adjust, q_interp_method
			do t = 1,q_data_num(i)
				read(pw_q_ser,*) rtemp1, rtemp2
				q_ser_time(t,i) = (rtemp1 + q_time_adjust) * q_time_conv	! julian time [s]
				q_ser_Q(t,i) = (rtemp2 + q_adjust) * q_unit_conv		! [m3/s]
			end do
		end do
	else if(q_ser_shape == 2) then
		allocate(q_buff(max_q_ser))
		
		! read q_ser2.inp
		open(pw_q_ser,file=id_q_ser2,form='formatted',status='old')
		
		! skip header lines
		call skip_header_lines(pw_q_ser,id_q_ser2)
		
		! read main body
		read(pw_q_ser,*) q_ser_id_dummy, itemp, &
		&	q_time_conv, q_time_adjust, q_unit_conv, q_adjust
		do i=1,num_q_ser
			q_data_num(i) = itemp
		end do
		
		do t=1,itemp
			read(pw_q_ser,*) rtemp1, (q_buff(i), i=1,num_q_ser)
			
			do i=1,num_q_ser
				q_ser_time(t,i) = (rtemp1 + q_time_adjust) * q_time_conv ! julian time [s]
				q_ser_Q(t,i) = (q_buff(i) + q_adjust) * q_unit_conv ! [m3/s]
			end do
		end do
		deallocate(q_buff)		
	end if

	close(pw_q_ser)

	! write Q series checking file ============================================!
	if(num_Qb_cell > 0 .and. check_Q == 1) then
		open(pw_check_Q, file=id_check_Q, form='formatted', status = 'replace')
		write(pw_check_Q,*) 'Title = "Q time series"'
		write(pw_check_Q,'(A)',advance = 'no') 'Variables = "Time (day)"'
		do i=1,num_Qb_cell
			write(pw_check_Q,'(A,I3,A)',advance = 'no') ', "Q(',i,')[m3/s]"'	
		end do
		write(pw_check_Q,*)
		
		! show the interpolation result
		if(q_interp_method == 1) then	! Linear interpolation		
			do elapsed_t=0,int(ndt*dt),(30*60) ! every 30 minute
				u2 = elapsed_t	! = elapsed time, time for looking
				write(pw_check_Q,'(F20.10)',advance = 'no') u2/86400.0_dp
				do i=1,num_Qb_cell
					do t=2,q_data_num(Q_ser_id(i))
						u1 = q_ser_time(t-1,Q_ser_id(i))	! [s], 	lower bound time
						v1 = q_ser_Q(t-1,Q_ser_id(i))		! [m3/s],lower bound Q
						u3 = q_ser_time(t,Q_ser_id(i))	! [s], 	upper bound time
						v3 = q_ser_Q(t,Q_ser_id(i))		! [m3/s],upper bound Q
						
						if(u2 >= u1 .and. u2 <= u3) then
							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated Q
							write(pw_check_Q,'(F20.5)',advance = 'no') v2*Q_portion(i)
							exit
						end if
					end do
				end do
				write(pw_check_Q,*)
			end do
		else if(q_interp_method == 2) then	! Lagrange interpolation
			order = 2	! order of polynomial
			allocate(x(order+1), y(order+1), difference(max_q_data_num))
			x = 0.0_dp
			y = 0.0_dp
			
			do elapsed_t=0,int(ndt*dt),(30*60)	! every 30 minute
				u2 = elapsed_t	! = elapsed time, time for looking
				write(pw_check_Q,'(F20.10)',advance = 'no') u2/86400.0_dp
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
						min_index = min_index - shift_num;
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
					write(pw_check_Q,'(F20.5)',advance = 'no') lagrange_sum*Q_portion(i)
				end do
				write(pw_check_Q,*)
			end do
			deallocate(x,y,difference)					
		end if
		close(pw_check_Q)
	end if
	! End of Q series checking ================================================!
end subroutine read_q_ser
