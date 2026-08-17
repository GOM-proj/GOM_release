!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine read_eta_ser
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, t
	integer :: eta_ser_id_dummy
	real(dp):: rtemp1, rtemp2, eta_time_conv, eta_time_adjust, eta_unit_conv, eta_adjust
	integer :: itemp
	real(dp),allocatable,dimension(:) :: eta_buff
	! jw
	
	if(eta_ser_shape == 1) then
		! jw
		open(pw_eta_ser,file=id_eta_ser1,form='formatted',status='old')
		
		! jw
		call skip_header_lines(pw_eta_ser, id_eta_ser1)
		
		! jw
		do i=1,num_eta_ser
			read(pw_eta_ser,*) eta_ser_id_dummy, eta_ser_data_num(i), &
			&	eta_time_conv, eta_time_adjust, eta_unit_conv, eta_adjust
			do t=1,eta_ser_data_num(i)
				read(pw_eta_ser,*) rtemp1, rtemp2
				eta_ser_time(t,i) = (rtemp1 + eta_time_adjust) * eta_time_conv 	! jw
				eta_ser_eta(t,i) = (rtemp2 + eta_adjust) * eta_unit_conv 			! jw
				! jw
			end do
		end do
	else if(eta_ser_shape == 2) then
		allocate(eta_buff(max_eta_ser))
		
		! jw
		open(pw_eta_ser,file=id_eta_ser2,form='formatted',status='old')
		
		! jw
		call skip_header_lines(pw_eta_ser, id_eta_ser2)
		
		! jw
		read(pw_eta_ser,*) eta_ser_id_dummy, itemp, &
		&	eta_time_conv, eta_time_adjust, eta_unit_conv, eta_adjust
		do i=1,num_eta_ser
			eta_ser_data_num(i) = itemp
		end do
		
		do t=1,itemp
			read(pw_eta_ser,*) rtemp1, (eta_buff(i), i=1,num_eta_ser)
			
			do i=1,num_eta_ser
				eta_ser_time(t,i) = (rtemp1 + eta_time_adjust) * eta_time_conv 	! jw
				eta_ser_eta(t,i) = (eta_buff(i) + eta_adjust) * eta_unit_conv 		! jw
			end do
		end do
		deallocate(eta_buff)
	end if
	
	close(pw_eta_ser)
	! jw
end subroutine read_eta_ser