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
	! End of local variables ==================================================!
	
	if(eta_ser_shape == 1) then
		! read eta_ser1.inp
		open(pw_eta_ser,file=id_eta_ser1,form='formatted',status='old')
		
		! skip header lines 
		call skip_header_lines(pw_eta_ser, id_eta_ser1)
		
		! read main body
		do i=1,num_eta_ser
			read(pw_eta_ser,*) eta_ser_id_dummy, eta_ser_data_num(i), &
			&	eta_time_conv, eta_time_adjust, eta_unit_conv, eta_adjust
			do t=1,eta_ser_data_num(i)
				read(pw_eta_ser,*) rtemp1, rtemp2
				eta_ser_time(t,i) = (rtemp1 + eta_time_adjust) * eta_time_conv 	! julian time [s]
				eta_ser_eta(t,i) = (rtemp2 + eta_adjust) * eta_unit_conv 			! [m]
				! write(*,*) eta_ser_time(t,i), eta_ser_eta(t,i)
			end do
		end do
	else if(eta_ser_shape == 2) then
		allocate(eta_buff(max_eta_ser))
		
		! read eta_ser2.inp
		open(pw_eta_ser,file=id_eta_ser2,form='formatted',status='old')
		
		! skip header lines 
		call skip_header_lines(pw_eta_ser, id_eta_ser2)
		
		! read main body
		read(pw_eta_ser,*) eta_ser_id_dummy, itemp, &
		&	eta_time_conv, eta_time_adjust, eta_unit_conv, eta_adjust
		do i=1,num_eta_ser
			eta_ser_data_num(i) = itemp
		end do
		
		do t=1,itemp
			read(pw_eta_ser,*) rtemp1, (eta_buff(i), i=1,num_eta_ser)
			
			do i=1,num_eta_ser
				eta_ser_time(t,i) = (rtemp1 + eta_time_adjust) * eta_time_conv 	! julian time [s]
				eta_ser_eta(t,i) = (eta_buff(i) + eta_adjust) * eta_unit_conv 		! [m]
			end do
		end do
		deallocate(eta_buff)
	end if
	
	close(pw_eta_ser)
	! stop 'jw'
end subroutine read_eta_ser