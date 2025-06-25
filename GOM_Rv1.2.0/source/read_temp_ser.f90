!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine read_temp_ser
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, t
	integer :: temp_ser_id_dummy
	real(dp):: rtemp1, rtemp2, temp_time_conv, temp_time_adjust, temp_unit_conv, temp_adjust
	integer :: itemp
	real(dp),allocatable,dimension(:) :: temp_buff
	! End of local variables ==================================================!
	
	if(temp_ser_shape == 1) then
		! read temp_ser1.inp
		open(pw_temp_ser,file=id_temp_ser1,form='formatted',status='old')
		
		! skip header lines
		call skip_header_lines(pw_temp_ser,id_temp_ser1)
		
		! read main body
		do i=1,num_temp_ser
			read(pw_temp_ser,*) temp_ser_id_dummy, temp_ser_data_num(i), &
			&	temp_time_conv, temp_time_adjust, temp_unit_conv, temp_adjust
			do t = 1,temp_ser_data_num(i)
				read(pw_temp_ser,*) rtemp1, rtemp2
				temp_ser_time(t,i) = (rtemp1 + temp_time_adjust) * temp_time_conv	! julian time [s]
				temp_ser_temp(t,i) = (rtemp2 + temp_adjust) * temp_unit_conv		! [C]
				! write(*,*) temp_ser_temp(t,i)
			end do
		end do
	else if(temp_ser_shape == 2) then
		allocate(temp_buff(max_temp_ser))
		
		! read temp_ser2.inp
		open(pw_temp_ser,file=id_temp_ser2,form='formatted',status='old')
		
		! skip header lines
		call skip_header_lines(pw_temp_ser,id_temp_ser2)
		
		! read main body
		read(pw_temp_ser,*) temp_ser_id_dummy, itemp, &
		&	temp_time_conv, temp_time_adjust, temp_unit_conv, temp_adjust
		do i=1,num_temp_ser
			temp_ser_data_num(i) = itemp
		end do
		
		do t=1,itemp
			read(pw_temp_ser,*) rtemp1, (temp_buff(i), i=1,num_temp_ser)
			
			do i=1,num_temp_ser
				temp_ser_time(t,i) = (rtemp1 + temp_time_adjust) * temp_time_conv ! julian time [s]
				temp_ser_temp(t,i) = (temp_buff(i) + temp_adjust) * temp_unit_conv ! [C]
				
				! write(*,*) temp_ser_temp(t,i)
			end do
		end do
		deallocate(temp_buff)		
	end if
	
	close(pw_temp_ser)
end subroutine read_temp_ser