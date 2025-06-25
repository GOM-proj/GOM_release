!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine read_salt_ser
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, t
	integer :: salt_ser_id_dummy
	real(dp):: rtemp1, rtemp2, salt_time_conv, salt_time_adjust, salt_unit_conv, salt_adjust
	integer :: itemp
	real(dp),allocatable,dimension(:) :: salt_buff
	! End of local variables ==================================================!
	
	if(salt_ser_shape == 1) then
		! read salt_ser1.inp
		open(pw_salt_ser,file=id_salt_ser1,form='formatted',status='old')
		
		! skip header lines
		call skip_header_lines(pw_salt_ser,id_salt_ser1)
		
		! read main body
		do i=1,num_salt_ser
			read(pw_salt_ser,*) salt_ser_id_dummy, salt_ser_data_num(i), &
			&	salt_time_conv, salt_time_adjust, salt_unit_conv, salt_adjust
			do t = 1,salt_ser_data_num(i)
				read(pw_salt_ser,*) rtemp1, rtemp2
				salt_ser_time(t,i) = (rtemp1 + salt_time_adjust) * salt_time_conv	! julian time [s]
				salt_ser_salt(t,i) = (rtemp2 + salt_adjust) * salt_unit_conv		! [psu]
			end do
		end do
	else if(salt_ser_shape == 2) then
		allocate(salt_buff(max_salt_ser))
		
		! read salt_ser2.inp
		open(pw_salt_ser,file=id_salt_ser2,form='formatted',status='old')
		
		! skip header lines
		call skip_header_lines(pw_salt_ser,id_salt_ser2)
		
		! read main body
		read(pw_salt_ser,*) salt_ser_id_dummy, itemp, &
		&	salt_time_conv, salt_time_adjust, salt_unit_conv, salt_adjust
		do i=1,num_salt_ser
			salt_ser_data_num(i) = itemp
		end do
		
		do t=1,itemp
			read(pw_salt_ser,*) rtemp1, (salt_buff(i), i=1,num_salt_ser)
			
			do i=1,num_salt_ser
				salt_ser_time(t,i) = (rtemp1 + salt_time_adjust) * salt_time_conv ! julian time [s]
				salt_ser_salt(t,i) = (salt_buff(i) + salt_adjust) * salt_unit_conv ! [psu]
			end do
		end do
		deallocate(salt_buff)		
	end if
	
	close(pw_salt_ser)
end subroutine read_salt_ser