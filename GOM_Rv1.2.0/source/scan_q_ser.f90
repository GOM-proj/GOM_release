!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine scan_q_ser
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, t, i1, i2
	! end of local variables --------------------------------------------------!

	! Read q_ser.inp ==========================================================!
	if(q_ser_shape == 1) then
		open(pw_q_ser, file=id_q_ser1, form='formatted', status = 'old')

		! skip header lines ----------------------------------------------------!
		call skip_header_lines(pw_q_ser,id_q_ser1)
		
		max_q_data_num = 1 ! to avoid allocation of 0 length array
		! max_q_ser = num_q_ser
		
		do i = 1, num_Q_ser
			read(pw_q_ser,*) i1, i2
			max_q_data_num = MAX(max_q_data_num,i2) ! update max_q_data_num
			do t=1,i2
				read(pw_q_ser,*)
			end do
		end do
	else if(q_ser_shape == 2) then
		open(pw_q_ser, file=id_q_ser2, form='formatted', status = 'old')

		! skip header lines ----------------------------------------------------!
		call skip_header_lines(pw_q_ser,id_q_ser2)
		
		max_q_data_num = 1 ! to avoid allocation of 0 length array
		! max_q_ser = num_q_ser
		
		read(pw_q_ser,*) i1, i2
		max_q_data_num = MAX(max_q_data_num,i2) ! update max_q_data_num		
	end if

	close(pw_q_ser)
end subroutine scan_q_ser