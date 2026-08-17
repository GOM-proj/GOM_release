!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_dump3D
   use mod_global_variables
   use mod_file_definition   
	implicit none
	
	integer :: i, k
	character(len=40) :: format_string
	! jw
	
	! jw
	write(format_string,'(A,I0,A)') '(',1 + maxlayer*6 + 3,'E15.5E4)' ! jw
	if(IS3D_dump_binary == 0) then
		! jw
		! jw
		! jw
		write(pw_dump3D,'(F15.5)') julian_day * IS3D_dump_time_conv
		do i=1,maxnod
			! jw
			write(pw_dump3D,format_string) eta_node(i), &
			&	(u_node(k,i), k=0,maxlayer), &
			&	(v_node(k,i), k=0,maxlayer), &
			&	(w_node(k,i), k=0,maxlayer), &
			&	(salt_node(k,i), k=1,maxlayer), &
			&	(temp_node(k,i), k=1,maxlayer), &
			&	(rho_node(k,i),  k=1,maxlayer)
		end do
	else if(IS3D_dump_binary == 1) then
		! jw
		write(pw_dump3D) julian_day * IS3D_dump_time_conv
		do i=1,maxnod
			write(pw_dump3D) eta_node(i),   &
			&	(u_node(k,i), k=0,maxlayer), &
			&	(v_node(k,i), k=0,maxlayer), &
			&	(w_node(k,i), k=0,maxlayer), &
			&	(salt_node(k,i), k=1,maxlayer), &
			&	(temp_node(k,i), k=1,maxlayer), &
			&	(rho_node(k,i),  k=1,maxlayer)
		end do		
	end if
end subroutine write_dump3D