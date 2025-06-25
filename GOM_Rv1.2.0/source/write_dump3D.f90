!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_dump3D
   use mod_global_variables
   use mod_file_definition   
	implicit none
	
	integer :: i, k
	character(len=40) :: format_string
	! End of local variables ==================================================!
	
	! write values at node for: eta, u, v, w, salt, temp, rho
	write(format_string,'(A,I0,A)') '(',1 + maxlayer*6 + 3,'E15.5E4)' ! ex, (10E15.5E4); ! here +1 for eta, +3 for uvw at bottom level
	if(IS3D_dump_binary == 0) then
		! ASCII format
		! Note:velocities are at vertical level, but salt and temp are at vertical layer
		! write(pw_dump3D,format_string) julian_day * IS3D_dump_time_conv
		write(pw_dump3D,'(F15.5)') julian_day * IS3D_dump_time_conv
		do i=1,maxnod
			! write(pw_dump3D,'(*(E15.5E4))') eta_node(i), &
			write(pw_dump3D,format_string) eta_node(i), &
			&	(u_node(k,i), k=0,maxlayer), &
			&	(v_node(k,i), k=0,maxlayer), &
			&	(w_node(k,i), k=0,maxlayer), &
			&	(salt_node(k,i), k=1,maxlayer), &
			&	(temp_node(k,i), k=1,maxlayer), &
			&	(rho_node(k,i),  k=1,maxlayer)
		end do
	else if(IS3D_dump_binary == 1) then
		! Binary format
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