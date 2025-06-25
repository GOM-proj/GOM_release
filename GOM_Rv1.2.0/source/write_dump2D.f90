!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_dump2D
   use mod_global_variables
   use mod_file_definition   
	implicit none
	
	integer :: i
	! End of local variables ==================================================!
	
	! write vertically averaged values at node for: eta, u, v, salt, and temp
	if(IS2D_dump_binary == 0) then
		! ASCII format
		write(pw_dump2D,'(F15.5)') julian_day * IS2D_dump_time_conv
		do i=1,maxnod
			write(pw_dump2D,'(6E15.5E4)') eta_node(i), ubar_node(i), vbar_node(i), sbar_node(i), tbar_node(i), rbar_node(i)
		end do
	else if(IS2D_dump_binary == 1) then
		! Bindary format
		write(pw_dump2D) julian_day * IS2D_dump_time_conv
		do i=1,maxnod
			write(pw_dump2D) eta_node(i), ubar_node(i), vbar_node(i), sbar_node(i), tbar_node(i), rbar_node(i)
		end do
	end if 	
end subroutine write_dump2D