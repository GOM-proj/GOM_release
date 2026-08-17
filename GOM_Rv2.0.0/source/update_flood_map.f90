!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine update_flood_map
	use mod_global_variables
	implicit none
	
	integer :: i
	! jw
	!$omp parallel do private(i)
	do i=1,maxnod
		if(eta_node(i) > max_eta_node(i)) then
			max_eta_node(i) = eta_node(i)
			max_flood_time(i) = julian_day
		end if
	end do
	!$omp end parallel do
end subroutine update_flood_map