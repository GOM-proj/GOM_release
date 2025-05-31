!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will find neighbour nodes for hurricane simulation
!! 
subroutine find_neighbour_nodes
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, k, ii
	real(dp):: dist_x, dist_y, min_dist, hurricane_search_range
	real(dp):: dist(maxnod)
	! jw
	allocate(hurricane_search_count(maxnod), hurricane_search_nodes(maxnod,maxnod))
	hurricane_search_count = 0
	hurricane_search_nodes = -999
	
	do i=1,maxnod
		do ii=1,maxnod
			if(i == ii) then
				dist(ii) = 1.0e15 ! jw
			else
				dist_x = x_node(i) - x_node(ii)
				dist_y = y_node(i) - y_node(ii)
				dist(ii) = sqrt(dist_x**2 + dist_y**2)
			end if
		end do
		
		min_dist = minval(dist)
		if(min_dist < 5000) then ! jw
			hurricane_search_range = 20*1000.0 ! jw
		else
			hurricane_search_range = min_dist * 4 ! jw
		end if
		
		! jw
		k = 0
		do ii=1,maxnod
			if(dist(ii) <= hurricane_search_range) then
				k = k+1
				hurricane_search_nodes(i,k) = ii
			end if
		end do
		hurricane_search_count(i) = k ! jw
	end do
	
	do i=1,maxnod
		if(hurricane_search_count(i) < 3) then
			write(*,*) i, hurricane_search_count(i), hurricane_search_nodes(i,1)
		end if
	end do
!	stop 'I amamam'
end subroutine find_neighbour_nodes