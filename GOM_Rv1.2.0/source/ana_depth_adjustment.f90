!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Re-define water depth analytically (i.e., overwrite initial water depth obtained from node.inp)
!! 
subroutine ana_depth_adjustment
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i
	! jw
	
	! jw
	do i=1,maxnod
		! jw
		! jw
		! jw
		! jw
		! jw
		! jw
		
		! jw
      ! jw
      ! jw
      ! jw
      
      if(i >= 15755 .and. i <= 15768) then
      	if(h_node(i) > 5.0) then
      		h_node(i) = 5.0
      	end if
      else if(i >= 15627 .and. i <= 15642) then
      	if(h_node(i) > 5.0) then
      		h_node(i) = 5.0
      	end if
      else if(i >= 15500 .and. i <= 15514) then
      	if(h_node(i) > 5.0) then
      		h_node(i) = 5.0
      	end if
      else if(i >= 15401 .and. i <= 15413) then
      	if(h_node(i) > 5.0) then
      		h_node(i) = 5.0
      	end if
      else if(i >= 15350 .and. i <= 15362) then	
      	if(h_node(i) > 5.0) then
      		h_node(i) = 5.0
      	end if
      else if(i >= 15303 .and. i <= 15316) then
      	if(h_node(i) > 5.0) then
      		h_node(i) = 5.0
      	end if
   	end if
      
      
	end do	
end subroutine ana_depth_adjustment