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
	! end of local variables ==================================================!
	
	! replace water depth -----------------------------------------------------!
	do i=1,maxnod
		! h_node(i) = 11.0
		! if(h_node(i) > 1.0) then
		! 	h_node(i) = h_node(i) + 0.969
		! end if
		! h_node(i) = h_node(i) + 0.969
		! h_node(i) = h_node(i) + 2.0
		
		! update dry node corresponding to the given analytical depth
      ! if(h_node(i) <= 0.0) then
      !   initial_wetdry_node(i) = 1
      ! end if
      
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