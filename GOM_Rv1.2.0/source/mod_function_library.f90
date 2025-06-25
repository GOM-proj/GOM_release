!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
module mod_function_library
   implicit none
	
	contains
	
	! =========================================================================!
	! jw, thise function has never been used. =================================!
	function lindex(node,i_element)
	   use mod_global_variables
	   use mod_file_definition
	
	   implicit none
	
	   integer :: lindex
	   integer,intent(in) :: node, i_element
	   integer :: j
	
	   lindex = 0
	   do j = 1, 3 
	      if(node == nodenum_at_cell(j,i_element))then
	         lindex = j
	      endif
	   enddo
	   if(lindex == 0)then
	     write(pw_run_log,*)node,' is not in element ', i_element
	   	stop
	   end if
	end function lindex
	
	! =========================================================================!
	! jw, this function has never been used. ==================================!
	function sums(x1,x2,x3,x4,y1,y2,y3,y4)
	   implicit none
	   real(8) :: sums
	   real(8),intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
	
	
	   sums = dabs( (x4-x3)*(y2-y3) + (x2-x3)*(y3-y4) ) / 2.0d0   &
	     &  + dabs( (x4-x1)*(y3-y1) - (y4-y1)*(x3-x1) ) / 2.0d0   &
	     &  + dabs( (y4-y1)*(x2-x1) - (x4-x1)*(y2-y1) ) / 2.0d0
	end function sums
	
	! =========================================================================!
	! calculate triangle area by assigned 3 nodes =============================!
	! i.e., (x1,y1), (x2,y2), (x3,y3)
	function calculate_area(x1,x2,x3,y1,y2,y3)
		use mod_global_variables, only : dp
	   implicit none
	
	   real(dp) :: calculate_area
	   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
	   ! end of local variables ===============================================!
	   
	   ! calculate_area = (((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3)) / 2.0_dp) ! jw, this was the original one.
	   calculate_area = (abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3)) / 2.0_dp) ! jw, this correct one, but the following equation is better.
	   ! calculate_area = 0.5_dp * abs(x1(y2-y3) + x2(y3-y1) + x3(y1-y2)) ! jw, new correct calculation; I don't understand why this does not work...
	end function calculate_area      
	
	! =========================================================================!
	! calculate the coriolis parameter from the latitude(in degrees) ==========!
	function coriolis_from_lat_deg(latitude)
		use mod_global_variables, only : dp
	   implicit none 
	
	   real(dp) :: coriolis_from_lat_deg
	   real(dp),parameter  :: pi1 = 3.141592654_dp, omega1=2.0_dp*pi1/86400.0_dp
	   real(dp),intent(in) :: latitude
	   ! end of local variables ===============================================!
	   
	   coriolis_from_lat_deg = 2.0_dp*omega1*dsin(latitude*pi1/180.0_dp) 
	end function coriolis_from_lat_deg 
end module mod_function_library

