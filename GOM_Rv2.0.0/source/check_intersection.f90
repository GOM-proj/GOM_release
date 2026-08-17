!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! program to detect if two segments (1,2) and (3,4) have common points (i.e., intersection)
!! assumption: the 4 pts are distinctive.
!! The eqs. for the 2 lines are: 
!!		x=x1+(x2-x1)*tt1 or
!! 	x=x3+(x4-x3)*tt2
!! output: 
!! 	iflag: 0: no intersection or colinear
!! 	iflag: 1: exactly 1 intersection.
!! if iflag=1, (xin,yin) is the intersection.
!! ===========================================================================!
subroutine check_intersection(x1,x2,x3,x4, y1,y2,y3,y4, iflag,xin,yin)
	use mod_global_variables, only : dp
	implicit none
	
	real(dp), parameter 	:: small2=0.0_dp ! jw
	real(dp), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
	integer, intent(out) :: iflag
	real(dp), intent(out):: xin,yin
	real(dp):: tt1,tt2 ! jw
	
	real(dp) :: delta, delta1, delta2 
	! jw
	
	tt1 = -1000.0_dp
	tt2 = -1000.0_dp
	iflag = 0
	
	! jw
	! jw
	! jw
	delta  = (x2-x1)*(y3-y4)-(y2-y1)*(x3-x4) ! jw
	delta1 = (x3-x1)*(y3-y4)-(y3-y1)*(x3-x4) ! jw
	delta2 = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1) ! jw
	
	if(delta/=0.0d0) then
		tt1 = delta1/delta ! jw
		tt2 = delta2/delta ! jw
		if(tt1>=-small2 .and. tt1<=1+small2 .and. tt2>=-small2 .and. tt2<=1+small2) then
			! jw
			iflag = 1 ! jw

			! jw
			xin = x1+(x2-x1)*tt1
			yin = y1+(y2-y1)*tt1
		end if
	end if
end subroutine check_intersection

