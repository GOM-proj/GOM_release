!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Inverse Distance Weighting (IDW) Interpolation - two inputs and outputs version:
!!		xloc    (input)   : data location x                
!!		yloc    (input)   : data location y                
!!		var_old1(input)   : wind_u
!!		var_old2(input)   : wind_v
!!		var_old3(input)	: airp
!!		var_new1(output)  : solution, interpolated value, wind_u
!!		var_new2(output)  : solution, interpolated value, wind_v
!!		var_new3(output)  : solution, interpolated value, airp
!! 
!! We will use only three closest points for the spatial interpolation
!! We also use power = 1
!! 
!! This version will work with read_hurricane_ser_2.f90
!! 
subroutine IDW_5(xloc, yloc, var_old1, var_old2, var_old3, var_new1, var_new2, var_new3)
   use mod_global_variables
   implicit none

   real(dp), dimension(maxnod), intent(in) :: xloc, yloc, var_old1, var_old2, var_old3
	real(dp), dimension(maxnod), intent(inout) :: var_new1, var_new2, var_new3
   
   integer :: i, k, num_stations, n1, n2, n3, ii
   integer :: power
   real(dp):: r1, r2, r3
   real(dp), dimension(3) :: min_dist
  	real(dp), allocatable :: dist(:)
   real(dp):: dist_x, dist_y
   real(dp):: xn, yn
   ! jw

   ! jw
   power = 1
	
	! jw
	! jw
	num_stations = 3

	! jw
	!$omp parallel do private(i,k,min_dist,dist,xn,yn,dist_x,dist_y,ii,n1,n2,n3,r1,r2,r3)
   do i = 1, maxnod
      ! jw
      min_dist(1) = 1.0e15	! jw
      min_dist(2) = 1.0e15 	! jw
      min_dist(3) = 1.0e15	! jw
      
      allocate(dist(hurricane_search_count(i)))
   	dist = 1.0e15
   	
   	! jw
   	xn = xloc(i)
   	yn = yloc(i)
   	
   	! jw
   	do k=1,hurricane_search_count(i)
   		dist_x = xn - x_node(hurricane_search_nodes(i,k))
   		dist_y = yn - y_node(hurricane_search_nodes(i,k))
   		dist(k) = sqrt(dist_x**2 + dist_y**2)
   	end do
   	
   	! jw
   	ii = minloc(dist,dim=1)
   	n1 = hurricane_search_nodes(i,ii)
   	min_dist(1) = dist(ii)
   	
   	! jw
   	dist(ii) = 1.0e15
   	ii = minloc(dist,dim=1)
   	n2 = hurricane_search_nodes(i,ii)
   	min_dist(2) = dist(ii)

   	! jw
   	dist(ii) = 1.0e15
   	ii = minloc(dist,dim=1)
   	n3 = hurricane_search_nodes(i,ii)
   	min_dist(3) = dist(ii)
   	
      ! jw
      r1 = (min_dist(2)*min_dist(3))**power
      r2 = (min_dist(1)*min_dist(3))**power
      r3 = (min_dist(1)*min_dist(2))**power

      var_new1(i) = ((r1*var_old1(n1)+r2*var_old1(n2)+r3*var_old1(n3)) / (r1 + r2 + r3))
      var_new2(i) = ((r1*var_old2(n1)+r2*var_old2(n2)+r3*var_old2(n3)) / (r1 + r2 + r3))
      var_new3(i) = ((r1*var_old3(n1)+r2*var_old3(n2)+r3*var_old3(n3)) / (r1 + r2 + r3))
      
      deallocate(dist)
! jw
! jw
   end do   ! jw
   !$omp end parallel do
end subroutine IDW_5
