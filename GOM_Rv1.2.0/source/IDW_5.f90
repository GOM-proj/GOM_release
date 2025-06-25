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
   ! End of local variables ==================================================!

   ! weighting factor: normally 1 or 2
   power = 1
	
	! interpolate with atmost 3 stations
	! num_stations = min(3,npoints)
	num_stations = 3

	! Inverse Distance Weighting (IDW) Interpolation with power = 1 ===========!
	!$omp parallel do private(i,k,min_dist,dist,xn,yn,dist_x,dist_y,ii,n1,n2,n3,r1,r2,r3)
   do i = 1, maxnod
      ! First, just set to a big reference number
      min_dist(1) = 1.0e15	! 1st closest distance
      min_dist(2) = 1.0e15 	! 2nd closest distance
      min_dist(3) = 1.0e15	! 3rd closest distance
      
      allocate(dist(hurricane_search_count(i)))
   	dist = 1.0e15
   	
   	! shifted new grid point
   	xn = xloc(i)
   	yn = yloc(i)
   	
   	! calculate distance between pre-defined neighbour nodes only
   	do k=1,hurricane_search_count(i)
   		dist_x = xn - x_node(hurricane_search_nodes(i,k))
   		dist_y = yn - y_node(hurricane_search_nodes(i,k))
   		dist(k) = sqrt(dist_x**2 + dist_y**2)
   	end do
   	
   	! 1st closest nodes
   	ii = minloc(dist,dim=1)
   	n1 = hurricane_search_nodes(i,ii)
   	min_dist(1) = dist(ii)
   	
   	! 2nd closest nodes
   	dist(ii) = 1.0e15
   	ii = minloc(dist,dim=1)
   	n2 = hurricane_search_nodes(i,ii)
   	min_dist(2) = dist(ii)

   	! 3rd closest nodes
   	dist(ii) = 1.0e15
   	ii = minloc(dist,dim=1)
   	n3 = hurricane_search_nodes(i,ii)
   	min_dist(3) = dist(ii)
   	
      ! Inverse Distance Weighting (IDW) Interpolation with power = 1 & 3 stations
      r1 = (min_dist(2)*min_dist(3))**power
      r2 = (min_dist(1)*min_dist(3))**power
      r3 = (min_dist(1)*min_dist(2))**power

      var_new1(i) = ((r1*var_old1(n1)+r2*var_old1(n2)+r3*var_old1(n3)) / (r1 + r2 + r3))
      var_new2(i) = ((r1*var_old2(n1)+r2*var_old2(n2)+r3*var_old2(n3)) / (r1 + r2 + r3))
      var_new3(i) = ((r1*var_old3(n1)+r2*var_old3(n2)+r3*var_old3(n3)) / (r1 + r2 + r3))
      
      deallocate(dist)
! 		write(*,*) i, n1, n2, n3, dist(1), dist(2), dist(3), var_old1(1), var_old2(1), var_new1(1), var_new2(1)
! 		stop 'I am here...'
   end do   ! do i = 1, maxnod
   !$omp end parallel do
end subroutine IDW_5
