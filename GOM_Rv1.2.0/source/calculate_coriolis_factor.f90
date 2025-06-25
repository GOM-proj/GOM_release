!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Calculation of coriolis factor according to Coriolis options
!! 
subroutine calculate_coriolis_factor
   use mod_global_variables
   use mod_file_definition
   implicit none
   
   integer :: j
   real(dp):: r, omega, f0, beta, y
   real(dp):: x2, y2
   ! End of local variables ==================================================!
   
   omega = 7.292E-5 	! earth's angular velocity [rad/sec]: (2pi/[23h 56min 4.09sec])

	!  initialize coriolis parameter for betaplane approximation
   if(Coriolis_option == 0) then			! No Coriolis
      coriolis_factor = 0.0_dp
   else if(Coriolis_option == 1) then 	! f-plane approximation
  		coriolis_factor = 2*omega*sin(lat_mid*deg2rad)
   else if(Coriolis_option == 2) then	! beta-plane approximation
   	if(coordinate_system == 2) then	! beta-plane approximation is only allowed when lon/lat is provided
			! prepare for the beta-plane approximation
			! Reference: p14, Introduction to Physical Oceanography (George Mellor)
			r = 6371*1000.0_dp 	! earth's radius
			f0 = 2*omega*sin(lat_mid*deg2rad)
			beta = 2*omega*cos(lat_mid*deg2rad)/r 	! note, this is the original equation
			! beta = 2*omega*cos(lat_mid*deg2rad) 	! note, I can cancle "/r" term, but I will keep the original equation
			
			! if the given node.inp is in UTM, then I need to find the mean lon/lat
			! Note:
			! 		I can do this step in read_node_inp.f90, 
			! 		but it is better to be here since this calculation is not required if coriolis_option == 1.
			!$omp parallel do private(j,x2,y2,y)
	   	do j=1,maxface
	   		! calculate face center position's latitude
				! conversion option in coordinate_conversion.f90:
				! 		1: lon/lat -> UTM
				! 		2: UTM -> lon/lat
	   		call coordinate_conversion(x_face(j),y_face(j),utm_projection_zone,2, x2,y2) ! (x2,y2) = (lon,lat) for this point
	   		y = r*(y2 - lat_mid)	! this is the original equation
	   		! y = (y2 - lat_mid) ! I can cancle "r" term, but I will keep the original equation
	   		coriolis_factor(j) = f0 + beta * y
	   	end do
	   	!$omp end parallel do
	   else
	   	write(pw_run_log,*) 'Error: beta-plane approximation is only allowed when lon/lat is provided in node.inp'
	   	write(pw_run_log,*) 'Stop'
	   	write(*,*) 'Error: beta-plane approximation is only allowed when lon/lat is provided in node.inp'
	   	write(*,*) 'Stop'
	   	stop
	   end if
   end if
end subroutine calculate_coriolis_factor
