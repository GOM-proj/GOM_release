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
   ! jw
   
   omega = 7.292E-5 	! jw

	! jw
   if(Coriolis_option == 0) then			! jw
      coriolis_factor = 0.0_dp
   else if(Coriolis_option == 1) then 	! jw
  		coriolis_factor = 2*omega*sin(lat_mid*deg2rad)
   else if(Coriolis_option == 2) then	! jw
   	if(coordinate_system == 2) then	! jw
			! jw
			! jw
			r = 6371*1000.0_dp 	! jw
			f0 = 2*omega*sin(lat_mid*deg2rad)
			beta = 2*omega*cos(lat_mid*deg2rad)/r 	! jw
			! jw
			
			! jw
			! jw
			! jw
			! jw
			!$omp parallel do private(j,x2,y2,y)
	   	do j=1,maxface
	   		! jw
				! jw
				! jw
				! jw
	   		call coordinate_conversion(x_face(j),y_face(j),utm_projection_zone,2, x2,y2) ! jw
	   		y = r*(y2 - lat_mid)	! jw
	   		! jw
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
