!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! Set analytical wind stress [N/m2] at every sides
!!
subroutine ana_windp
	use mod_global_variables
	implicit none
	integer :: j
	real(dp):: ana_wind_stress_x, ana_wind_stress_y
	real(dp):: rho_water, rho_water2
	! End of local variables ==================================================!
	
	! Set analytical wind stress [N/m2] at every sides
	! Provide wind stress in true xy-coordinate
	! 1 dyne/cm2 = 0.1 N/m2
	ana_wind_stress_x = 0.1_dp
	ana_wind_stress_y = 0.0_dp
	
	! provide the reference water density:
	! Use water reference water density rho_o or the ideal water denisty 1000.0
	! rho_water = rho_o ! or
	rho_water = 1000.0_dp
	
	rho_water2 = 1.0/rho_water ! I am using this to avoid using division.
	
	! Do not change below this line ===========================================!
	! Since we have to set wind stress at each face (not each element), 
	! we have to calculate [wind_stress_normal] and [wind_stress_tangnt] at each face.
	!$omp parallel do private(j)
	do j=1,maxface
		wind_stress_normal(j) =  (ana_wind_stress_x * cos_theta(j) + ana_wind_stress_y * sin_theta(j))*rho_water2
		wind_stress_tangnt(j) = (-ana_wind_stress_x * sin_theta(j) + ana_wind_stress_y * cos_theta(j))*rho_water2
	end do
	!$omp end parallel do
		
end subroutine ana_windp