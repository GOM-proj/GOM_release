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
	! jw
	
	! jw
	! jw
	! jw
	ana_wind_stress_x = 0.1_dp
	ana_wind_stress_y = 0.0_dp
	
	! jw
	! jw
	! jw
	rho_water = 1000.0_dp
	
	rho_water2 = 1.0/rho_water ! jw
	
	! jw
	! jw
	! jw
	!$omp parallel do private(j)
	do j=1,maxface
		wind_stress_normal(j) =  (ana_wind_stress_x * cos_theta(j) + ana_wind_stress_y * sin_theta(j))*rho_water2
		wind_stress_tangnt(j) = (-ana_wind_stress_x * sin_theta(j) + ana_wind_stress_y * cos_theta(j))*rho_water2
	end do
	!$omp end parallel do
		
end subroutine ana_windp