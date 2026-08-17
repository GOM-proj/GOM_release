!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! compute wind stress components
!! jw, check it again.
!! not yet correctly included heat_model: see tauxz & tauyz
subroutine calculate_wind_stress 
   use mod_global_variables
   implicit none
   
   integer :: j, n1, n2 
   real(dp) :: wind_speed, wind_normal, wind_tangent   
   real(dp) :: Cda, Cda_min, Cda_max, rho_ratio
   real(dp),dimension(maxnod) :: tauxz, tauyz
	! jw
	
   ! jw
   rho_ratio = rho_a/rho_o
	
	! jw
	if(wind_formular == 1) then
		! jw
		! jw
	   Cda_min=0.001_dp*(0.75+0.067*4.0)	! jw
	   Cda_max=0.001_dp*(0.75+0.067*50.0)	! jw
	else if(wind_formular == 2) then
		! jw
		! jw
	   Cda_min=0.001_dp*(0.61+0.063*6.0)	! jw
	   Cda_max=0.001_dp*(0.61+0.063*50.0)	! jw
	else if(wind_formular == 3) then
		! jw
		! jw
	   Cda_min=0.001_dp*(0.80+0.065*1.0)	! jw
	   Cda_max=0.001_dp*(0.80+0.065*50.0)	! jw
	else
		! jw
		! jw
		! jw
	   Cda_min=0.001_dp*(0.61+0.063*6.0)	! jw
	   Cda_max=0.001_dp*(0.61+0.063*50.0)	! jw
	end if
	

	! jw
   if(wind_flag == 1 .and. heat_model_flag == 0) then
   	!$omp parallel do private(j,wind_speed,wind_normal,wind_tangent,Cda)
   	do j=1,maxface
   		! jw
   		! jw
         wind_speed 	 =  dsqrt(wind_u_at_face(j)**2 + wind_v_at_face(j)**2)
         wind_normal  =  wind_u_at_face(j)*cos_theta(j)   &
             &        +  wind_v_at_face(j)*sin_theta(j)
         wind_tangent = -wind_u_at_face(j)*sin_theta(j)   &
             &        +  wind_v_at_face(j)*cos_theta(j)
         
         ! jw
         if(wind_formular == 1) then
         	! jw
				Cda = 0.001_dp*(0.75_dp + 0.067_dp*wind_speed)				
			else if(wind_formular == 2) then
				! jw
         	Cda = 0.001_dp*(0.61_dp + 0.063_dp*wind_speed)
         else if(wind_formular == 3) then
         	! jw
         	Cda = 0.001_dp*(0.80_dp + 0.065_dp*wind_speed)
         end if
         
         Cda = min(max(Cda,Cda_min),Cda_max) ! jw
         Gamma_T(j) = Cda*rho_ratio*wind_speed
         wind_stress_normal(j) = Gamma_T(j)*wind_normal*spinup_function_wind 	! jw
         wind_stress_tangnt(j) = Gamma_T(j)*wind_tangent*spinup_function_wind	! jw
      end do
      !$omp end parallel do
   end if
	
	! jw
   if(wind_flag == 1 .and. heat_model_flag /= 0) then
   	!$omp parallel do private(j,n1,n2)
   	do j=1,maxface
         n1 = nodenum_at_face(1,j)
         n2 = nodenum_at_face(2,j)

         if(top_layer_at_node(n1) == 0 .or. top_layer_at_node(n2) == 0) then
            wind_stress_normal(j) = 0.0
            wind_stress_tangnt(j) = 0.0
         else
            wind_stress_normal(j) =  (tauxz(n1)+tauxz(n2))/2.*cos_theta(j)   &
                    &             +  (tauyz(n1)+tauyz(n2))/2.*sin_theta(j)
            wind_stress_tangnt(j) = -(tauxz(n1)+tauxz(n2))/2.*sin_theta(j)   &
                    &             +  (tauyz(n1)+tauyz(n2))/2.*cos_theta(j)
				! jw
				! jw
            ! jw
            ! jw
            ! jw
            wind_stress_normal(j) = - wind_stress_normal(j)*rho_ratio*spinup_function_wind
            wind_stress_tangnt(j) = - wind_stress_tangnt(j)*rho_ratio*spinup_function_wind
         end if
      end do
      !$omp end parallel do
   end if

   if(hurricane_flag == 1) then
   	!$omp parallel do private(j,wind_speed,wind_normal,wind_tangent,Cda)
   	do j=1,maxface
         wind_speed = dsqrt(wind_u_at_face(j)**2 + wind_v_at_face(j)**2)
         wind_normal  =  wind_u_at_face(j)*cos_theta(j)   &
             &        +  wind_v_at_face(j)*sin_theta(j)
         wind_tangent = -wind_u_at_face(j)*sin_theta(j)   &
             &        +  wind_v_at_face(j)*cos_theta(j)
         
         ! jw
         if(wind_formular == 1) then
         	! jw
				Cda = 0.001_dp*(0.75_dp + 0.067_dp*wind_speed)				
			else if(wind_formular == 2) then
				! jw
         	Cda = 0.001_dp*(0.61_dp + 0.063_dp*wind_speed)
         else if(wind_formular == 3) then
         	! jw
         	Cda = 0.001_dp*(0.80_dp + 0.065_dp*wind_speed)
         else
         	! jw
         	! jw
         	Cda = 0.001_dp*(0.61_dp + 0.063_dp*wind_speed)
         end if
						
			! jw
			! jw
			Cda = MIN(0.003,Cda)

         Gamma_T(j) = Cda*rho_ratio*wind_speed
         wind_stress_normal(j) = Gamma_T(j)*wind_normal*spinup_function_wind
         wind_stress_tangnt(j) = Gamma_T(j)*wind_tangent*spinup_function_wind
      end do
      !$omp end parallel do
   end if
end subroutine calculate_wind_stress
