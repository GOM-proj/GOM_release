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
	! End of local variables ==================================================!
	
   ! pre calculate the density ratio between air and water:
   rho_ratio = rho_a/rho_o
	
	! define minimum and maximum drag coefficient for each formula
	if(wind_formular == 1) then
		! Garrat's (1977) equation
		! originally 4 [m/s] < wind speed < 21 [m/s]
	   Cda_min=0.001_dp*(0.75+0.067*4.0)	! minimum wind speed =  4 m/s
	   Cda_max=0.001_dp*(0.75+0.067*50.0)	! maximum wind speed = 50 m/s		
	else if(wind_formular == 2) then
		! Smith's (1980) equation
		! originally 6 [m/s] < wind speed < 22 [m/s]
	   Cda_min=0.001_dp*(0.61+0.063*6.0)	! minimum wind speed =  6 m/s
	   Cda_max=0.001_dp*(0.61+0.063*50.0)	! maximum wind speed = 50 m/s
	else if(wind_formular == 3) then
		! Wu's (1982) equation
		! originally 1 [m/s] < wind speed
	   Cda_min=0.001_dp*(0.80+0.065*1.0)	! minimum wind speed =  1 m/s
	   Cda_max=0.001_dp*(0.80+0.065*50.0)	! maximum wind speed = 50 m/s
	else
		! use default as Smith's equation
		! Smith's (1980) equation
		! originally 6 [m/s] < wind speed < 22 [m/s]
	   Cda_min=0.001_dp*(0.61+0.063*6.0)	! minimum wind speed =  6 m/s
	   Cda_max=0.001_dp*(0.61+0.063*50.0)	! maximum wind speed = 50 m/s
	end if
	

	! Calculate wind stress at each face (side) ===============================!
   if(wind_flag == 1 .and. heat_model_flag == 0) then
   	!$omp parallel do private(j,wind_speed,wind_normal,wind_tangent,Cda)
   	do j=1,maxface
   		! convert wind speed [m/s] in xy (true east and true north) coordinate to 
   		! face normal and face tangent wind speed [m/s]
         wind_speed 	 =  dsqrt(wind_u_at_face(j)**2 + wind_v_at_face(j)**2)
         wind_normal  =  wind_u_at_face(j)*cos_theta(j)   &
             &        +  wind_v_at_face(j)*sin_theta(j)
         wind_tangent = -wind_u_at_face(j)*sin_theta(j)   &
             &        +  wind_v_at_face(j)*cos_theta(j)
         
         ! choose wind formula
         if(wind_formular == 1) then
         	! Garrat's (1977) equation (also used in CH3D)
				Cda = 0.001_dp*(0.75_dp + 0.067_dp*wind_speed)				
			else if(wind_formular == 2) then
				! Smith's (1980) equation
         	Cda = 0.001_dp*(0.61_dp + 0.063_dp*wind_speed)
         else if(wind_formular == 3) then
         	! Wu's (1982) equation
         	Cda = 0.001_dp*(0.80_dp + 0.065_dp*wind_speed)
         end if
         
         Cda = min(max(Cda,Cda_min),Cda_max) ! (Cda_min <= Cda < Cda_max)
         Gamma_T(j) = Cda*rho_ratio*wind_speed
         wind_stress_normal(j) = Gamma_T(j)*wind_normal*spinup_function_wind 	! Cd*rho*u*|u|
         wind_stress_tangnt(j) = Gamma_T(j)*wind_tangent*spinup_function_wind	! Cd*rho*v*|u|         
      end do
      !$omp end parallel do
   end if
	
	! jw, I have to update this part using Gamma_T
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
				! sign and scale difference between stresses tauxz and wind_stress_normal
				! original
            ! wind_stress_normal(j) = - wind_stress_normal(j)/rho_o*spinup_function_wind
            ! wind_stress_tangnt(j) = - wind_stress_tangnt(j)/rho_o*spinup_function_wind
            ! jw, correction
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
         
         ! choose wind formula
         if(wind_formular == 1) then
         	! Garrat's (1977) equation (also used in CH3D)
				Cda = 0.001_dp*(0.75_dp + 0.067_dp*wind_speed)				
			else if(wind_formular == 2) then
				! Smith's (1980) equation
         	Cda = 0.001_dp*(0.61_dp + 0.063_dp*wind_speed)
         else if(wind_formular == 3) then
         	! Wu's (1982) equation
         	Cda = 0.001_dp*(0.80_dp + 0.065_dp*wind_speed)
         else
         	! use default as Smith's equation
         	! Smith's (1980) equation
         	Cda = 0.001_dp*(0.61_dp + 0.063_dp*wind_speed)
         end if
						
			! jw, Where this comes from?
			! Maybe from 0.002 < Cd < 0.003 from high wind speeds between 20 and 77 m/s
			Cda = MIN(0.003,Cda)

         Gamma_T(j) = Cda*rho_ratio*wind_speed
         wind_stress_normal(j) = Gamma_T(j)*wind_normal*spinup_function_wind
         wind_stress_tangnt(j) = Gamma_T(j)*wind_tangent*spinup_function_wind
      end do
      !$omp end parallel do
   end if
end subroutine calculate_wind_stress
