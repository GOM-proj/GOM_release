!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! jw, not yet corrected.
!! 
!! This subroutine calculates hypothetical wind field in each mesh node 
!! over the entire nodes using holland's model
!! 
!! ===========================================================================!
subroutine holland_storm_surge        
   use mod_global_variables
   use mod_function_library, only : coriolis_from_lat_deg
   
   implicit none
   integer :: i
   real(dp):: ratio  
   real(dp),allocatable,dimension(:) :: distance_from_node, distance_from_element  
   real(dp),parameter :: inward_rotation =18.0 
   real(dp)::    x, y 
   real(dp):: hurric_constant_1, gamma, beta_plus_theta, theta_hur 
   real(dp):: cyclostrophic_velocity, geostrophic_velocity 
   real(dp):: gradient_velocity, angle 
   real(dp):: delta_x_hurric, delta_y_hurric 
   real(dp):: max_cyclostrophic_wind_speed 
   real(dp):: hurric_translation_velocity 
   real(dp):: max_surface_wind_speed
   real(dp),parameter :: air_density = 1.293d0  ! jw
   real(dp),parameter :: gradient_to_surface_wind = 0.95d0 
   real(dp),parameter :: air_p_freestream = 1013.0d0 ! jw

   integer :: temp_ob_element, n1, n2

   real(dp):: holland_temp1, holland_temp2
   real(dp),allocatable,dimension(:) :: 		&
   &			wind_speed_x_holland_hurricane,  &
	&			wind_speed_y_holland_hurricane,  &
	&			air_pressure_holland_hurricane
	! jw
	
   do while (current_jday_1900 > hurric_julian_day_1900_2   &
     &                       .and.                          &
     &                   .not. is_iostat_end(650)                   ) 
        
      hurric_year_1   = hurric_year_2 
      hurric_month_1  = hurric_month_2 
      hurric_day_1    = hurric_day_2 
      hurric_hour_1   = hurric_hour_2 
      hurric_minute_1 = hurric_minute_2 
        
      hurric_julian_day_1      = hurric_julian_day_2 
      hurric_julian_day_1900_1 = hurric_julian_day_1900_2 
    
      hurric_x_1              = hurric_x_2 
      hurric_y_1              = hurric_y_2 
      hurric_latitude_1       = hurric_latitude_2 
      hurric_delta_pressure_1 = hurric_delta_pressure_2 
      hurric_mwr_1            = hurric_mwr_2 

      read (650,*) hurric_year_2, hurric_month_2,                   & 
       &           hurric_day_2, hurric_hour_2, hurric_minute_2, 	&
       &           hurric_x_2, hurric_y_2, hurric_latitude_2,       &
       &           hurric_delta_pressure_2, hurric_mwr_2 
        
      call julian(hurric_year_2, hurric_month_2, hurric_day_2,    &
       &          hurric_hour_2, hurric_minute_2, 0,            &
       &          hurric_julian_day_2, hurric_julian_day_1900_2) 

   enddo  ! jw
!*************************************************************************************************
!*****          calculate the ratio to multiply the variables by based upon                  *****
!*****                       the time between the two readings.                              *****
!*************************************************************************************************
   allocate(wind_speed_x_holland_hurricane(maxnod))
   allocate(wind_speed_y_holland_hurricane(maxnod))
   allocate(air_pressure_holland_hurricane(maxnod))
   allocate(distance_from_node(maxnod))
   allocate(distance_from_element(maxele))

   if(  julian_day >= holland_start_jday   &
     &            .and.                          &
     &  julian_day <= holland_end_jday     &
     &            .and.                          &
     &            .not. is_iostat_end(650) ) then     	
      ratio = (current_jday_1900        - hurric_julian_day_1900_1)    &
        &   / (hurric_julian_day_1900_2 - hurric_julian_day_1900_1) 
!*************************************************************************************************
!*****                      get the values at the current timestep                           *****
!*************************************************************************************************
      hurric_x = hurric_x_1 + (hurric_x_2 - hurric_x_1) * ratio 
      hurric_y = hurric_y_1 + (hurric_y_2 - hurric_y_1) * ratio 
  
      hurric_delta_pressure =  hurric_delta_pressure_1   &
               &            + (hurric_delta_pressure_2 - hurric_delta_pressure_1) * ratio 
  
      hurric_latitude =  hurric_latitude_1    &
               &      + (hurric_latitude_2 - hurric_latitude_1) * ratio 
      
      hurric_mwr = hurric_mwr_1 + (hurric_mwr_2 - hurric_mwr_1) * ratio 
!*************************************************************************************************
!*****                                    convert units                                      *****
!*************************************************************************************************
      hurric_mwr  = hurric_mwr  * 1.852000e+03 
      hurric_delta_pressure = hurric_delta_pressure * 3.386379e+03 
!*************************************************************************************************
!*****                               calculate some constants                                *****
!*************************************************************************************************
      hurric_coriolis = coriolis_from_lat_deg(hurric_latitude) 
      delta_x_hurric = hurric_x_2-hurric_x_1 
      delta_y_hurric = hurric_y_2-hurric_y_1 
      if(delta_x_hurric /= 0.0) then 
         hurric_beta = atan(delta_y_hurric/delta_x_hurric) ! jw
      else 
         if(delta_y_hurric >= 0.0) then 
            hurric_beta =  pi/2.0d0 
         else 
            hurric_beta = -pi/2.0d0 
         endif 
      endif 

      if(delta_x_hurric < 0.0 ) then 
         hurric_beta = pi + hurric_beta 
      endif 

      hurric_translation_velocity = sqrt( (hurric_x_2-hurric_x_1)**2      & 
                  &               +       (hurric_y_2-hurric_y_1)**2 )    &
                  &               / (hurric_julian_day_1900_2-hurric_julian_day_1900_1)   &
                  &               /  86400.0d0 
!*************************************************************************************************
!*****                     calculate the maximum wind speed                                  *****
!*************************************************************************************************
      max_cyclostrophic_wind_speed = sqrt(-hurric_delta_pressure/air_density/exp(1.0))  !m/s 
      max_surface_wind_speed = gradient_to_surface_wind    &
                    &        * max_cyclostrophic_wind_speed + 0.5*hurric_translation_velocity 
!*************************************************************************************************
!*****                  calculate the wind speed in the x- and y-directions                  ***** 
!*************************************************************************************************
      do i = 1, maxnod
         x = x_node(i)
         y = y_node(i)

         distance_from_node(i) = sqrt( (hurric_x-x)**2 + (hurric_y-y)**2 ) 
         if(distance_from_node(i) == 0.0) then
            write(*,*) ' distance_from_node(i) = zero at ', 'i = ', i
            write(*,*) ' check distance '
            write(*,*) 'press enter to continue'		  
         endif         
         if(x >= 0.9e+15 .and. y >= 0.9e+15) then ! jw
            gradient_velocity = 0.0d0 
            distance_from_node(i) = 1.0d0 
            x = 0.9e+19 
            y = 0.9e+19 
         else 
            delta_x_hurric = x - hurric_x 
            delta_y_hurric = y - hurric_y 
            if(delta_x_hurric == 0.0)then
               delta_x_hurric = 0.00001d0
            endif
            beta_plus_theta = atan(delta_y_hurric/delta_x_hurric) ! jw
            if(delta_x_hurric < 0.0 ) then 
               beta_plus_theta = pi + beta_plus_theta 
            endif 
            theta_hur = beta_plus_theta - hurric_beta 
            hurric_constant_1 = (-hurric_delta_pressure*hurric_mwr)   &
                  &           *   exp(-hurric_mwr/distance_from_node(i))      &
                  &           / ( distance_from_node(i)*air_density) 
               
            cyclostrophic_velocity = sqrt(hurric_constant_1)      !unit = m/s
              
            if(cyclostrophic_velocity < 0.001d0) then 
               gradient_velocity = 0.0d0 
            else 
               geostrophic_velocity =  hurric_constant_1   &
                          &         / (distance_from_node(i)*hurric_coriolis) 
               gamma = 0.5d0*(   hurric_translation_velocity                  &
                 &           * sin(theta_hur)/cyclostrophic_velocity          &
                 &           + cyclostrophic_velocity/geostrophic_velocity)   
              
               gradient_velocity = cyclostrophic_velocity*(sqrt(gamma*gamma + 1) - gamma)   &
                       &         * gradient_to_surface_wind 
            endif 
         endif 
              
         angle = pi*0.5d0 + theta_hur + hurric_beta + inward_rotation*pi/180.0d0 
            
         wind_speed_x_holland_hurricane(i) = gradient_velocity*cos(angle)   ! jw
         wind_speed_y_holland_hurricane(i) = gradient_velocity*sin(angle)   ! jw

      enddo  ! jw
!*************************************************************************************************
!*****            calculate distance fron element center to hurricane center                 *****
!*************************************************************************************************
      do i = 1, maxele
         x = x_cell(i)
         y = y_cell(i)

         distance_from_element(i) = sqrt( (hurric_x-x)**2 + (hurric_y-y)**2 ) 
         if(distance_from_element(i) == 0.0) then
            write(*,*) ' distance_from_element(i) = zero at ', 'i = ', i
            write(*,*) ' check distance '
            write(*,*) 'press enter to continue'		  
         endif         
            
         if(x >= 0.9e+15 .and. y >= 0.9e+15) then ! jw
            gradient_velocity = 0.0d0 
            distance_from_element(i) = 1.0d0 
            x = 0.9e+19 
            y = 0.9e+19 
         endif 
      enddo  
!*************************************************************************************************
!*****    calculate atmospheric pressure in all the elements including open boundaries       *****
!*************************************************************************************************
      do i = 1, maxnod
         if(x >= 0.9e+15 .and. y >= 0.9e+15) then ! jw
            air_pressure_holland_hurricane(i) = 0.0d0 
         else 
            air_pressure_holland_hurricane(i) =  air_p_freestream + hurric_delta_pressure   &
               &     * (1.0d0 - exp(-hurric_mwr/distance_from_node(i)))/100.0d0
                                                                    ! jw
         endif 
      enddo 
!*************************************************************************************************
!*****          calculate elevation due to pressure deficit along open boundaries            ***** 
!*************************************************************************************************		
		if(holland_flag == 1) then ! jw
      	do i = 1, num_ob_cell
         	temp_ob_element = ob_cell_id(i)
         	eta_from_Holland_at_ob(i) =                                                 &
            &      -  hurric_delta_pressure                                            &
            &      * ( 1.0d0 - exp(-hurric_mwr/distance_from_element(temp_ob_element)) &
            &        )                                                                 &
            &      / (rho_o)/gravity   ! jw
         enddo
      endif      
   else
      do i = 1, maxnod
            wind_speed_x_holland_hurricane(i) = 0.0d0 
            wind_speed_y_holland_hurricane(i) = 0.0d0 
! jw
            air_pressure_holland_hurricane(i) = 0.0d0
      enddo 

		do i = 1, num_ob_cell
			temp_ob_element = ob_cell_id(i)
			eta_from_Holland_at_ob(i) = 0.0d0
		end do
	end if   ! jw


!*************************************************************************************************
!*****  wind velocity for storm surge model(holland model wind)
!*************************************************************************************************

   do i = 1, maxnod
      wind_u1(i) = wind_u2(i)
      wind_v1(i) = wind_v2(i)
      wind_u2(i) = wind_speed_x_holland_hurricane(i)
      wind_v2(i) = wind_speed_y_holland_hurricane(i)
   enddo
   do i = 1, maxnod
      wind_u0(i) = wind_u2(i)
      wind_v0(i) = wind_v2(i)
   enddo !i

   do i = 1, maxface
      n1 = nodenum_at_face(1,i)
      n2 = nodenum_at_face(2,i)
      holland_temp1 = (wind_u2(n1) + wind_u2(n2))/2.0
      holland_temp2 = (wind_v2(n1) + wind_v2(n2))/2.0
      wind_u_at_face(i) = holland_temp1
      wind_v_at_face(i) = holland_temp2
   enddo

!***************************************************************************************************
!*****                include storm surge model for atmospheric pressure term                  *****
!***************************************************************************************************
   do i = 1, maxnod 
      air_p1(i) = air_pressure_holland_hurricane(i)*100.0d0
   enddo

   deallocate(distance_from_node)
   deallocate(distance_from_element)

   deallocate(wind_speed_x_holland_hurricane,wind_speed_y_holland_hurricane)
   deallocate(air_pressure_holland_hurricane)

   return
!*************************************************************************************************
end subroutine holland_storm_surge
!*************************************************************************************************



