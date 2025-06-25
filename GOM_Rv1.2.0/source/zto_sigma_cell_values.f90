!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine zto_sigma_cell_values(salt_sigma, temp_sigma, rho_sigma)
	use mod_global_variables
	implicit none
	
	real(dp),dimension(maxlayer,maxele),intent(inout) :: salt_sigma, temp_sigma, rho_sigma

	integer :: i, k, k1, k2
	integer :: num_vertical_layer
	real(dp):: total_water_depth, dz2
	real(dp):: u1,u2,u3
	real(dp):: v1_salt, v1_temp, v1_rho
	real(dp):: v3_salt, v3_temp, v3_rho

	! End of local variables ==================================================!
	
	! interpolate cell values into sigma layer
	!$omp parallel do private(i,k,k1,k2, &
	!$omp &						  total_water_depth,dz2,num_vertical_layer, &
	!$omp &						  u1,u2,u3, v1_salt,v1_temp,v1_rho, v3_salt,v3_temp,v3_rho)	
	do i=1,maxele
		! if this is a dry element, put "0" value for variables, otherwise interpolate values.
		if(top_layer_at_element(i) ==0 ) then
			do k=1,maxlayer
				salt_sigma(k,i) = 0.0_dp
				temp_sigma(k,i) = 0.0_dp
				rho_sigma(k,i) = 0.0_dp
			end do
			continue ! skip rest part and go to the next element
		end if
		
		! rest part is for wet element...
		! calculate total water depth.
      total_water_depth = h_cell(i) + eta_cell(i)
      
      ! calculate dz in both sigma and z-grid systems
      ! dz in sigma system:
		dz2 = total_water_depth/maxlayer ! dz in sigma system
		
		! dz in z-grid system
		num_vertical_layer = top_layer_at_element(i) - bottom_layer_at_element(i) + 1 ! number of vertical layers at this element for z-grid system
		
		! variables should be written from bottom to top order since I created each brick from bottom to top.
		do k2=1,maxlayer ! sigma-layer
			u2 = -h_cell(i) + 0.5*dz2 + dz2*(k2-1) ! cell center position for each sigma layer
			
			do k1=bottom_layer_at_element(i),top_layer_at_element(i) ! z-layer system		
				
				if(k1 == bottom_layer_at_element(i)) then
					u1 = -h_cell(i)
					u3 = -h_cell(i) + dz_cell(k1,i)*0.5 ! cell center position for bottom z-layer
					
					if(u2 >= u1 .and. u2 <= u3) then
						v1_salt = salt_cell(k1,i)
						v3_salt = salt_cell(k1,i)					

						v1_temp = temp_cell(k1,i)
						v3_temp = temp_cell(k1,i)

						v1_rho = rho_cell(k1,i)
						v3_rho = rho_cell(k1,i)
						exit ! exit z-layer loop
					end if
				else if(k1 == top_layer_at_element(i)) then
					u1 = -(MSL-(z_level(k1-1) + dz_cell(k1,i)*0.5)) ! layer center is located at 0.5*dz above from the bottom level
					u3 = -(MSL-(z_level(k1-1) + dz_cell(k1,i))) ! location of free surface
					
					if(u2 >= u1 .and. u2 <= u3) then
						v1_salt = salt_cell(k1,i)
						v3_salt = salt_cell(k1,i)
						
						v1_temp = temp_cell(k1,i)
						v3_temp = temp_cell(k1,i)

						v1_rho = rho_cell(k1,i)
						v3_rho = rho_cell(k1,i)
						exit ! exit z-layer loop
					end if
				else
					u1 = -(MSL-(z_level(k1-1) - 0.5*dz_cell(k1-1,i))) ! lower z-layer center position
					u3 = -(MSL-(z_level(k1-1) + 0.5*dz_cell(k1,i))) ! current z-layer center position
					if(u2 >= u1 .and. u2 <= u3) then
						v1_salt = salt_cell(k1-1,i)
						v3_salt = salt_cell(k1,i)
						
						v1_temp = temp_cell(k1-1,i)
						v3_temp = temp_cell(k1,i)

						v1_rho = rho_cell(k1-1,i)
						v3_rho = rho_cell(k1,i)
						exit ! exit z-layer loop
					end if
				end if				
			end do ! z-layer loop
			! linear interpolation:
			! 		these are for v2 values (i.e., interpolated values in sigma layer)
			salt_sigma(k2,i) = (v3_salt-v1_salt)*(u2-u1)/(u3-u1) + v1_salt
			temp_sigma(k2,i) = (v3_temp-v1_temp)*(u2-u1)/(u3-u1) + v1_temp
			rho_sigma(k2,i) = (v3_rho-v1_rho)*(u2-u1)/(u3-u1) + v1_rho
		end do ! do k2=1,maxlayer ! sigma-layer
	end do ! do i=1,maxele
	!$omp end parallel do
end subroutine zto_sigma_cell_values
