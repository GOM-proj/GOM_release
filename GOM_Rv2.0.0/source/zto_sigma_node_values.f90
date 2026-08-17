!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine zto_sigma_node_values(u_sigma, v_sigma, w_sigma)
	use mod_global_variables
	implicit none
	
	real(dp),dimension(0:maxlayer,maxnod),intent(inout) :: u_sigma, v_sigma, w_sigma

	integer :: i, k, k1, k2
	integer :: num_vertical_layer
	real(dp):: total_water_depth, dz2
	real(dp):: u1,u2,u3
	real(dp):: v1_u, v1_v, v1_w
	real(dp):: v3_u, v3_v, v3_w
	
	! jw
	
	! jw
	!$omp parallel do private(i,k,k1,k2, &
	!$omp &						  total_water_depth,dz2,num_vertical_layer, &
	!$omp &						  u1,u2,u3, v1_u,v1_v,v1_w, v3_u,v3_v,v3_w)
	do i=1,maxnod
		! jw
		if(top_layer_at_node(i) ==0 ) then
			do k=0,maxlayer
				u_sigma(k,i) = 0.0_dp
				v_sigma(k,i) = 0.0_dp
				w_sigma(k,i) = 0.0_dp
			end do
			continue ! jw
		end if
		
		! jw
		! jw
      total_water_depth = h_node(i) + eta_node(i)
      
      ! jw
      ! jw
		dz2 = total_water_depth/maxlayer ! jw
		
		! jw
		num_vertical_layer = top_layer_at_node(i) - bottom_layer_at_node(i) + 1 ! jw
		
		! jw
		do k2=0,maxlayer ! jw
			u2 = -h_node(i) + dz2*k2 ! jw
			
			do k1=bottom_layer_at_node(i),top_layer_at_node(i) ! jw
				if(k1 == bottom_layer_at_node(i)) then
					u1 = -h_node(i) ! jw
					u3 = -h_node(i) + dz_node(k1,i) ! jw
					
					if(u2 >= u1 .and. u2 <= u3) then
						v1_u = u_node(k1-1,i)
						v3_u = u_node(k1,i)			

						v1_v = v_node(k1-1,i)
						v3_v = v_node(k1,i)

						v1_w = w_node(k1-1,i)
						v3_w = w_node(k1,i)
						exit ! jw
					end if
				else if(k1 == top_layer_at_node(i)) then
					u1 = -(MSL-(z_level(k1-1))) ! jw
					u3 = -(MSL-(z_level(k1-1) + dz_node(k1,i))) ! jw
					
					if(u2 >= u1 .and. u2 <= u3) then
						v1_u = u_node(k1-1,i)
						v3_u = u_node(k1,i)
						
						v1_v = v_node(k1-1,i)
						v3_v = v_node(k1,i)

						v1_w = w_node(k1-1,i)
						v3_w = w_node(k1,i)
						exit ! jw
					end if
				else
					u1 = -(MSL-z_level(k1-1)) ! jw
					u3 = -(MSL-(z_level(k1-1) + dz_node(k1,i))) ! jw
					if(u2 >= u1 .and. u2 <= u3) then
						v1_u = u_node(k1-1,i)
						v3_u = u_node(k1,i)
						
						v1_v = v_node(k1-1,i)
						v3_v = v_node(k1,i)

						v1_w = w_node(k1-1,i)
						v3_w = w_node(k1,i)
						exit ! jw
					end if
				end if				
			end do ! jw
			! jw
			! jw
			u_sigma(k2,i) = (v3_u-v1_u)*(u2-u1)/(u3-u1) + v1_u
			v_sigma(k2,i) = (v3_v-v1_v)*(u2-u1)/(u3-u1) + v1_v
			w_sigma(k2,i) = (v3_w-v1_w)*(u2-u1)/(u3-u1) + v1_w
		end do ! jw
		! jw
	end do ! jw
	!$omp end parallel do
end subroutine zto_sigma_node_values
