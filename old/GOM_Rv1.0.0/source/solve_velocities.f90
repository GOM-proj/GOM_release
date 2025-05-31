!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This is the second part of the solve momentum equation which solves face normal and tangential velocities (horizontal velocities), Eq(44)
!! and solve vertical velocities at element's center (top and bottom at each level) from the freesurface Eq(47)
!!  
subroutine solve_velocities
   use mod_global_variables
   use mod_file_definition
   implicit none
   
   integer :: i, j, k, l, kk
   integer :: j_temp, n1, n2, ibnd
   integer :: t_layer, b_layer
   real(dp):: vnorm
   real(dp):: sum1
   real(dp):: rtemp0, rtemp1, rtemp2
	! jw

	! jw
	! jw
 	!$omp parallel
 	!$omp do private(j,k,kk,n1,n2,rtemp0,rtemp1,rtemp2,ibnd,vnorm)
   do j = 1, maxface
      n1 = nodenum_at_face(1,j)
      n2 = nodenum_at_face(2,j)
		
      if(top_layer_at_face(j) == 0) then	
      	! jw
         do k = 1, maxlayer
            un_face_new(k,j) = 0.0_dp
            vn_face_new(k,j) = 0.0_dp
         end do
      else 
      	! jw
      	! jw
      	! jw
      	! jw
      	! jw
      	! jw
    	   ! jw
         do k = 1, top_layer_at_face(j) - bottom_layer_at_face(j) + 1
         	! jw
            kk = top_layer_at_face(j) + 1 - k 
            
            ! jw
    	      ! jw
    	      if(boundary_type_of_face(j) == -1) then
    	      	! jw
    	      	un_face_new(kk,j) = 0.0_dp    	      	
    	      else if(boundary_type_of_face(j) == 0) then 
    	      	! jw
    	      	! jw
    	      	! jw
    	      	rtemp0 = theta*gravity*dt/delta_j(j) &
               &	* ( eta_cell_new(adj_cellnum_at_face(2,j)) - eta_cell_new(adj_cellnum_at_face(1,j)) ) * AinvDeltaZ1(k,j)
               
               ! jw
               ! jw
               
               ! jw
               rtemp1 = rtemp0 * cos_theta2(j) ! jw
               rtemp2 = rtemp0 * sin_theta2(j) ! jw
               
               ! jw
               ! jw
               un_face_new(kk,j) = 	AinvG1(k,j) - rtemp1               
            else if(boundary_type_of_face(j) > 0) then
		 			! jw
		 			! jw
		 			! jw
		 			rtemp0 = theta*gravity*dt/delta_j(j) &
               &	* ( eta_at_ob_new(ob_element_flag(adj_cellnum_at_face(1,j))) - eta_cell_new(adj_cellnum_at_face(1,j)) )	&
               &	* AinvDeltaZ1(k,j)
               
               un_face_new(kk,j) = AinvG1(k,j) - rtemp0
            end if
            
		 		! jw
		 		! jw
		 		! jw
		 		! jw
            if(wetdry_node(n1) == 0 .and. wetdry_node(n2) == 0) then
            	! jw
               ! jw
               ! jw
               ! jw
               ! jw
               ! jw
               
               ! jw

					! jw
               rtemp1 = theta*gravity*dt/face_length(j) &
               &		* ( eta_node(n2) - eta_node(n1) )	&
               &		* AinvDeltaZ2(k,j)
               
               ! jw
               ! jw
               ! jw
	         	if(boundary_type_of_face(j) == 0) then
	         		! jw
	         		! jw
	         		! jw
	         		vn_face_new(kk,j) = AinvG2(k,j) - (rtemp1 + rtemp2)
		         else if(boundary_type_of_face(j) > 0) then
		         	! jw
		         	! jw
		         	vn_face_new(kk,j) = AinvG2(k,j) - rtemp1
		         else if(boundary_type_of_face(j) == -1) then
		         	! jw
		         	! jw
		         	vn_face_new(kk,j) = AinvG2(k,j) - rtemp1
		         end if
            end if

				! jw
            if(initial_wetdry_node(n1) == 1 .or. initial_wetdry_node(n2) == 1) then
               vn_face_new(kk,j) = dmax1(-5.0_dp, dmin1(vn_face_new(kk,j), 5.0_dp))
            end if
         end do ! jw


         ! jw
         if(isflowside3(j) > 0) then
         	ibnd = isflowside3(j)
         	do k = bottom_layer_at_face(j), top_layer_at_face(j)
         		un_face_new(k,j) = Qu_boundary(ibnd)
         	end do
         end if
         
         ! jw
			if(isflowside4(j) > 0) then
				ibnd = isflowside4(j)
				if(WR_layer(ibnd) == 999) then
					! jw
					k = top_layer_at_face(j)
					un_face_new(k,j) = WRu_boundary(ibnd)
				else if(WR_layer(ibnd) == 0) then
					! jw
					k = bottom_layer_at_face(j)
					! jw
					un_face_new(k,j) = WRu_boundary(ibnd)
				else
					k = WR_layer(ibnd)
					un_face_new(k,j) = WRu_boundary(ibnd)
				end if
			end if


			! jw
         if(isflowside2(j) > 0) then
         	ibnd  = isflowside2(j)
            vnorm = u_boundary(ibnd)*cos_theta(j)   &
              &   + v_boundary(ibnd)*sin_theta(j)

				write(*,*) ' in calculate velo at face'
				write(*,*) u_boundary(ibnd), v_boundary(ibnd)
				write(*,*) 'press enter to continue'

            do k = bottom_layer_at_face(j), top_layer_at_face(j)
               vn_face_new(k,j) = 0.0_dp
               un_face_new(k,j) = vnorm                                       &
                       &            + dsqrt(gravity/h_face(j))   &
                       &            * eta_cell_new(adj_cellnum_at_face(j,1))
            end do
         end if
         
                  
			! jw
			! jw
         do k = 1, bottom_layer_at_face(j)-1
            un_face_new(k,j) = 0.0_dp  ! jw
            vn_face_new(k,j) = 0.0_dp 
         end do

			! jw
         do k = top_layer_at_face(j)+1, maxlayer
            un_face_new(k,j) = un_face_new(top_layer_at_face(j),j)
            vn_face_new(k,j) = vn_face_new(top_layer_at_face(j),j)
         end do
      end if ! jw
   end do ! jw
 	!$omp end do
	
	! jw
	! jw
	! jw
 	!$omp do private(i,k,l,sum1,j_temp, t_layer, b_layer)
   do i = 1, maxele
   	t_layer = top_layer_at_element(i)
   	b_layer = bottom_layer_at_element(i)
   	
      if(t_layer == 0) then
         do k = 0, maxlayer
            wn_cell_new(k,i) = 0.0_dp
         end do
      else 
      	! jw
      	! jw
      	! jw
! jw
! jw
         	! jw
         	! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
			
			! jw
			wn_cell_new(b_layer-1,i) = 0.0_dp ! jw
			
         do k = b_layer, t_layer
            sum1 = 0.0_dp
            do l = 1, tri_or_quad(i)
               j_temp = facenum_at_cell(l,i)
               if(k >= bottom_layer_at_face(j_temp) .and. k <= top_layer_at_face(j_temp)) then
                  sum1 = sum1 + sign_in_outflow(l,i) * face_length(j_temp) * dz_face(k,j_temp) * un_face_new(k,j_temp)
               end if
            end do
            wn_cell_new(k,i) = wn_cell_new(k-1,i) - sum1/area(i)
         end do
			
			! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
			
			! jw
         do k = 0, b_layer-2
            wn_cell_new(k,i) = 0.0_dp
         end do

			! jw
         do k = t_layer+1, maxlayer
            wn_cell_new(k,i) = wn_cell_new(t_layer,i)
         end do
      end if
   end do
 	!$omp end do
 	!$omp end parallel

	! jw
	if(dia_face_velocity == 1) then
		write(pw_dia_face_velocity_uv,*) 'it = ', it, ', elapsed_time = ', elapsed_time		
		write(pw_dia_face_velocity_uv,*) 'un_face_new(1,j), vn_face_new(1,j)'
		do j=1,maxface
			write(pw_dia_face_velocity_uv,'(A3, I5, 2E15.5)') 'j=', j, un_face_new(1,j), vn_face_new(1,j)
		end do
	end if

	if(dia_face_velocity == 1) then
		write(pw_dia_face_velocity_w,*) 'it = ', it, ', elapsed_time = ', elapsed_time		
		write(pw_dia_face_velocity_w,*) 'wn_cell_new(1,i)'
		do i=1,maxele
			write(pw_dia_face_velocity_w,'(A3, I5, E15.5)') 'i=', i, wn_cell_new(1,i)
		end do
	end if
	
	! jw
   
end subroutine solve_velocities
