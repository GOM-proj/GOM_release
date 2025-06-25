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
	! End of local variables ==================================================!

	! Calculate new horizontal velocities u & v at face at each layer =========!
	! solve the momentum equation for normal & tangent velocities , equation (13) in page 337
 	!$omp parallel
 	!$omp do private(j,k,kk,n1,n2,rtemp0,rtemp1,rtemp2,ibnd,vnorm)
   do j = 1, maxface
      n1 = nodenum_at_face(1,j)
      n2 = nodenum_at_face(2,j)
		
      if(top_layer_at_face(j) == 0) then	
      	! if this face is dry, set face normal & tangent velocities to 0.
         do k = 1, maxlayer
            un_face_new(k,j) = 0.0_dp
            vn_face_new(k,j) = 0.0_dp
         end do
      else 
      	! if this face is wet, calculate face normal & tangent velocities at each vertical layer (i.e., at vertical center of each layer) 
      	! from the momentum equation: Eq(5-1) in GOM manual
      	! note: it is better to locate if statement outside the vertical do loop (for 3D model setup),
      	! but I just keep this approach for the simplicity in cooding.
      	! In the future, I have to correct this for the efficiency. 
      	
         ! Eq(44) in our manuscript
      	! Normal velocity ------------------------------------------------!
    	   ! if(adj_cellnum_at_face(2,j) /= 0) then ! if this is not a boundary face
         do k = 1, top_layer_at_face(j) - bottom_layer_at_face(j) + 1
         	! note, "kk" is the vertical layer in reverse order (surface to bottom)
            kk = top_layer_at_face(j) + 1 - k 
            
            ! Normal velocity ------------------------------------------------!
    	      ! if(adj_cellnum_at_face(2,j) /= 0) then ! if this is not a boundary face
    	      if(boundary_type_of_face(j) == -1) then
    	      	! if this is a land face, no face normal velocity
    	      	un_face_new(kk,j) = 0.0_dp    	      	
    	      else if(boundary_type_of_face(j) == 0) then 
    	      	! if this is a normal water face,
    	      	
    	      	! let's kill this now...
    	      	! this is the new approach assuming nonorthogonal grid
    	      	! note: only true face-normal portion must be used for implicit barotropic gradient term,
    	      	! that is why we need cos_theta2(j)
!    	      	rtemp0 = theta*gravity*dt/delta_j(j) &
!               &	* ( eta_cell_new(adj_cellnum_at_face(2,j)) - eta_cell_new(adj_cellnum_at_face(1,j)) ) * AinvDeltaZ1(k,j)
               
               ! Dr. Sun's approach
               ! rtemp1 = rtemp0 / sin_theta2(j) ! true face-normal portion
               
               ! jw's approach
!               rtemp1 = rtemp0 * cos_theta2(j) ! true face-normal portion
!               rtemp2 = rtemp0 * sin_theta2(j) ! true face tangential portion, this will be used later
               
               ! note: The barotropic gradient term in the explicit term, i.e., AinvG1(j,k), 
               ! already includs only true face-normal portion (it is done in solve_momentum_equation_omp.f90). 
!               un_face_new(kk,j) = 	AinvG1(k,j) - rtemp1

					! this is the original approach assuming orthogonal grid
					un_face_new(kk,j) = 	&
					&		AinvG1(k,j) - theta*gravity*dt/delta_j(j) &
               &	* ( eta_cell_new(adj_cellnum_at_face(2,j)) - eta_cell_new(adj_cellnum_at_face(1,j)) ) &
               &	* AinvDeltaZ1(k,j)
            else if(boundary_type_of_face(j) > 0) then
		 			! if this is a open boundary face
		 			
    	      	! let's kill this now...
    	      	! this is the new approach assuming nonorthogonal grid
    	      	! Note: both methods should be identical since we are just using true-face normal velocities at open boundary face
		 	 		! note, at open boundary all portions are in the true-face normal direction, 
		 			! and that is why we must not include cos_theta2(j) here
! 		 			rtemp0 = theta*gravity*dt/delta_j(j) &
!               &	* ( eta_at_ob_new(ob_element_flag(adj_cellnum_at_face(1,j))) - eta_cell_new(adj_cellnum_at_face(1,j)) )	&
!               &	* AinvDeltaZ1(k,j)               
!               un_face_new(kk,j) = AinvG1(k,j) - rtemp0

					! this is the original approach assuming orthogonal grid
					un_face_new(kk,j) = &
					&		AinvG1(k,j) - theta*gravity*dt/delta_j(j) &
					&	* ( eta_at_ob_new(ob_element_flag(adj_cellnum_at_face(1,j))) - eta_cell_new(adj_cellnum_at_face(1,j)) )	&
					&	* AinvDeltaZ1(k,j)  
            end if
            
		 		! Tangential velocity --------------------------------------------!
		 		! note: 
		 		! 		wetdry_node == 0 (wet node)
		 		! 		wetdry_node == 1 (dry node)
            if(wetdry_node(n1) == 0 .and. wetdry_node(n2) == 0) then
            	! if both nodes are not dry (i.e., wet node)
               ! full portion of node to node tangential velocity

					! This is the new approach assuming non-orthogonal grid:
					! let's kill this now...
               ! Dr. Sun's approach
               ! rtemp0 = theta*gravity*dt/face_length(j) &
               ! &		* ( eta_node(n2) - eta_node(n1) )	&
               ! &		* AinvDeltaZ2(k,j)
               
               ! vn_face_new(kk,j) = AinvG2(k,j) - rtemp0

					! jw's approach
!               rtemp1 = theta*gravity*dt/face_length(j) &
!               &		* ( eta_node(n2) - eta_node(n1) )	&
!               &		* AinvDeltaZ2(k,j)
               
               ! Now, I have to include tangential portion which generated from cell to cell barotropic gradient
               ! note: the tangential portion, which generated from the explicit barotropic gradient term on the fake face-normal axis, 
               ! is already included in AinvG2(k,j) (it is don in solve_momentum_equation_omp.f90)
! 	         	if(boundary_type_of_face(j) == 0) then
	         		! if this is a normal water face,         		
	         		! tangential portion from the fake face-normal component must be included.
	         		! note: rtemp2 is previously calculated.
! 	         		vn_face_new(kk,j) = AinvG2(k,j) - (rtemp1 + rtemp2)
! 		         else if(boundary_type_of_face(j) > 0) then
		         	! if this is an open boundary face, 
		         	! no tangential portion from the fake face-normal component
! 		         	vn_face_new(kk,j) = AinvG2(k,j) - rtemp1
! 		         else if(boundary_type_of_face(j) == -1) then
		         	! if this is land boundary face,
		         	! only node to node tangential component is included
! 		         	vn_face_new(kk,j) = AinvG2(k,j) - rtemp1
! 		         end if
		         
		         ! This is the original approach assuming orthogonal grid:
		         vn_face_new(kk,j) = &
		         &		AinvG2(k,j) - theta*gravity*dt/face_length(j) &
               &	* ( eta_node(n2) - eta_node(n1) )	&
               &	* AinvDeltaZ2(k,j)
               
            end if

				! impose bounds for tang. vel. for initially dry nodes (h_node < 0)
            if(initial_wetdry_node(n1) == 1 .or. initial_wetdry_node(n2) == 1) then
            	! Shouldn't it be zero???
               vn_face_new(kk,j) = dmax1(-5.0_dp, dmin1(vn_face_new(kk,j), 5.0_dp))
            end if
         end do ! k = 1, top_layer_at_face(j)-bottom_layer_at_face(j)+1


         ! impose river b.c.
         if(isflowside3(j) > 0) then
         	ibnd = isflowside3(j)
         	do k = bottom_layer_at_face(j), top_layer_at_face(j)
         		un_face_new(k,j) = Qu_boundary(ibnd)
         	end do
         end if
         
         ! impose withdrawal/return b.c.
			if(isflowside4(j) > 0) then
				ibnd = isflowside4(j)
				if(WR_layer(ibnd) == 999) then
					! anytime surface layer
					k = top_layer_at_face(j)
					un_face_new(k,j) = WRu_boundary(ibnd)
				else if(WR_layer(ibnd) == 0) then
					! anytime bottom layer
					k = bottom_layer_at_face(j)
					! write(*,*) 'k,j,WRu_boundary(ibnd) = ', k,j,WRu_boundary(ibnd)
					un_face_new(k,j) = WRu_boundary(ibnd)
				else
					k = WR_layer(ibnd)
					un_face_new(k,j) = WRu_boundary(ibnd)
				end if
			end if


			! impose radiation b.c.
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
                       &            * eta_cell_new(adj_cellnum_at_face(1,j))
            end do
         end if
         
                  
			! set velocity to zero below the bottom layer (check this again)
			! do k = 0, bottom_layer_at_face(j)-1
         do k = 1, bottom_layer_at_face(j)-1
            un_face_new(k,j) = 0.0_dp  ! un_face_new(bottom_layer_at_face(i),i)
            vn_face_new(k,j) = 0.0_dp 
         end do

			! set above top layer velocities to top layer velocity 
         do k = top_layer_at_face(j)+1, maxlayer
            un_face_new(k,j) = un_face_new(top_layer_at_face(j),j)
            vn_face_new(k,j) = vn_face_new(top_layer_at_face(j),j)
         end do
      end if ! if(top_layer_at_face(j) == 0) then else       
   end do ! do j = 1, maxface
 	!$omp end do
	
	! Calculate w =============================================================!
	! vertical velocity at element center
	! solve eq(47), Lee et al. (2020)
 	!$omp do private(i,k,l,sum1,j_temp, t_layer, b_layer)
   do i = 1, maxele
   	t_layer = top_layer_at_element(i)
   	b_layer = bottom_layer_at_element(i)
   	
      if(t_layer == 0) then
         do k = 0, maxlayer
            wn_cell_new(k,i) = 0.0_dp
         end do
      else 
      	! if wet element, calculate velocity at both lower and upper prism face of each water layer.
      	! Eq(5-4) in GOM manual
      	! This is a original explicit method
!         wn_cell_new(bottom_layer_at_element(i)-1,i) = 0.0_dp ! no vertical flow at lower face of the bottom layer (i.e., boundary condition)         
!         do k = bottom_layer_at_element(i), top_layer_at_element(i)
         	! Water will stack up not stack down, that is why following relation is true:
         	! new vertical velocity at upper face = old vertical velocity at lower face + vertical velocity created by horizontal inflows through each face
!            wn_cell_new(k,i) = wn_cell_new(k-1,i)
!            do l = 1, tri_or_quad(i)
!               j_temp = facenum_at_cell(l,i)
!               if(k >= bottom_layer_at_face(j_temp) .and. k <= top_layer_at_face(j_temp)) then
!                  wn_cell_new(k,i) = wn_cell_new(k,i) &
!                  &           - sign_in_outflow(l,i) * face_length(j_temp) &
!                  &           * dz_face(k,j_temp) * un_face_new(k,j_temp) &
!                  &           / area(i)
!               end if
!            end do
!         end do
			
			! second method, same as the first method
			wn_cell_new(b_layer-1,i) = 0.0_dp ! no vertical flow at lower face of the bottom layer (i.e., boundary condition)
			
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
			
			! This is another method with semi-implicit approach....
! 			wn_cell_new(bottom_layer_at_element(i)-1,i) = 0.0_dp ! no vertical flow at lower face of the bottom layer (i.e., boundary condition)
! 			do k=bottom_layer_at_element(i),top_layer_at_element(i)
! 				sum1 = 0.0_dp
! 				do l=1,tri_or_quad(i)
! 					j_temp = facenum_at_cell(l,i)
! 					if(k >= bottom_layer_at_face(j_temp) .and. k <= top_layer_at_face(j_temp)) then
! 						sum1 = sum1 + sign_in_outflow(l,i)*face_length(j_temp) * dz_face(k,j_temp) &
! 						&		* (theta*un_face_new(k,j_temp) + (1.0-theta)*un_face(k,j_temp))
! 					end if
! 				end do
! 				wn_cell_new(k,i) = (1.0 - 1.0/theta)*wn_cell(k,i) + wn_cell(k-1,i) + (1.0/theta - 1.0)*wn_cell(k-1,i) &
! 				&						- (1.0/(theta*area(i)))*sum1
! 			end do
			
			! 	set velocity to zero below the bottom layer
         do k = 0, b_layer-2
            wn_cell_new(k,i) = 0.0_dp
         end do

			! set above top layer velocities to top layer velocity 
         do k = t_layer+1, maxlayer
            wn_cell_new(k,i) = wn_cell_new(t_layer,i)
         end do
      end if
   end do
 	!$omp end do
 	!$omp end parallel

	! write diagonostic file ==================================================!
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
	
	! write(*,'(*(F15.10))') (un_face(k,150), k=1,maxlayer)
   
end subroutine solve_velocities
