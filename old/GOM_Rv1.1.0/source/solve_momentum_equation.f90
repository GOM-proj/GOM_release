!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Calculate matrices and vectors for momentum equation, along each face.
!! This subroutine will calculate:
!! 		A: matrix A, 											Eq(44) in Lee et al., 2020 (or, (Eq(5-1) and (5-2) in GOM Manual) using Thomas Algorithm:)
!! 		G1(*): 	vector G for normal velocity, 		Eq(45) in Lee et al., 2020
!! 		G2(*): 	vector G for tangential velocity, 	Eq(45) in Lee et al., 2020 for tangential
!! 		AinvG1: 	normal velocity component,				Eq(46), just part
!! 		AinvG2: 	tangential velocity component,		Eq(46), just part
!! 		AinvDeltaZ1: normal velocity component,		Eq(46), just part
!! 		AinvDeltaZ2:tangential velocity component, 	Eq(46), just part
!! 
subroutine solve_momentum_equation
   use mod_global_variables
   use mod_file_definition
   implicit none

   integer :: i, j, k, kk, k2, l, ll, ie, nd, jf, ibnd
   integer :: ncyc, node1, node2, i_temp_element, num_vertical_layer
   real(dp):: hdiff ! jw
   real(dp):: arg, detp_dx, detp_dy, utmp, vtmp, sum1, sum2
   real(dp):: dudx, dvdx, dudy, dvdy, d2udy, d2vdy, d2udx, d2vdx, rl
   real(dp),dimension(12)  			:: temp 		! jw
   real(dp),dimension(maxlayer+1)  	:: G1, G2 	! jw
   real(dp),dimension(maxlayer+1)  	:: a_lower_mat, b_diagonal_mat, c_upper_mat, gam ! jw
   real(dp),dimension(maxlayer+1,5)	:: solution, rhs
   real(dp):: rtemp0, rtemp1, rtemp2, rtemp3, rtemp4
   ! jw
   
	! jw
	! jw
	! jw
	! jw
	AinvDeltaZ1 = 0.0_dp		! jw
	AinvDeltaZ2 = 0.0_dp		! jw
	AinvG1  		= 0.0_dp 	! jw
	AinvG2		= 0.0_dp		! jw
	
	! jw
	G1					= 0.0_dp
   G2					= 0.0_dp
	temp				= 0.0_dp
	a_lower_mat 	= 0.0_dp
	b_diagonal_mat	= 0.0_dp
	c_upper_mat 	= 0.0_dp
	gam 				= 0.0_dp
	solution			= 0.0_dp
	rhs				= 0.0_dp

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
	!$omp parallel do private(j,k,kk,l,ll,jf,node1,node2,num_vertical_layer) &
	!$omp& private(a_lower_mat,b_diagonal_mat,c_upper_mat,rhs,solution,gam) &
	!$omp& private(ie,temp,ncyc,arg,detp_dx,detp_dy,nd) &
	!$omp& private(G1,G2,utmp,vtmp,hdiff,dudx,dvdx,dudy,dvdy,d2udy,d2vdy,rl,d2udx,d2vdx) &
	!$omp& private(k2,sum1,sum2,rtemp0,rtemp1,rtemp2,rtemp3,rtemp4,i_temp_element,ibnd)
   do j = 1, maxface
		! jw
      if(top_layer_at_face(j) == 0) then
      	cycle ! jw
      end if
      
      ! jw
		! jw
		! jw
      node1 = nodenum_at_face(1,j) ! jw
      node2 = nodenum_at_face(2,j) ! jw
      num_vertical_layer = top_layer_at_face(j) - bottom_layer_at_face(j) + 1 ! jw

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
      do k = bottom_layer_at_face(j)+1, top_layer_at_face(j)-1
      	! jw
      	! jw
      	! jw
      	! jw
      	! jw
      	! jw
         kk = top_layer_at_face(j) - k + 1 									! jw
			
         a_lower_mat(kk)    = -dt * Av(k,  j)/dzhalf_face(k  ,j) 	! jw
         c_upper_mat(kk)    = -dt * Av(k-1,j)/dzhalf_face(k-1,j) 	! jw
         b_diagonal_mat(kk) = -a_lower_mat(kk) + dz_face(k,j) - c_upper_mat(kk)
      end do

		! jw
		! jw
		! jw
		! jw
      if(top_layer_at_face(j) == bottom_layer_at_face(j)) then ! jw
      	! jw
      	! jw
      	b_diagonal_mat(1) = dz_face(top_layer_at_face(j),j) + Gamma_B(j)*dt
      	
         ! jw
         ! jw
      else ! jw
			! jw
			! jw
         c_upper_mat(1) = -dt * Av(top_layer_at_face(j)-1,j) /	dzhalf_face(top_layer_at_face(j)-1,j)

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
         b_diagonal_mat(1) =  Gamma_T(j)*dt + dz_face(top_layer_at_face(j),j) - c_upper_mat(1)

			! jw
			! jw
         a_lower_mat(num_vertical_layer) = -dt * Av(bottom_layer_at_face(j),j) / dzhalf_face(bottom_layer_at_face(j),j)
         
         b_diagonal_mat(num_vertical_layer) = &
         &	-a_lower_mat(num_vertical_layer) + dz_face(bottom_layer_at_face(j),j) + Gamma_B(j)*dt ! jw
      end if
      ! jw
      
		! jw
      if(no_Etide_species > 0 .and. h_face(j) >= Etide_cutoff_depth .and. adj_cellnum_at_face(2,j) /= 0) then
         do l = 1, 2
            ie = adj_cellnum_at_face(l,j)
            temp(l) = 0.0
            do jf = 1, no_Etide_species
               ncyc  = int(Etide_frequency(jf)*elapsed_time/2.0d0/pi)
               arg   = Etide_frequency(jf)*elapsed_time-ncyc*2*pi   		&
               &     + Etide_species(jf)*lon_cell(ie)   						&
               &     + Etide_astro_arg_degree(jf)
               
               temp(l) = temp(l)   &
               &    + spinup_function_tide * Etide_amplitude(jf)      	&
               &                           * Etide_nodal_factor(jf)   	&
               &    * Etide_species_coef_at_element(Etide_species(jf),ie)   &
               &    * dcos(arg)	
            end do
         end do
         ! jw
         detp_dx = (temp(2)-temp(1))/delta_j(j)
      else
         detp_dx = 0.0d0
      end if

      if(no_Etide_species > 0 .and. h_face(j) >= Etide_cutoff_depth) then
         do l = 1, 2
            nd = nodenum_at_face(l,j)
            temp(l) = 0
            do jf = 1, no_Etide_species
               ncyc = int(Etide_frequency(jf)*elapsed_time/2.0d0/pi)
               arg  = Etide_frequency(jf)*elapsed_time-ncyc*2*pi   	&
               &    + Etide_species(jf)*lon_node(nd)*deg2rad			&
               &    + Etide_astro_arg_degree(jf)
               
               temp(l) = temp(l)   														&
               &       + spinup_function_tide * Etide_amplitude(jf)   		&
               &                              * Etide_nodal_factor(jf)  	&
               &       * Etide_species_coef_at_node(Etide_species(jf),nd)	&
               &       * dcos(arg)	
            end do
         end do
         detp_dy = (temp(2)-temp(1))/face_length(j)
      else
         detp_dy = 0.0
      end if
      
		! jw
		! jw
		! jw
		! jw
		! jw
		do k = 1, num_vertical_layer ! jw
      	! jw
         kk = top_layer_at_face(j) - k + 1	! jw

			! jw
			! jw
			! jw
			! jw

         !====================================================================!
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
			!====================================================================!
			! jw
			! jw
			! jw
			! jw
         ! jw
         ! jw
         ! jw

			! jw
			! jw
         G1(k) = dz_face(kk,j) * (un_ELM(kk,j) + coriolis_factor(j) * dt * vn_face(kk,j))
         G2(k) = dz_face(kk,j) * (vn_ELM(kk,j) - coriolis_factor(j) * dt * un_face(kk,j))
			
			! jw
			! jw
         if(adj_cellnum_at_face(2,j) /= 0) then ! jw
         	! jw
         	! jw
         	! jw
         	
         	! jw
         	! jw
         	! jw
         	! jw
         	! jw
            do l = 1, 2
               nd = nodenum_at_face(l,j)
               utmp = (u_node(kk,nd)+u_node(kk-1,nd))*0.5
               vtmp = (v_node(kk,nd)+v_node(kk-1,nd))*0.5
               ! jw
               temp(2*l-1) =  utmp*cos_theta(j)   &
                 &         +  vtmp*sin_theta(j) ! jw
               temp(2* l ) = -utmp*sin_theta(j)   &
                 &         +  vtmp*cos_theta(j) ! jw
            end do
            
            ! jw
         	! jw
         	! jw
         	! jw
         	! jw
            do l = 1, 2
               ie = adj_cellnum_at_face(l,j)
               utmp = 0.0
               vtmp = 0.0
               do ll = 1, tri_or_quad(ie)
                  nd   = nodenum_at_cell(ll,ie)
                  utmp =  utmp   &
                  &    + (u_node(kk,nd)+u_node(kk-1,nd))*0.5/tri_or_quad(ie) ! jw
                  vtmp =  vtmp   &
                  &    + (v_node(kk,nd)+v_node(kk-1,nd))*0.5/tri_or_quad(ie) ! jw
               end do
               ! jw
               temp(2*l+3) =  utmp*cos_theta(j)   &
                 &         +  vtmp*sin_theta(j) ! jw
               temp(2*l+4) = -utmp*sin_theta(j)   &
                 &         +  vtmp*cos_theta(j) ! jw
            end do

            dudx = (temp(7)-temp(5)) / delta_j(j)
            dvdx = (temp(8)-temp(6)) / delta_j(j)
            dudy = (temp(3)-temp(1)) / face_length(j)
            dvdy = (temp(4)-temp(2)) / face_length(j)
            
				! jw
            ! jw
            ! jw
            ! jw
            ! jw
            ! jw
            ! jw
            ! jw
            hdiff = smagorinsky_parameter   &
              &   * sqrt(dudx**2+dvdy**2+0.5*(dvdx+dudy)**2)
            
				! jw
				! jw
            Kh(k,j) = hdiff * area(adj_cellnum_at_face(1,j)) ! jw
            
            ! jw
            hdiff = hdiff + Ah_0
            Kh(k,j) = Kh(k,j) + Kh_0

				! jw
            d2udy = 4.0*(temp(3)+temp(1)-2*un_face(kk,j))/face_length(j)**2
            d2vdy = 4.0*(temp(4)+temp(2)-2*vn_face(kk,j))/face_length(j)**2

            do l = 1, 2
               ie = adj_cellnum_at_face(l,j)
               rl = sqrt(   (x_cell(ie)-x_face(j))**2   &
               &          + (y_cell(ie)-y_face(j))**2  )
               temp(2*l+7) = (temp(2*l+3)-un_face(kk,j))/rl ! jw
               temp(2*l+8) = (temp(2*l+4)-vn_face(kk,j))/rl ! jw
            end do

            d2udx = 2.0*(temp(9 )+temp(11)) / delta_j(j)
            d2vdx = 2.0*(temp(10)+temp(12)) / delta_j(j)
				
				G1(k) = G1(k)+dz_face(kk,j)*dt*hdiff*(d2udx+d2udy)
            G2(k) = G2(k)+dz_face(kk,j)*dt*hdiff*(d2vdx+d2vdy)            
         end if ! jw
         
			! jw
			! jw
			! jw
			! jw
         if(boundary_type_of_face(j) == 0) then ! jw
         	! jw
         	! jw
         	! jw
         	! jw
         	rtemp0= gravity * dt / delta_j(j) & 
         	&		* (1.0-theta) * (eta_cell(adj_cellnum_at_face(2,j)) - eta_cell(adj_cellnum_at_face(1,j)))
         	
         	rtemp1= rtemp0 * cos_theta2(j) ! jw
         	rtemp2= rtemp0 * sin_theta2(j) ! jw
            
            ! jw
            G1(k) = G1(k) - dz_face(kk,j) * rtemp1
         else if(boundary_type_of_face(j) > 0) then
         	! jw
         	! jw
         	! jw
         	! jw
         	! jw
         	ie = ob_element_flag(adj_cellnum_at_face(1,j)) ! jw
         	rtemp1= gravity * dt / delta_j(j) &
         	&		* (1.0-theta) * (eta_at_ob_old(ie) - eta_cell(adj_cellnum_at_face(1,j)))
            G1(k) = G1(k) - dz_face(kk,j) * rtemp1
         else if(boundary_type_of_face(j) == -1) then
         	! jw
         	G1(k) = 0.0_dp ! jw
         end if

         ! jw
         ! jw
         if(wetdry_node(node1) == 0 .and. wetdry_node(node2) == 0) then ! jw
         	! jw
         	! jw
        		rtemp0= gravity * dt / face_length(j) &
        		&		* (1.0-theta) * (eta_node(node2) - eta_node(node1))

         	if(boundary_type_of_face(j) == 0) then
         		! jw
         		! jw
         		G2(k) = G2(k) - dz_face(kk,j) * (rtemp0 + rtemp2)
	         else if(boundary_type_of_face(j) > 0) then
	         	! jw
	         	G2(k) = G2(k) - dz_face(kk,j) * rtemp0
	         else if(boundary_type_of_face(j) == -1) then
	         	! jw
	         	G2(k) = G2(k) - dz_face(kk,j) * rtemp0
	         end if
         end if         
         
			! jw
         if(baroclinic_flag == 1) then
          	sum1 = 0.0_dp ! jw
         	sum2 = 0.0_dp ! jw
         	
            do k2=kk,top_layer_at_face(j)
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
               if(adj_cellnum_at_face(2,j) /= 0 .and.   &
            	&  k2 >= bottom_layer_at_element(adj_cellnum_at_face(1,j)) .and.   &
            	&  k2 <= top_layer_at_element(adj_cellnum_at_face(1,j))    .and.   &
            	&  k2 >= bottom_layer_at_element(adj_cellnum_at_face(2,j)) .and.   &
            	&  k2 <= top_layer_at_element(adj_cellnum_at_face(2,j))) then						
						if(k2 == kk) then
							! jw
							sum1 = sum1 + &
							&	0.5_dp*(rho_cell(k2,adj_cellnum_at_face(2,j)) - rho_cell(k2,adj_cellnum_at_face(1,j)))*dz_face(k2,j)
							! jw
						else
							! jw
							sum1 = sum1 + &
							&	1.0_dp*(rho_cell(k2,adj_cellnum_at_face(2,j)) - rho_cell(k2,adj_cellnum_at_face(1,j)))*dz_face(k2,j)
							! jw
						end if
               end if
               
               ! jw
               ! jw
               if(k2 >= bottom_layer_at_node(node1)  .and.   &
               &  k2 <= top_layer_at_node(node1) .and.   &
               &  k2 >= bottom_layer_at_node(node2)  .and.   &
               &  k2 <= top_layer_at_node(node2)) then						
						if(k2 == kk) then
							! jw
							sum2 = sum2 + 0.5_dp*(rho_node(k2,node2)-rho_node(k2,node1))*dz_face(k2,j)
							! jw
						else
							! jw
							sum2 = sum2 + 1.0_dp*(rho_node(k2,node2)-rho_node(k2,node1))*dz_face(k2,j)
							! jw
						end if
               end if 
            end do ! jw

				! jw
 				rtemp1 = dz_face(k,j)*(-gravity*dt/(rho_o*delta_j(j)))*sum1 ! jw
 				rtemp2 = dz_face(k,j)*(-gravity*dt/(rho_o*face_length(j)))*sum2 ! jw
				
				! jw
! jw
! jw
! jw
! jw

				! jw
 				rtemp3 = rtemp1 * cos_theta2(j) ! jw
 				rtemp4 = rtemp1 * sin_theta2(j) ! jw
				
				! jw
 				G1(k) = G1(k) + rtemp3 * spinup_function_baroclinic ! jw
 				G2(k) = G2(k) + (rtemp2 + rtemp4) * spinup_function_baroclinic
         end if ! jw

         ! jw
         if(k == 1) then
         	G1(k) = G1(k) + dt*wind_stress_normal(j)	! jw
            G2(k) = G2(k) + dt*wind_stress_tangnt(j)	! jw
         end if

			! jw
			! jw
			! jw
			! jw
			! jw
			! jw
         G1(k) = G1(k) + 0.69*dz_face(kk,j)*gravity*dt*detp_dx
         G2(k) = G2(k) + 0.69*dz_face(kk,j)*gravity*dt*detp_dy

			! jw
         ! jw
         if(adj_cellnum_at_face(2,j) /= 0) then ! jw
            do l = 1, 2
               temp(l) = 0.0
               i_temp_element = adj_cellnum_at_face(l,j)
               do ll = 1, tri_or_quad(i_temp_element) ! jw
                  temp(l) = temp(l)                                                  &
                    &     + air_p1(nodenum_at_cell(ll,i_temp_element))   &
                    &     / tri_or_quad(i_temp_element)
               end do
            end do
				! jw
				rtemp0 = -dt/rho_o*(temp(2)-temp(1))/delta_j(j)
         end if
         rtemp1 = rtemp0 * cos_theta2(j) ! jw
         rtemp2 = rtemp0 * sin_theta2(j) ! jw
         G1(k) = G1(k) + dz_face(kk,j) * rtemp1
         
			! jw
			rtemp3 = -dt/rho_o*((air_p1(node2)-air_p1(node1)) / face_length(j))
			G2(k) = G2(k) + dz_face(kk,j) * (rtemp3 + rtemp2)
         ! jw

         
         !====================================================================!
         ! jw
         ! jw
         !====================================================================!

			! jw
         if(isflowside(j) > 0) then ! jw
				! jw
            ibnd    =  isflowside(j) ! jw
            ! jw
            
            ! jw
            ! jw
            ! jw
            G1(k) =  u_boundary(ibnd)*cos_theta(j)   &
            &     +  v_boundary(ibnd)*sin_theta(j)
            G2(k) = -u_boundary(ibnd)*sin_theta(j)   &
            &     +  v_boundary(ibnd)*cos_theta(j)
         else if(boundary_type_of_face(j) == -1) then 
				! jw
         	! jw
            G1(k) = 0.0 ! jw
         end if
         
         ! jw
         ! jw
         ! jw
         ! jw
         ! jw
         ! jw
         ! jw
         if(isflowside3(j) > 0) then
 			 	ibnd = isflowside3(j)
 				
 				! jw
 				! jw
            G1(k) = G1(k) + dz_face(kk,j) * Qu_boundary(ibnd)
            !! G2(k) =  G2(k) + dz_face(kk,j) * Qv_boundary(ibnd)
			end if
			
			! jw
			if(isflowside4(j) > 0) then
				ibnd = isflowside4(j)
				if(WR_layer(ibnd) == 999) then
					! jw
					if(kk == top_layer_at_face(j)) then
						G1(k) = G1(k) + dz_face(kk,j) * WRu_boundary(ibnd)
					end if
				else if(WR_layer(ibnd) == 0) then
					! jw
					if(kk == bottom_layer_at_face(j)) then
						! jw
						G1(k) = G1(k) + dz_face(kk,j) * WRu_boundary(ibnd)
					end if
				else
					if(kk == WR_layer(ibnd)) then
						G1(k) = G1(k) + dz_face(kk,j) * WRu_boundary(ibnd)
					end if
				end if
			end if


			! jw
         rhs(k,1) = G1(k)
         rhs(k,2) = G2(k)
         rhs(k,3) = dz_face(kk,j)
		end do ! jw

		! jw
		! jw
		! jw
		! jw
		
		! jw
      call tridiagonal_solver(maxlayer+1, num_vertical_layer, 3, a_lower_mat, b_diagonal_mat, c_upper_mat, rhs, solution, gam)
		
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
		do k=1,num_vertical_layer
			AinvG1(k,j) = solution(k,1)
			AinvG2(k,j) = solution(k,2)
		end do
		
		! jw
		! jw
		if(isflowside(j) > 0) then
			do k=1,num_vertical_layer
 				AinvG1(k,j) = G1(k)
 				AinvG2(k,j) = G2(k)				
			end do
		end if
		
		! jw
		! jw
		! jw
		if(boundary_type_of_face(j) == -1 .and. isflowside3(j) == 0 .and. isflowside4(j) == 0) then
			do k=1,num_vertical_layer
 				AinvG1(k,j) = 0.0_dp
 				! jw
			end do			
		end if
		
		
 		! jw
		! jw
		! jw
		! jw
		do k=1,num_vertical_layer
			AinvDeltaZ1(k,j) = solution(k,3)
 			AinvDeltaZ2(k,j) = solution(k,3)
 		end do
 		
		! jw
		! jw
 		if(isflowside(j) > 0) then
 			do k=1,num_vertical_layer
				AinvDeltaZ1(k,j) = 0.0_dp
 				AinvDeltaZ2(k,j) = 0.0_dp
 			end do
 		end if
 		
		! jw
		! jw
		! jw
		if(boundary_type_of_face(j) == -1) then
			do k=1,num_vertical_layer
 				AinvDeltaZ1(k,j) = 0.0_dp
 				! jw
			end do
		end if
		
	end do ! jw
	!omp end parallel do

	! jw
	if(dia_momentum == 1) then
		write(pw_dia_momentum,*) 'it = ', it, ', elapsed_time = ', elapsed_time		
		write(pw_dia_momentum,'(A)') &
		&	'j, AinvG1(maxlayer,maxface), AinvG2(maxlayer,maxface), AinvDeltaZ1(maxlayer,maxface), AinvDeltaZ2(maxlayer,maxface)'
		do j=1,maxface
	 		write(pw_dia_momentum,'(A, I5, *(E15.5))') 'j= ', j, &
	 		&	(AinvG1(k,j), k=1,maxlayer), &
	 		&	(AinvG2(k,j), k=1,maxlayer), &
	 		&	(AinvDeltaZ1(k,j), k=1,maxlayer), &
	 		&	(AinvDeltaZ2(k,j), k=1,maxlayer)
	 	end do
	end if	
end subroutine solve_momentum_equation
