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
   real(dp):: hdiff ! Ah; i.e., horizontal diffusion coefficient
   real(dp):: arg, detp_dx, detp_dy, utmp, vtmp, sum1, sum2
   real(dp):: dudx, dvdx, dudy, dvdy, d2udy, d2vdy, d2udx, d2vdx, rl
   real(dp),dimension(12)  			:: temp 		! array to store partial derivative terms for horizontal diffusion term
   real(dp),dimension(maxlayer+1)  	:: G1, G2 	! G vector for normal & tangential velocity
   real(dp),dimension(maxlayer+1)  	:: a_lower_mat, b_diagonal_mat, c_upper_mat, gam ! in A matrix
   real(dp),dimension(maxlayer+1,5)	:: solution, rhs
   real(dp):: rtemp0, rtemp1, rtemp2, rtemp3, rtemp4
   ! End of local variables ==================================================!
   
	! re-initialize following global variables:
	! Note: 
	! 		I defined big (maxface related) arrays as global variables, 
	! 		but small (only maxlayer related) arrays as local variables.
	AinvDeltaZ1 = 0.0_dp		! inv(A)dz_face, [maxlayer,maxface], defined in mod_global_variables for normal velocity
	AinvDeltaZ2 = 0.0_dp		! inv(A)dz_face, [maxlayer,maxface], defined in mod_global_variables for tangential velocity
	AinvG1  		= 0.0_dp 	! inv(A)G1, [maxlayer,maxface], for normal velocity, defined in mod_global_variables
	AinvG2		= 0.0_dp		! inv(A)G2, [maxlayer,maxface], for tangential velocity, defined in mod_global_variables
	
	! initialize local variables:
	G1					= 0.0_dp
   G2					= 0.0_dp
	temp				= 0.0_dp
	a_lower_mat 	= 0.0_dp
	b_diagonal_mat	= 0.0_dp
	c_upper_mat 	= 0.0_dp
	gam 				= 0.0_dp
	solution			= 0.0_dp
	rhs				= 0.0_dp

	! include air pressure model for atmospheric pressure term ================!
	! convert unit from [milibar] to [N/m^2], 1[mb] = 100[N/m^2] = 100[kg/m/s^2]
	! in the equation, the unit of air pressure is N/m^2, so need to convert unit
!   if(airp_flag == 1) then
!   	do i = 1, maxnod	
!         air_p1(i) = air_p0(i)*100.0_dp
!      end do
!   else if(airp_flag == 2) then
!      do i = 1, maxnod
!         air_p1(i) = air_p0(i)*100.0_dp
!      end do
!   end if

	! momentum equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
	! momentum equation will be calculated at each face (with vertical layers in the face)
	!$omp parallel do private(j,k,kk,l,ll,jf,node1,node2,num_vertical_layer) &
	!$omp& private(a_lower_mat,b_diagonal_mat,c_upper_mat,rhs,solution,gam) &
	!$omp& private(ie,temp,ncyc,arg,detp_dx,detp_dy,nd) &
	!$omp& private(G1,G2,utmp,vtmp,hdiff,dudx,dvdx,dudy,dvdy,d2udy,d2vdy,rl,d2udx,d2vdx) &
	!$omp& private(k2,sum1,sum2,rtemp0,rtemp1,rtemp2,rtemp3,rtemp4,i_temp_element,ibnd)
   do j = 1, maxface
		! If it is a dry face, do not calculate momentum equation.
      if(top_layer_at_face(j) == 0) then
      	cycle ! do nothing and go to the next face
      end if
      
      ! following is only for wet face: i.e., if(top_layer_at_face(j) /= 0) then
		! Only if it is a wet face, calcuate momentum equation.
		! Find two nodes and total number of vertical layers at this face.
      node1 = nodenum_at_face(1,j) ! starting (1st) node for this face
      node2 = nodenum_at_face(2,j) ! ending (2nd) node for this face
      num_vertical_layer = top_layer_at_face(j) - bottom_layer_at_face(j) + 1 ! number of vertical layers at this face

		! ======================================================================!
		!      								 n                                   
		! setup coefficient  matrix  A           in Casulli & Walters, 2000, p338
		!                              j,k+-1/2                            
		! ======================================================================!
		! 
		! First, setup interior elements of matrix A:
		! 		note for numbering for vertical direction
		!     matrix number is reversed for top and bottom layer; see more details in Casulli & Walters, 2000, p337-338 (jw)
		!     top layer is always 1 and bottom layer is equal to actual surface layer number
		!     i.e. bottom layer = actual vertical layer number, 1 = actual top layer
		!     this is for just convinience for solving tridiagonal matrix
		
		! ----------------------------------------------------------------------!
		! setup diagonal elements of matrix A: a, b, & c
		! A = [b c 0 0 0 0 		-> surface layer
		!      a b c 0 0 0 
		!      0 a b c 0 0
		!      - - - - - -
		!      0 0 0 a b c
		!      0 0 0 0 a b]		-> bottom layer
		! 
		! For example, if there are 10 z-layers, and the top layer is 10 and the bottom layer is 9, for the current cell,
		! Then, the matrix will be:
		! A = [b c
		!      a b]
		! b & x = [10th z_layer vector/solution
		!          9th  z_layer vector/solution]
		! ----------------------------------------------------------------------!
		! Note, we are constructing this matrix at vertical layers not vertical levels.
		! First, setup A matrix at the middle of the water column (mid-layers)
		! If it is a 2D or if it only has two layers, 
		! this part will not be performed since it only has top & bottom layers (not mid-layers)
		
		! From here, I will use "kk" for matrix ordering, and "k" for corresponding variable's layer ordering
      do k = bottom_layer_at_face(j)+1, top_layer_at_face(j)-1
      	! Here, Matrix A is in the order as we see water colum (i.e.):
      	! 		[A] = [top_layer				-> kk = 1
      	! 				 top_layer - 1 		-> kk = 2
      	! 				 ...						-> ...
      	!            bottom_layer + 1		-> kk = top_layer_at_face(j) - 1
      	!            bottom_layer    ]	-> kk = top_layer_at_face(j)
         kk = top_layer_at_face(j) - k + 1 									! this is the vertical layer in variables, from top to bottom (reverse order)
			
         a_lower_mat(kk)    = -dt * Av(k,  j)/dzhalf_face(k  ,j) 	! jw, Av(k,j)   -> Av(k+1/2) in equation (at upper level of the layer)
         c_upper_mat(kk)    = -dt * Av(k-1,j)/dzhalf_face(k-1,j) 	! jw, Av(k-1,j) -> Av(k-1/2) in equation (at bottom level of the layer)
         b_diagonal_mat(kk) = -a_lower_mat(kk) + dz_face(k,j) - c_upper_mat(kk)
      end do

		! Second, setup the top and bottom layer's elements of matrix A
		! if the layer is just 2d, use bottom layer only; Casulli & Walters, 2000, p339
		! That is why we have to exclude 'Gamma_T * dt' term for 2D simulation
		! And, the wind stress term will be located in G vector (which is the known values)
      if(top_layer_at_face(j) == bottom_layer_at_face(j)) then ! if it is 2D
      	! if it is 2D, we only need the diagonal matrix term, b_diagonal_mat
      	! This is correct.
      	b_diagonal_mat(1) = dz_face(top_layer_at_face(j),j) + Gamma_B(j)*dt
      	
         ! Do not use this...
         ! b_diagonal_mat(1) = Gamma_T(j)*dt + dz_face(top_layer_at_face(j),j) + Gamma_B(j)*dt
      else ! if it is 3D
			! at the top(surface) layer -----------------------------------------!
			! Only b and c terms are required
         c_upper_mat(1) = -dt * Av(top_layer_at_face(j)-1,j) /	dzhalf_face(top_layer_at_face(j)-1,j)

			! jw, check this again.... where is gamma_t term???
			! see more details in Eq(4-42) & (5-2) in GOM manual
			! original
         ! b_diagonal_mat(1) =  dz_face(top_layer_at_face(j),j) - c_upper_mat(1)
         
         ! jw's correction; for the 3D model we have to include 'Gamma_T*dt' term 
         ! since we have to use the difference between wind speed and surface water speed.
         ! So, if we include Gamma_T*dt term here, the exerting wind stress will be slightly weaker 
         ! than the result with excluding this term (I confirmed it).
         ! However, if analytical wind stress is used, there is no way to include Gamma_T*dt term.
         ! Thus, simply ignore Gamma_T*dt term for analytical solution (it will be very small difference).
         b_diagonal_mat(1) =  Gamma_T(j)*dt + dz_face(top_layer_at_face(j),j) - c_upper_mat(1)

			! at the bottom layer -----------------------------------------------!
			! Only a and b terms are requried
         a_lower_mat(num_vertical_layer) = -dt * Av(bottom_layer_at_face(j),j) / dzhalf_face(bottom_layer_at_face(j),j)
         
         b_diagonal_mat(num_vertical_layer) = &
         &	-a_lower_mat(num_vertical_layer) + dz_face(bottom_layer_at_face(j),j) + Gamma_B(j)*dt ! -a_low_mat = +(a_term)			
      end if
      ! End of matrix A ======================================================!
      
		! pre-compute tidal potential gradients for each side to save time =====!
      if(no_Etide_species > 0 .and. h_face(j) >= Etide_cutoff_depth .and. adj_cellnum_at_face(2,j) /= 0) then
         do l = 1, 2
            ie = adj_cellnum_at_face(l,j)
            temp(l) = 0.0_dp
            do jf = 1, no_Etide_species
               ncyc  = int(Etide_frequency(jf)*elapsed_time/2.0_dp/pi)
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
         ! here temp(2) and temp(1) are the unit of eta, i.e, [m]; i.e., eta = a*cos(theta) like term.
         detp_dx = (temp(2)-temp(1))/delta_j(j)
      else
         detp_dx = 0.0_dp
      end if

      if(no_Etide_species > 0 .and. h_face(j) >= Etide_cutoff_depth) then
         do l = 1, 2
            nd = nodenum_at_face(l,j)
            temp(l) = 0.0_dp
            do jf = 1, no_Etide_species
               ncyc = int(Etide_frequency(jf)*elapsed_time/2.0_dp/pi)
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
         detp_dy = 0.0_dp
      end if
      
		! ======================================================================!
		!                         n          n
		!     calculate vectors Fu     and  G  in page 337 in (Casulli and Walters, 2000) paper
		!                         j,m        j
		! ======================================================================!		
		do k = 1, num_vertical_layer ! do k = 1, top_layer_at_face(j)-bottom_layer_at_face(j)+1
      	! surface to bottom order...
         kk = top_layer_at_face(j) - k + 1	! this is the vertical layer for variables (from top to bottom)

			! b.c. to be imposed at the end
			! coriolis, advection and wind stress:
			! 		G1 is the normal velocity component (for u*)
			! 		G2 is the tangential velocity component (for v*)

         !====================================================================!
         ! Now we will include following forcing terms in G1 & G2 vectors:
         ! 		(0) Previous velocity (explicit velocity)
         ! 		(1) Nonlinear advection
         ! 		(2) Coriolis
         ! 		(3) Horizontal diffusion
         ! 		(4) Barotropic gradient
         ! 		(5) Baroclinic gradient
         ! 		(6) Wind stress
         ! 		(7) Earth tidal potential
         ! 		(8) Atmospheric pressure gradient
			!====================================================================!
			! include (0)previous velocity (i.e. tendency term) & (1)nonlinear advection & (2)Coriolis term -----------------!
			! Originally, it should be like below, but explicit velocity (tendency term) is already included (coupled) in un_ELM & vn_ELM 
			! in solve_nonlinear_advection.f90 & main.f90.
			! So, when including other terms, the tendency term should not be included again.
         ! G1(k) = dz_face(kk,j) * (un_face(kk,j) + un_ELM(kk,j) + coriolis_factor(j) * dt * vn_face(kk,j))
         ! G2(k) = dz_face(kk,j) * (vn_face(kk,j) + vn_ELM(kk,j) - coriolis_factor(j) * dt * un_face(kk,j))
         ! i.e., un_vace(j,kk) + un_ELM(kk,j) -> is already in un_LEM(kk,j)

			! That is why we have to use this...
			! In GOM, explicit velocity is included in un_ELM & vn_ELM already
         G1(k) = dz_face(kk,j) * (un_ELM(kk,j) + coriolis_factor(j) * dt * vn_face(kk,j))
         G2(k) = dz_face(kk,j) * (vn_ELM(kk,j) - coriolis_factor(j) * dt * un_face(kk,j))
			
			! include (3)horizontal diffusion term ------------------------------!
			! smagorinsky diffusivity at half levels 
         if(adj_cellnum_at_face(2,j) /= 0) then ! if adjacent cell exists, include horizontal diffusion. Otherwise, do not include horizontal diffusion.
         	! Note, in this approach, we are not including horizontal diffusion term between open boundary element and the ghost cell
         	
         	! we have to find velocity information at 4 points:
         	! 		2 at face nodes (beginning and ending nodes) at each face
         	! 		2 at element centers (which share this face)
         	
         	! first, find velocities at 2 nodes at this face
         	! 		temp(1) = normal velocity at node#1
         	! 		temp(2) = tangent velocity at node#1
         	! 		temp(3) = normal velocity at node#2
         	! 		temp(4) = tangent velocity at node#2
            do l = 1, 2
               nd = nodenum_at_face(l,j)
               utmp = (u_node(kk,nd)+u_node(kk-1,nd))*0.5_dp
               vtmp = (v_node(kk,nd)+v_node(kk-1,nd))*0.5_dp
               ! change to face nomal/tangent velocities
               temp(2*l-1) =  utmp*cos_theta(j)   &
                 &         +  vtmp*sin_theta(j) ! temp(1) & (3); normal velocity at this vertical layer
               temp(2* l ) = -utmp*sin_theta(j)   &
                 &         +  vtmp*cos_theta(j) ! temp(2) & (4); tangent velocity at this vertical layer
            end do
            
            ! second, find velocities at 2 adjacent element centers which share this face
         	! 		temp(5) = normal velocity at element#1
         	! 		temp(6) = tangent velocity at element#1
         	! 		temp(7) = normal velocity at element#2
         	! 		temp(8) = tangent velocity at element#2
            do l = 1, 2
               ie = adj_cellnum_at_face(l,j)
               utmp = 0.0_dp
               vtmp = 0.0_dp
               do ll = 1, tri_or_quad(ie)
                  nd   = nodenum_at_cell(ll,ie)
                  utmp =  utmp   &
                  &    + (u_node(kk,nd)+u_node(kk-1,nd))*0.5_dp/tri_or_quad(ie) ! average nodal velocities at this vertical layer
                  vtmp =  vtmp   &
                  &    + (v_node(kk,nd)+v_node(kk-1,nd))*0.5_dp/tri_or_quad(ie) ! average nodal velocities at this vertical layer
               end do
               ! change to face nomal/tangent velocities
               temp(2*l+3) =  utmp*cos_theta(j)   &
                 &         +  vtmp*sin_theta(j) ! temp(5) & (7); normal velocity at cell center at this vertical layer
               temp(2*l+4) = -utmp*sin_theta(j)   &
                 &         +  vtmp*cos_theta(j) ! temp(6) & (8); tangent velocity at cell center at this vertical layer
            end do

            dudx = (temp(7)-temp(5)) / delta_j(j)
            dvdx = (temp(8)-temp(6)) / delta_j(j)
            dudy = (temp(3)-temp(1)) / face_length(j)
            dvdy = (temp(4)-temp(2)) / face_length(j)
            
				! Now, we are ready to solve Smagorinsky model.            
            ! Note: here, we have to multiply "area" to get Ah.
            ! However, when we calculate the second derivative terms later, we also have to divide it with the "area".
            ! And, that is why I omit the "area" here.
            ! i.e., this is the original Ah:
            ! hdiff = smagorinsky_parameter*area(adj_cellnum_at_face(1,j))   &
            ! &			* dsqrt(dudx**2+dvdy**2+0.5*(dvdx+dudy)**2)
            ! but, use this instead:
            hdiff = smagorinsky_parameter   &
              &   * sqrt(dudx**2+dvdy**2+0.5_dp*(dvdx+dudy)**2)
            
				! And, prepare for the eddy diffusivity
				! Note, for Kh, I have to keep the "area" term            
            Kh(k,j) = hdiff * area(adj_cellnum_at_face(1,j)) ! i.e., Ah = Kh
            
            ! Now, add background viscosity & diffusivity
            hdiff = hdiff + Ah_0
            Kh(k,j) = Kh(k,j) + Kh_0

				! Now, we have to calculate the second derivative terms: d^2u/dx^2 & d^2u/dy^2
            d2udy = 4.0_dp*(temp(3)+temp(1)-2*un_face(kk,j))/face_length(j)**2
            d2vdy = 4.0_dp*(temp(4)+temp(2)-2*vn_face(kk,j))/face_length(j)**2

            do l = 1, 2
               ie = adj_cellnum_at_face(l,j)
               rl = sqrt(   (x_cell(ie)-x_face(j))**2   &
               &          + (y_cell(ie)-y_face(j))**2  )
               temp(2*l+7) = (temp(2*l+3)-un_face(kk,j))/rl ! temp(9) & (11)
               temp(2*l+8) = (temp(2*l+4)-vn_face(kk,j))/rl ! temp(10) & (12)
            end do

            d2udx = 2.0_dp*(temp(9 )+temp(11)) / delta_j(j)
            d2vdx = 2.0_dp*(temp(10)+temp(12)) / delta_j(j)
				
				G1(k) = G1(k)+dz_face(kk,j)*dt*hdiff*(d2udx+d2udy)
            G2(k) = G2(k)+dz_face(kk,j)*dt*hdiff*(d2vdx+d2vdy)            
         end if ! if(adj_cellnum_at_face(2,j) /= 0) then
         
			! include (4)barotropic gradient ------------------------------------!
			! normal gradient (cell to cell between two adjacent cells which share this side)
			! Note: we don't need to calculate face normal velocity at land boundary face.
			! this is for faces which two water cells share
			! boundary_type_of_face: -1 (land face), 0 (inner water), 1 (open boundary)
         ! if(boundary_type_of_face(j) == 0) then ! or you can use: if(adj_cellnum_at_face(2,j) /= 0) then
         if(adj_cellnum_at_face(2,j) /= 0) then
         	! if an adjacent cell exists (i.e., water face), calculate barotropic gradient between two cells (only for normal component)
         	! note: when using fake face-normal axis, we must separate the fake-normal component to the true face-normal component.
         	! that is why I include cos_theta2(j) at the end of this equation.
         	! the sin_theta2(j) portion will be included in the tangential component, i.e., G2
         	
         	! this is the new approach with fake face-normal axis:
!         	rtemp0= gravity * dt / delta_j(j) & 
!         	&		* (1.0_dp-theta) * (eta_cell(adj_cellnum_at_face(2,j)) - eta_cell(adj_cellnum_at_face(1,j)))
         	
!         	rtemp1= rtemp0 * cos_theta2(j) ! true face-normal component
!         	rtemp2= rtemp0 * sin_theta2(j) ! tangential component, this will be added in G2(k)
            
            ! only include face-normal omponent
!           G1(k) = G1(k) - dz_face(kk,j) * rtemp1
            
            ! this is the original approach assuming orthogonal grid
         	G1(k) = G1(k) - dz_face(kk,j) * gravity * dt / delta_j(j) & 
         	&		* (1.0_dp-theta) * (eta_cell(adj_cellnum_at_face(2,j)) - eta_cell(adj_cellnum_at_face(1,j)))
            
         else if(boundary_type_of_face(j) > 0) then
         	! if this is an open boundary face
         	! calculate barotropic gradient between this cell and the ghost boundary cell (only for normal component)
         	! Here, eta_at_ob_old(ie) is the correcponding water surface elevation at the ghost boundary element (at the previous time step)
         	! note: at the boundary face, all force acts perpendicular to the face (i.e., true face-normal direction).
         	! thus, we should not include cos_theta2(j) here; i.e., all portion is perpendicular to the face
         	ie = ob_element_flag(adj_cellnum_at_face(1,j)) ! ie = ith open boundary element
            G1(k) = G1(k) - dz_face(kk,j) * gravity * dt / delta_j(j) &
         	&		* (1.0_dp-theta) * (eta_at_ob_old(ie) - eta_cell(adj_cellnum_at_face(1,j)))
         	
         ! the following statement is redundant, and thus deactivated.
!          else if(boundary_type_of_face(j) == -1) then
         	! land boundary face
!         	G1(k) = G1(k) + 0.0_dp ! no face-normal barotropic gradient exists. 
         end if

         ! tangential gradient (node to node between two nodes in this side)
         ! Note: tangential gradient can be calculated at all face.
         if(wetdry_node(node1) == 0 .and. wetdry_node(node2) == 0) then ! wetdry_node = 0 means it is a wet node
         	! if both beginning and ending nodes for this side are wet,
         	! calculate barotropic gradient between two nodes (only for tangential component)
         	
         	! this is the original method assuming orthogonal mesh:
            G2(k) = G2(k) - dz_face(kk,j) * gravity * dt / face_length(j) &
            &     * (1.0_dp-theta) * (eta_node(node2) - eta_node(node1))
         	
         	! this is the new method assuming non-orthogonal mesh:
!        		rtemp0= gravity * dt / face_length(j) &
!        		&		* (1.0-theta) * (eta_node(node2) - eta_node(node1))

!         	if(boundary_type_of_face(j) == 0) then
         		! if this is a normal water face,         		
         		! tangential portion from the fake face-normal component must be included.
!         		G2(k) = G2(k) - dz_face(kk,j) * (rtemp0 + rtemp2)
! 	         else if(boundary_type_of_face(j) > 0) then
	         	! if this is an open boundary face, there is no fake-normal component
! 	         	G2(k) = G2(k) - dz_face(kk,j) * rtemp0
! 	         else if(boundary_type_of_face(j) == -1) then
	         	! there is no fake-normal component
! 	         	G2(k) = G2(k) - dz_face(kk,j) * rtemp0
! 	         end if
         end if         
         
			! include (5)baroclinic gradient term -------------------------------!
         if(baroclinic_flag == 1) then
          	sum1 = 0.0_dp ! for face normal gradient
         	sum2 = 0.0_dp ! for tangential gradient
         	
            do k2=kk,top_layer_at_face(j)
            	! face normal baroclinic gradient; i.e., cell_center to cell_center
            	! If adjacent cell exists, and if this layer is located between top and bottom layers of both cells,
            	! i.e, only through same vertical layers, calculate baroclinic gradient between two cells.
					
	         	! note: as in barotropic gradient term,
	         	! when using fake face-normal axis, we must separate the fake-normal component to the true face-normal component.
	         	! that is why I include cos_theta2(j) at the end of this equation.
	         	! Note: we are calculating baroclinic gradient through this vertical layer of the face 
	         	! only if two jth-side sharing cells have the same vertical layer.
	         	! i.e., if the first cell's bottom_layer_at_element(j,1) == 3,
	         	!          the second cell's bottom_layer_at_element(j,2) == 5,
	         	! then we only calculate density gradient from vertical layer #5 ~ surface layer.
	         	! i.e., there will be no horizontal baroclinic density gradient at bottom two layers at this side (i.e., vertical layer#3 & #4)
               if(adj_cellnum_at_face(2,j) /= 0 .and.   &
            	&  k2 >= bottom_layer_at_element(adj_cellnum_at_face(1,j)) .and.   &
            	&  k2 <= top_layer_at_element(adj_cellnum_at_face(1,j))    .and.   &
            	&  k2 >= bottom_layer_at_element(adj_cellnum_at_face(2,j)) .and.   &
            	&  k2 <= top_layer_at_element(adj_cellnum_at_face(2,j))) then						
						if(k2 == kk) then
							! include a half portion
							sum1 = sum1 + &
							&	0.5_dp*(rho_cell(k2,adj_cellnum_at_face(2,j)) - rho_cell(k2,adj_cellnum_at_face(1,j)))*dz_face(k2,j)
							! sum1 = sum1 + 0.5_dp*dz_face(k2,j)
						else
							! include a full portion
							sum1 = sum1 + &
							&	1.0_dp*(rho_cell(k2,adj_cellnum_at_face(2,j)) - rho_cell(k2,adj_cellnum_at_face(1,j)))*dz_face(k2,j)
							! sum1 = sum1 + dz_face(k2,j)
						end if
               end if
               
               ! tangential baroclinic gradient; i.e., node to node gradient
               ! both nodes must be wet.
               if(k2 >= bottom_layer_at_node(node1)  .and.   &
               &  k2 <= top_layer_at_node(node1) .and.   &
               &  k2 >= bottom_layer_at_node(node2)  .and.   &
               &  k2 <= top_layer_at_node(node2)) then						
						if(k2 == kk) then
							! include a half portion
							sum2 = sum2 + 0.5_dp*(rho_node(k2,node2)-rho_node(k2,node1))*dz_face(k2,j)
							! sum2 = sum2 + 0.5_dp*dz_face(k2,j)
						else
							! include a full portion
							sum2 = sum2 + 1.0_dp*(rho_node(k2,node2)-rho_node(k2,node1))*dz_face(k2,j)
							! sum2 = sum2 + dz_face(k2,j)
						end if
               end if 
            end do ! do k2=kk,top_layer_at_face(j)

				! this was the original equation
 				rtemp1 = dz_face(k,j)*(-gravity*dt/(rho_o*delta_j(j)))*sum1 ! full baroclinic gradient on the fake face-normal axis
 				rtemp2 = dz_face(k,j)*(-gravity*dt/(rho_o*face_length(j)))*sum2 ! full baroclinic gradient on the true face-tangential axis
				
				! this is the new equation for rtemp1 & rtemp2
! 			 	rtemp1 = dz_face(k,j)*(-gravity*dt/(rho_o*delta_j(j)))* &
! 			 	&			sum1*(rho_cell(kk,adj_cellnum_at_face(2,j)) - rho_cell(kk,adj_cellnum_at_face(1,j))) ! full baroclinic gradient on the fake face-normal axis
! 				rtemp2 = dz_face(k,j)*(-gravity*dt/(rho_o*face_length(j)))* &
! 				&			sum2*(rho_node(kk,node2)-rho_node(kk,node1)) ! full baroclinic gradient on the true face-tangential axis

				! this is original...
 				rtemp3 = rtemp1 * cos_theta2(j) ! true face-normal portion from rtemp1
 				rtemp4 = rtemp1 * sin_theta2(j) ! true face-tangential portion from rtemp1
				
				! face-tangential portion which generated from rtemp1 must be added to the true baroclinic tangential gradient, i.e., rtemp2 + rtemp4
 				G1(k) = G1(k) + rtemp3 * spinup_function_baroclinic ! use only true face-normal portion
 				G2(k) = G2(k) + (rtemp2 + rtemp4) * spinup_function_baroclinic
         end if ! baroclinic_flag == 1

         ! include (6)wind stress at the surface -----------------------------!
         if(k == 1) then
         	G1(k) = G1(k) + dt*wind_stress_normal(j)	! wind_stress_normal(j) = gamma_t*u_air in the equation
            G2(k) = G2(k) + dt*wind_stress_tangnt(j)	! wind_stress_tangnt(j) = gamma_t*v_air in the equation            
         end if

			! include (7)earth tidal potential ----------------------------------!
			! detp_dx = (eta(2) - eta(1))/delta_j : horizontal presure gradient due to earth tidal potential
			! here eta(2) and eta(1) are earth tidal elevation, i.e., a*cos(theta) like term.
			! Note: 0.69 is the "effective Earth elasticity factor"; see more details in:
			! 			Numerical Models of Oceans and Oceanic Processes by Kantha and Clayson
			! 			also Schwiderski(1980), Hendershott(1981)
         G1(k) = G1(k) + 0.69_dp*dz_face(kk,j)*gravity*dt*detp_dx
         G2(k) = G2(k) + 0.69_dp*dz_face(kk,j)*gravity*dt*detp_dy

			! include (8)atmosperic pressure gradient ---------------------------!
         ! pressure gradient in normal direction (between elements)
         if(adj_cellnum_at_face(2,j) /= 0) then ! jw, check v56
            do l = 1, 2
               temp(l) = 0.0
               i_temp_element = adj_cellnum_at_face(l,j)
               do ll = 1, tri_or_quad(i_temp_element) ! jw, this is not eleven but alphabet ll
                  temp(l) = temp(l)                                                  &
                    &     + air_p1(nodenum_at_cell(ll,i_temp_element))   &
                    &     / tri_or_quad(i_temp_element)
               end do
            end do
            ! This is for the orthogonal mesh
				G1(k) = G1(k) - dz_face(kk,j)*dt/rho_o*(temp(2)-temp(1))/delta_j(j)
				
				! The following commented lines are for the unorthogonal mesh
! 				rtemp0 = -dt/rho_o*(temp(2)-temp(1))/delta_j(j)
         end if
!         rtemp1 = rtemp0 * cos_theta2(j) ! true face-normal portion
!         rtemp2 = rtemp0 * sin_theta2(j) ! true face-tangential portion
!         G1(k) = G1(k) + dz_face(kk,j) * rtemp1
         
			! pressure gradient in tangential direction
! 			rtemp3 = -dt/rho_o*((air_p1(node2)-air_p1(node1)) / face_length(j))
! 			G2(k) = G2(k) + dz_face(kk,j) * (rtemp3 + rtemp2)
         ! End of forcing term inclusion =====================================!

         
         !====================================================================!
         ! Now, include flux boundary conditions, such as tidal flow (not river inflow), in the momentum equation
         ! Boundary velocities should be asigned directly.
         !====================================================================!

			! impose boundary condition =========================================!
         if(isflowside(j) > 0) then ! actually /= 0; if this side is a tidal flow boundary
				! impose flux b.c.
            ibnd    =  isflowside(j) ! ibnd = ob_segment number
            ! write(*,*) 'j= ',j, ',   isflowside(j) = ', isflowside(j)
            
            ! here u_boundary(ibnd) = u velocity at ibnd's ob segment
            ! Here u_boundary(ibnd)  is u velocity in xy coordinate (not u*v* coordinate), jw, make this sure
            ! That is why we have to convert uv to u*v*
            G1(k) =  u_boundary(ibnd)*cos_theta(j)   &
            &     +  v_boundary(ibnd)*sin_theta(j)
            G2(k) = -u_boundary(ibnd)*sin_theta(j)   &
            &     +  v_boundary(ibnd)*cos_theta(j)
         else if(boundary_type_of_face(j) == -1) then 
				! land b.c.
         	! if this side is a land boundary, normal velocity should be zero. 
            G1(k) = 0.0_dp ! velocity normal should be zero at land boundary
         end if
         
         ! impose river boundary condition ===================================!
         ! We can include river inflow information either in the momentum equation (Eq (44) ~ (45)) 
         ! or in the free-surface equation (which is the continuity equation) (Eq(46)).
         ! Both approach has pros and cons:
         ! 		Mass conservation is better achieved when using the continuity equatioin,
         ! 		while we can get realistic flow vectors at boundary face when using the momentum equation. 
         ! Note: if river is activated here, you also must correct carefully next step shown after "call tridiagonal_solver".
         if(isflowside3(j) > 0) then
 			 	ibnd = isflowside3(j)
 				
 				! river flow acts perpendicular to the face,
 				! thus we only have to update G1 not G2:
            G1(k) = G1(k) + dz_face(kk,j) * Qu_boundary(ibnd)
            !! G2(k) =  G2(k) + dz_face(kk,j) * Qv_boundary(ibnd)
			end if
			
			! impose Withdraw/Return boundary condtion ==========================!
			if(isflowside4(j) > 0) then
				ibnd = isflowside4(j)
				if(WR_layer(ibnd) == 999) then
					! anytime surface layer
					if(kk == top_layer_at_face(j)) then
						G1(k) = G1(k) + dz_face(kk,j) * WRu_boundary(ibnd)
					end if
				else if(WR_layer(ibnd) == 0) then
					! anytime bottom layer
					if(kk == bottom_layer_at_face(j)) then
						! write(*,*) 'j, k, kk, WRu_boundary(ibnd) = ', j, k, kk, WRu_boundary(ibnd)
						G1(k) = G1(k) + dz_face(kk,j) * WRu_boundary(ibnd)
					end if
				else
					if(kk == WR_layer(ibnd)) then
						G1(k) = G1(k) + dz_face(kk,j) * WRu_boundary(ibnd)
					end if
				end if
			end if


			! right hand side 
         rhs(k,1) = G1(k)
         rhs(k,2) = G2(k)
         rhs(k,3) = dz_face(kk,j)
		end do ! k = 1, top_layer_at_face(j)-bottom_layer_at_face(j)+1  (k=1,num_vertical_layer)

		! gam & solution are output results, we don't need gam, which is the pivot value for tridiagonal matrix solver
		! solution(:,1) = AinvG1
		! solution(:,2) = AinvG2
		! solution(:,3) = AinvDZ (AinvDeltaZ1 & AinvDeltaZ2)
		
		! Actually, we just need to give "maxlayer", which is the maximum dimension of the matrix, for the first item of "tridiagonal_solver"
      call tridiagonal_solver(maxlayer+1, num_vertical_layer, 3, a_lower_mat, b_diagonal_mat, c_upper_mat, rhs, solution, gam)
		
		! Caution: the eq. is different at b.c. ================================!
		! Thus, those values should be explicitly given.      
		! Note, following part is extremely important.
		! get AinvG1, AinvG2, AinvDeltaZ1, and AinvDeltaZ2 from Thomas algorithm:
		! 		do k=1,num_vertical_layer
		! 			AinvG1(k,j) = solution(k,1) 		! , normal velocity contribution from explicit terms
		! 			AinvG2(k,j) = solution(k,2) 		! , tangential velocity contribution from explicit terms
		! 			AinvDeltaZ1(k,j) = solution(k,3) ! , normal velocity contribution from implicit barotropic gradient between cell centers
		! 			AinvDeltaZ2(k,j) = solution(k,3) ! , tangential velocity contribution from implicit barotropic gradient between two nodes
		! 		end do
		! ======================================================================!
		! (1) Contribution from explicit forcing terms -------------------------!
 		if(isflowside(j) > 0) then
			! if this side is a tidal flow boundary,
			! use pre-calculated (explicit) values assuming that the tidal flow is dominant.		
 			do k=1,num_vertical_layer
  				AinvG1(k,j) = G1(k)
  				AinvG2(k,j) = G2(k)				
 			end do			
 		else if(boundary_type_of_face(j) == -1 .and. isflowside3(j) == 0 .and. isflowside4(j) == 0) then
			! if this side is a land boundary side and not a river boundary side and not a withdraw/return boundary,
			! normal velocity should be zero, but tangential velocity is not zero at the land boundary.
			! i.e., for a river boundary side, we keep the solution.			
 			do k=1,num_vertical_layer
  				AinvG1(k,j) = 0.0_dp
  				AinvG2(k,j) = solution(k,2)
 			end do
 		else
			! if this side is a normal water cell side,
			! use calculated value.			
 			do k=1,num_vertical_layer
 				AinvG1(k,j) = solution(k,1)
 				AinvG2(k,j) = solution(k,2)
 			end do			
 		end if

		! this is the REM3D_v6:
! 		do k=1,num_vertical_layer
! 			if(isflowside(j) > 0 .or. boundary_type_of_face(j) == -1) then
! 				AinvG1(k,j) = G1(k)
! 			else
! 				AinvG1(k,j) = solution(k,1)
! 			end if
			
! 			if(isflowside(j) > 0) then
! 				AinvG2(k,j) = G2(k)
! 			else
! 				AinvG2(k,j) = solution(k,2)
! 			end if
! 		end if
		
 		! (2) Contribution from implicit barotropic gradient -------------------!
 		if(isflowside(j) > 0) then
			! if this side is a tidal flow boundary,
			! there will be no contribution from implicit barotropic gradient (neither normal nor tangential)	 		 			
 			do k=1,num_vertical_layer
 				AinvDeltaZ1(k,j) = 0.0_dp
  				AinvDeltaZ2(k,j) = 0.0_dp
  			end do
		else if(boundary_type_of_face(j) == -1) then 
		! here, if I use (adj_cellnum_at_face(j,2) == 0), it will also include the open boundary side, 
		! and this process will delete the normal barotropic gradient term from the ghost cell.
		! Thus, I must use "boundary_type_of_face" here
		! Note, I used "adj_cellnum_at_face(j,2) == 0" in without ghost cell approach.
		
			! If this side is a land boundary side,
			! there will be no contribution from implicit normal barotropic gradient,
			! but there will be contribution from implicit tangential barotropic gradient between two nodes
 			do k=1,num_vertical_layer
  				AinvDeltaZ1(k,j) = 0.0_dp
  				AinvDeltaZ2(k,j) = solution(k,3)
 			end do
 		else
			! If this side is a normal water side & open boundary side,
			! there will be contribution from normal and tangential barotropic gradient.			
 			do k=1,num_vertical_layer
 				AinvDeltaZ1(k,j) = solution(k,3)
 	 			AinvDeltaZ2(k,j) = solution(k,3)
 	 		end do			
 		end if

		! this is the REM3D_v6:
! 		do k=1,num_vertical_layer
! 			if(isflowside(j) > 0 .or. adj_cellnum_at_face(2,j) == 0) then
! 				AinvDeltaZ1(k,j) = 0.0_dp
! 			else
! 				AinvDeltaZ1(k,j) = solution(k,3)
! 			end if
			
! 			if(isflowside(j) > 0) then
! 				AinvDeltaZ2(k,j) = 0.0_dp
! 			else
! 				AinvDeltaZ2(k,j) = solution(k,3)
! 			end if
! 		end do


! 		if(dia_momentum == 1) then
! 			if(j == 1) then
! 				write(pw_dia_momentum,*) 'it = ', it, ', elapsed_time = ', elapsed_time		
! 				write(pw_dia_momentum,'(A)') &
! 				&	'j, a_lower_mat(maxlayer+1), b_diagonal_mat(maxlayer+1), c_upper_mat(maxlayer+1), &
! 				&	rhs(maxlayer+1,1), rhs(maxlayer+1,2), rhs(maxlayer+1,3)'
! 			end if
! 	 		write(pw_dia_momentum,'(A, I5, *(E15.5))') &
! 	 		&	'j= ', j, &
! 	 		&	(a_lower_mat(k), k=1,maxlayer), &
! 	 		&	(b_diagonal_mat(k), k=1,maxlayer), &
! 	 		&	(c_upper_mat(k), k=1,maxlayer), &
! 	 		&	(rhs(k,1), k=1,maxlayer), &
! 	 		&	(rhs(k,2), k=1,maxlayer), &
! 	 		&	(rhs(k,3), k=1,maxlayer)
! 		end if

	end do ! do j=1,maxface
	!omp end parallel do

	! write diagonostic file
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
