!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! wave-continuity equations
!! 
!! ob_eta_type(i)	= -1 : radiation bc
!!                       no input in this file; elevations are computed
!!                       as average of surrounding elevations
!!                =  1 : this boundary is forced by tidal harmonic constituents
!!                       need tidal constituent, amplitude and 
!!                       phase angle at each element
!! 					=  2 : time history of elevation on this boundary
!!                       time history of elevation(real tidal data) is read
!!                       from data file fort.250 which is "real_tide.dat"
!! 
subroutine solve_free_surface_equation
   use mod_global_variables
   use mod_file_definition
   implicit none

   integer :: i, j, k, l, kk
   integer :: ii, icount, ie, nd, ibnd

   real(dp):: sum0, sum1, dzT_Ai_dz, const1, const2, vnorm, rel
	
	real(dp),dimension(maxele) 	:: rhs_1, eta_guess
	real(dp),dimension(maxele*5) 	:: a1
	real(dp),dimension(0:4,maxele):: sparsem
	integer,dimension(maxele+1) 	:: ia
	integer,dimension(maxele*5) 	:: ja
	integer :: n
   integer :: nz
   integer :: array_check
   ! End of local variables ==================================================!

	! initialize variables
   rhs_1 	 	= 0.0_dp
   eta_guess 	= 0.0_dp
   a1 		 	= 0.0_dp
   sparsem 		= 0.0_dp
   ia 			= 0
   ja 			= 0

   
   ! =========================================================================!
   ! First, we have to prepare the surface elevation at open boundary elements,
   ! which is externally given, i.e., known value.
   ! this is done at compute_boundary_eta.f90
   ! End of surface open boundary preparation ================================!

	! setup coefficient matrix, sparsem, for the wave Eq(44) ~ (46) ===========!
	! no boundary conditions are involved yet.

   ! The original sparse matrix should be a [maxele, maxele] matrix, and
   ! it should be a five-diagonal matrix since it should include informaion of
   ! the current element, and 4 sides when we use rectangular grids.
   ! However, if we reconstruct the original matrix only including non-zero elements, 
   ! the sparse matrix can be reduced to [maxele,5] since we have maximum 4 sides:
   ! 		sparsem(0,i) = information of ith element itself: (area of ith element + dzT_Ai_dz term), this is the diagonal value and it should be dominant
   ! 		sparsem(1,i) = information of the 1st side's neighbor element
   ! 		sparsem(2,i) = information of the 2nd side's neighbor element
   ! 		sparsem(3,i) = information of the 3rd side's neighbor element
   ! 		sparsem(4,i) = information of the 4th side's neighbor element (only for quadrilateral element)
  	
	!$omp parallel do private(i,j,l,k,ii,dzT_Ai_dz,const1,const2,kk,ibnd,vnorm)
	do i = 1, maxele
      sparsem(0,i) = area(i) ! P_i in Eq(46)
      
      ! calculate left hand side of wave equation ============================!
      ! Eq (46) in Lee et al. or Eq (15) in Casulli & Walters
      do l = 1, tri_or_quad(i) ! we have to update equation for each side of this element 
         ii = adj_cellnum_at_cell(l,i)
         j = facenum_at_cell(l,i)
         
         ! if(top_layer_at_face(j) /= 0 .and. ii /= 0) then ! this is wrong
         if(top_layer_at_face(j) /= 0 .and. boundary_type_of_face(j) /= -1) then ! this is correct
         	! if it is not a dry element, i.e., if this element is a wet element,
         	! and this side is not a land boundary face (i.e., adjacent cell exists or this is an open boundary face),
         	! claculate the contribution of barotropic gradient between elements, i.e., g*dt^2*theta^2*sigma(~~)*(eta2 - eta1) term.
         	! Here, dzT_Ai_dz is the [(deltaZ)_transpose*Ainv*DeltaZ] in the eqaution
            dzT_Ai_dz = 0.0_dp
            
            ! integrate [dzT_Ai_dz] over vertical column in this face
            do k = 1, top_layer_at_face(j) - bottom_layer_at_face(j) + 1
            	! Note, the only face-normal velocity contributes to surface water variation.
            	! That is why only AinvDeltaZ1 is used.
            	! Note: AinvDeltaZ1 is already calculated in solve_momentum_equation_omp.f90
               dzT_Ai_dz = dzT_Ai_dz + dz_face(top_layer_at_face(j)+1-k,j) * AinvDeltaZ1(k,j)
            end do
				
				! update dzT_Ai_dz as g*dt^2*theta^2*lamda/delta*dzT_Ai_dz for this side
            dzT_Ai_dz = gravity * dt**2 * theta**2 * face_length(j) / delta_j(j) * dzT_Ai_dz	

            ! now we have to only include true face-normal component from the fake face-normal value
            ! Dr. Sun's approach
            ! dzT_Ai_dz = dzT_Ai_dz / sin_theta2(j)
            
            ! jw's approach
            ! let's kill this now...
            ! dzT_Ai_dz = dzT_Ai_dz * cos_theta2(j)
            
            ! check if it is positive definite...
!            if(dzT_Ai_dz < 0.0_dp) then
!               write(pw_run_log,*) 'Not positive definite at element number', i, dzT_Ai_dz
!               write(*,*) 'Not positive definite at element number', i, dzT_Ai_dz
!               stop
!            end if
            
				! sparsem(i,0) is summation of dzT_Ai_dz term over all sides plus area(i)
				!              in the left-hand-side of Eq(46) 
				! sparsem(i,l) is summation of dzT_Ai_dz term over the vertical direction
				!              in the left-hand-side of Eq(46) i.e. all 2nd term of Eq(46) 
            sparsem(0,i) = sparsem(0,i) + dzT_Ai_dz 	! Area * eta_j1^(n+1) + g * dt^2 * theta^2 * ~~~ * eta_j1^(n+1) term; for current cell; this should be dominant since it has the area term
            
            if(ii /= 0) then
            	! only if it is the side which sharing with a neighbor cell, 
            	! i.e., not open boundary or land boundary side, include barotropic contribution.
            	! Note that open boundary side information will be added to the right-hand-side vector, rhs. 
            	sparsem(l,i) = -dzT_Ai_dz	! ~~~eta_j2^(n+1) term; for current cell's each neighbour cell; this should be small comparing the main element
            end if
          
          ! not yet correctly implemented for the radiation boundary face...  
!         else if(isflowside2(j) > 0) then ! if this side is the radiation boundary face
!         	! radiation bc
!            dzT_Ai_dz = 0.0_dp
!            do k = bottom_layer_at_face(j), top_layer_at_face(j)
!            	dzT_Ai_dz = dzT_Ai_dz               &
! 					&   			+ dt*face_length(j)   &
! 					&   			* dsqrt(gravity/h_face(j))*dz_face(k,j)
!            end do
!            sparsem(i,0) = sparsem(i,0) + dzT_Ai_dz
         end if ! if(top_layer_at_face(j) /= 0) then
      end do ! l=1,tri_or_quad(i)
      ! End of left hand side of wave Eq(46) =================================!

		! calculate right hand side of wave Eq(46) =============================!
		! Note: I can merge these two (left & right side of equations) section, but I keep this way to easy to understand...
      rhs_1(i) = area(i) * eta_cell(i)
      
      ! Note, the only face-normal velocity contributes to surface water variation.
      do l = 1, tri_or_quad(i)
         j = facenum_at_cell(l,i)
         if(top_layer_at_face(j) /= 0 .and. boundary_type_of_face(j) /= -1) then ! this is correct
            const1 = 0.0_dp
            const2 = 0.0_dp
            do k = 1, top_layer_at_face(j) - bottom_layer_at_face(j) + 1
               kk = top_layer_at_face(j) + 1 - k
               const1 = const1 + dz_face(kk,j)*un_face(kk,j) 	! DeltaZ_transpose * U term in the equation
               const2 = const2 + dz_face(kk,j)*AinvG1(k,j) 		! DeltaZ_transpose * AinvG term in the equation
            end do

            rhs_1(i) = rhs_1(i)                                                      &
            &    	  - (1.0_dp-theta)*dt*sign_in_outflow(l,i)*face_length(j)*const1   &
            &    	  -         theta *dt*sign_in_outflow(l,i)*face_length(j)*const2
         end if
			
			! impose radiation boundary condition (passive boundary =============!
         if(isflowside2(j) > 0) then
            ibnd  = isflowside2(j)
            vnorm = u_boundary(ibnd)*cos_theta(j)   &
				&   	+ v_boundary(ibnd)*sin_theta(j)

            do k = bottom_layer_at_face(j),top_layer_at_face(j)
               rhs_1(i) = rhs_1(i) - dt*face_length(j)*dz_face(k,j)*vnorm
            end do
         end if   ! if(isflowside2(j) > 0) then
         
         ! include river boundary ============================================!
         ! Note: river discharge is now included in "solve_momentum_equation_omp.f90"
         ! so, deactivate the followings
         ! Both will give idential results, but this will use boundary face information, and this will be slightly slower.
         ! You can put more than 1 river boundary faces in one river boundary element if you use this approach.
         ! if(isflowside3(j) > 0) then ! if this side is river boundary side
 			! 	ibnd = isflowside3(j) ! ith river
 			! 	rhs2(i) = rhs2(i) + Q_add(ibnd)*dt
			! end if
      end do   ! l=1,tri_or_quad
      
      ! include river boundary -----------------------------------------------!
      ! Note: river discharge is now included in "solve_momentum_equation_omp.f90"
      ! so, deactivate the followings
      ! Both will give idential results, but this will use boundary element information, and this will be slightly faster.
      ! if(Qb_element_flag(i) /= 0) then ! if this element is set to river boundary
      ! 	ibnd = Qb_element_flag(i) ! ith river
      ! 	rhs(i) = rhs(i) + Q_add(ibnd)*dt
      ! end if
      
      ! include boundary ghost element information if this element is an open boundary element
      if(ob_element_flag(i) > 0) then
      	do l=1,tri_or_quad(i)
      		j = facenum_at_cell(l,i)
      		if(boundary_type_of_face(j) > 0) then
      			if(top_layer_at_face(j) /= 0) then
		         	! if it is not a dry element, i.e., if this element is a wet element; boundary element must be wet all time.
		         	! dzT_Ai_dz is the [(deltaZ)_transpose*Ainv*DeltaZ] in the eqaution
		            dzT_Ai_dz = 0.0_dp 
		            
		            ! integrate [dzT_Ai_dz] over vertical column in this face
		            do k = 1, top_layer_at_face(j) - bottom_layer_at_face(j) + 1
		            	! Note, the only face-normal velocity contributes to surface water variation.
		            	! That is why only AinvDeltaZ1 is used.
		               dzT_Ai_dz = dzT_Ai_dz + dz_face(top_layer_at_face(j)+1-k,j) * AinvDeltaZ1(k,j) 
		            end do
						
						! update dzT_Ai_dz as g*dt^2*theta^2*lamda/delta*dzT_Ai_dz for this side
		            dzT_Ai_dz = gravity * dt**2 * theta**2 * face_length(j) / delta_j(j) * dzT_Ai_dz	
		            
		            ! now, move (dzT_Ai_dz * eta) term to the right hand side since this term is explicitly given at ghost boundary element
		            ! note: eta_at_ob_new() is already calculated at compute_boundary_eta.f90
		            rhs_1(i) = rhs_1(i) + dzT_Ai_dz*eta_at_ob_new(ob_element_flag(i))
		            ! write(*,*) ob_element_flag(i), eta_at_ob_new(ob_element_flag(i))
		            ! write(*,*) dzT_Ai_dz, rhs2(i)
		         end if
      		end if
      	end do
      end if
      ! End of right hand side of wave Eq(46) ================================!
   end do   ! i=1,maxele  
	!$omp end parallel do
	! End of initial setup for the coefficient matrix in the wave Eq(46):
	! 		sparsem: which is the left side matrix
	! 		rhs2: which is the right side of matrix
	! =========================================================================!
   
	! ITPACK2C ================================================================!
	! create a mapping between element index and actual eq. index =============!
	! Assemble sparse matrix format to slove new water surface elevation eta ==!
  	! If an element is assigned as an open boundary (tidal boundary), 
  	! then the element information can be excluded since the boundary condition is
  	! directly (explicitly) applied.
  	! Here, a, ia, and ja are the matrixs used in ITPACK2C; so check ITPACK2C manual for detail. 
!   allocate(a(maxele*5), ia(maxele+1), ja(maxele*5))	! since sparsm(maxele,0:4), and this is the maximum array size before excluding open boundary side information
   
	! The matrix is symmetric, so the new information will only have the upper diagonal elements.
	! n: 	original order of the linear system [n, n] or [maxele, maxele], 
	! nz: total number of non-zero element in the [n, 5] sparse matrix; and in symmetric sparse matrix (our case)
   n = 0		! for the size of IA in ITPACK2C
   nz = 0  	! for the size of A & JA in ITPACK2C
   do i = 1, maxele
		n 		= n +1 						! actual element number (order of the linear system)
      nz 	= nz+1 						! number of non-zero elements
      ia(n) = nz							! position number of the first non-zero value at each row
      ja(nz)= i      					! column number for non-zero values
      a1(nz) = sparsem(0,i) 			! Left-hand-side matrix, including diagonal value
      eta_guess(n)= eta_cell_new(i) ! initial guess

      ! Now, we have to include lower/upper diagonal values 
      ! Note, only this part is different with ITPACK2C version.
      ! ITPACK2C only requires upper diagonal part, but this version requires all.
      do l = 1, tri_or_quad(i)
         ii = adj_cellnum_at_cell(l,i)
         if(ii /= 0) then
				nz = nz+1 ! update non-zero element number
				ja(nz) = ii
				a1(nz) = sparsem(l,i)
         end if
      end do ! l=1,tri_or_quad(i)
   end do ! i=1,maxele
   
   ! position # of the first non-zero value at each row, 
   ! and it requires the position # at the next line (see more details in ITPACK manual)
   ia(n+1) = nz+1 

   ! array check -------------------------------------------------------------!
   ! this is for easy debugging
   ! 0: no check
   ! 1: check
   array_check = 0
   if(array_check == 1) then
	   open(10,file='./check.txt',form='formatted',status='replace')
	   write(10,*) 'area(maxele):'
	   do i=1,maxele
	   	write(10,*) 'area(',i,') = ', area(i)
	   end do
	
	   write(10,*) 'sparsem(maxele,0:4):'
	   do i=1,maxele
	   	write(10,*) 'sparsem(l,',i,') = ', (sparsem(l,i), l=0,4)
	   end do
	
	   write(10,*) 'AinvDeltaZ1(maxlayer,maxface):'
	   do j=1,maxface
	   	write(10,*) 'AinvDeltaZ1(k,',j,') = ', (AinvDeltaZ1(k,j), k=1,maxlayer)
	   end do
	   
		write(10,*)
	   write(10,*) 'n = ', n
	   write(10,*) 'ia = ', size(ia)
	   write(10,*) 'ja = ', size(ja)
	   write(10,*) 'a1 = ', size(a1)
	   write(10,*) 'rhs = ', size(rhs_1)
	   write(10,*) 'eta_guess = ', size(eta_guess)   
	   write(10,*)
	   write(10,*) 'ia = '
	   do i=1,size(ia)
	   	write(10,*) 'ia(',i,') = ', ia(i)
	   end do
	   write(10,*)
	   
	   write(10,*) 'ja = '
	   do i=1,size(ja)
	   	write(10,*) 'ja(',i,') = ', ja(i)
	   end do
	   write(10,*)
	   
	   write(10,*)	'a1 = '
	   do i=1,size(a1)
	   	write(10,*) 'a1(',i,') = ', a1(i)
	   end do
	   write(10,*)
	   
	   write(10,*)	'rhs= '
	   do i=1,size(rhs_1)
	   	write(10,*) 'rhs(',i,') = ', rhs_1(i)
	   end do
	   write(10,*)
	   
	   write(10,*) 'eta_guess = '
	   do i=1,size(eta_guess)
	   	write(10,*) 'eta_guess(',i,') = ', eta_guess(i)
	   end do
	   close(10)
	   stop 'dddd'
	end if

   ! End of sparse matrix preparation ========================================!
   ! =========================================================================!
	
	! Following part is different from ITPACK2C version =======================!   
	! call jw's cg algorithm
	call pre_conj_grad(n,nz,a1,ja,ia,rhs_1,eta_guess) ! Jacobi preconditioner version
	
	! =========================================================================!
	!                       re-assemble new elevations
	! =========================================================================!
	!$omp parallel
	!$omp do private(i)
   do i = 1, maxele
		eta_cell_new(i) = eta_guess(i)
   end do
	!$omp end do

	! impose sponge layer boundary condition ==================================!
	! jw, currently not activated
	if(i_sponge_layer_flag /= 0) then
		!$omp do private(i,l, rel)
		do i = 1, maxele
      	rel = 0.0_dp
			do l = 1, tri_or_quad(i)
				rel = rel + sponge_relax(nodenum_at_cell(l,i))/tri_or_quad(i)
			end do
			eta_cell_new(i) = eta_cell_new(i)*rel + etaic(i)*(1.0_dp-rel)
		end do
		!$omp end do		
	end if

	! compute nodal elevations for tangential velocity ========================!
	! nodal elevation will be used for tangential velocity calculation
	! We cannot simply calculate as following (e.g.) since each element has a different area: 
	! 		eta_node(i) = (eta1 + eta2 + eta3)/3 (if a node has three adjacent elements), we should not use this method 
	! That is why we have to calculate as:
	! 		eta_node(i) = (total volume of adjacent elements)/(total area of adjacent elements)
	!$omp do private(i,j, ie,sum1,sum0)
   do i = 1, maxnod
      wetdry_node(i) = 0   ! reliable elevation

      sum1  = 0.0_dp 	! sum0 of elemental elevations
      sum0  = 0.0_dp 	! sum0 of area

      do j = 1, adj_cells_at_node(i)
         ie = adj_cellnum_at_node(j,i)
         if(top_layer_at_element(ie) == 0) then
            wetdry_node(i) = 1 ! if element is dry, set wetdry_node to 1, otherwise, wetdry_node = 0 (for all non-dry elements)
         else
            sum1 = sum1 + area(ie)*eta_cell_new(ie) 	! total_volume
            sum0  = sum0  + area(ie)						! total_area
         end if
      end do
		
      if(wetdry_node(i) == 0) then
      	! if this is a wet node:
      	eta_node(i) = sum1/sum0 ! total_volume / total_area
      end if
   end do   !  i=1,maxnod
	!$omp end do

	! estimate bed elevations for rewetting ===================================!
	!$omp do private(i,j, nd,sum0,icount)
   do i = 1, maxnod
      if(wetdry_node(i) == 1) then ! if this is a dry node
         sum0 = 0.0_dp
         icount = 0
         do j = 1, adj_nodes_at_node(i)
            nd = adj_nodenum_at_node(j,i)
            if(wetdry_node(nd) == 0) then ! if the neighbor node is a wet node
               icount = icount+1
               sum0 = sum0 + eta_node(nd)
            end if
         end do
         if(icount /= 0) then
         	eta_node(i) = sum0 / icount
         end if
      end if
   end do
	!$omp end do
	!$omp end parallel
	! end of calculating water surface elevation at all nodes =================!	

	! write diagnostic file
	if(dia_freesurface == 1) then
		write(pw_dia_freesurface,'(A,I5,A,E15.5)') 'it = ', it, ', elapsed_time = ', elapsed_time
		write(pw_dia_freesurface,'(A,I5)') 'n  = ', n
		write(pw_dia_freesurface,'(A,I5)') 'nz = ', nz
		write(pw_dia_freesurface,*)
		
		write(pw_dia_freesurface,'(A)') 'a(maxele*5), ja(maxele*5):'
		do i=1,maxele*5
			write(pw_dia_freesurface,'(E15.5, I5)') a1(i), ja(i)
		end do
		write(pw_dia_freesurface,*)
		
		write(pw_dia_freesurface,'(A)') 'ia(maxele+1):'
		do i=1,maxele+1
			write(pw_dia_freesurface,'(I5)') ia(i)
		end do
		write(pw_dia_freesurface,*)
		
		write(pw_dia_freesurface,'(A)') 'rhs(maxele), eta_guess(maxele), eta_cell_new(i):'
		do i=1,maxele
			write(pw_dia_freesurface,'(3E15.5)') rhs_1(i), eta_guess(i), eta_cell_new(i)
		end do
		
		write(pw_dia_freesurface,*)
		write(pw_dia_freesurface,'(A)') 'i, eta_node(maxnod)'
		do i=1,maxnod
			write(pw_dia_freesurface,'(I5,E15.5)') i, eta_node(i)
		end do		
	end if

	! deallocate temporary array for matrix solver 
!   deallocate(sparsem)
!   deallocate(rhs)
!   deallocate(ia, ja)
!   deallocate(a)
!   deallocate(eta_guess)
end subroutine solve_free_surface_equation
