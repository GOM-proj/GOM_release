!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! wave-continuity equations
!! 
!! ob_eta_type(i)	= -1 : radiation bc
!!                       no input in this file; elevations are computed
!!                       as average of surrounding elevations
!! 					=  1 : time history of elevation on this boundary
!!                       time history of elevation(real tidal data) is read
!!                       from data file fort.250 which is "real_tide.dat"
!!                =  2 : this boundary is forced by tidal harmonic constituents
!!                       need tidal constituent, amplitude and 
!!                       phase angle at each element
!!                =  3 : Same as #2 but with simaple sine wave
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
   ! jw

	! jw
   rhs_1 	 	= 0.0_dp
   eta_guess 	= 0.0_dp
   a1 		 	= 0.0_dp
   sparsem 		= 0.0_dp
   ia 			= 0
   ja 			= 0

   
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
  	
	!$omp parallel do private(i,j,l,k,ii,dzT_Ai_dz,const1,const2,kk,ibnd,vnorm)
	do i = 1, maxele
      sparsem(0,i) = area(i) ! jw
      
      ! jw
      ! jw
      do l = 1, tri_or_quad(i) ! jw
         ii = adj_cellnum_at_cell(l,i)
         j = facenum_at_cell(l,i)
         
         if(top_layer_at_face(j) /= 0) then
         	! jw
         	! jw
         	! jw
            dzT_Ai_dz = 0.0_dp
            
            ! jw
            do k = 1, top_layer_at_face(j) - bottom_layer_at_face(j) + 1
            	! jw
            	! jw
            	! jw
               dzT_Ai_dz = dzT_Ai_dz + dz_face(top_layer_at_face(j)+1-k,j) * AinvDeltaZ1(k,j)
            end do
				
				! jw
            dzT_Ai_dz = gravity * dt**2 * theta**2 * face_length(j) / delta_j(j) * dzT_Ai_dz	

            ! jw
            ! jw
            ! jw
            
            ! jw
            dzT_Ai_dz = dzT_Ai_dz * cos_theta2(j)
            
            ! jw
            if(dzT_Ai_dz < 0.0_dp) then
               write(pw_run_log,*) 'Not positive definite at element number', i, dzT_Ai_dz
               write(*,*) 'Not positive definite at element number', i, dzT_Ai_dz
               stop
            end if
            
				! jw
				! jw
				! jw
				! jw
            sparsem(0,i) = sparsem(0,i) + dzT_Ai_dz 	! jw
            
            if(ii /= 0) then
            	! jw
            	! jw
            	! jw
            	sparsem(l,i) = -dzT_Ai_dz	! jw
            end if
          
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
         end if ! jw
      end do ! jw
      ! jw

		! jw
      rhs_1(i) = area(i) * eta_cell(i)
      
      ! jw
      do l = 1, tri_or_quad(i)
         j = facenum_at_cell(l,i)
         if(top_layer_at_face(j) /= 0) then ! jw
            const1 = 0.0_dp
            const2 = 0.0_dp
            do k = 1, top_layer_at_face(j) - bottom_layer_at_face(j) + 1
               kk = top_layer_at_face(j) + 1 - k
               const1 = const1 + dz_face(kk,j)*un_face(kk,j) 	! jw
               const2 = const2 + dz_face(kk,j)*AinvG1(k,j) 		! jw
            end do

            rhs_1(i) = rhs_1(i)                                                     &
            &    	  - (1.0-theta)*dt*sign_in_outflow(l,i)*face_length(j)*const1   &
            &    	  -      theta *dt*sign_in_outflow(l,i)*face_length(j)*const2
         end if
			
			! jw
         if(isflowside2(j) > 0) then
            ibnd  = isflowside2(j)
            vnorm = u_boundary(ibnd)*cos_theta(j)   &
				&   	+ v_boundary(ibnd)*sin_theta(j)

            do k = bottom_layer_at_face(j),top_layer_at_face(j)
               rhs_1(i) = rhs_1(i) - dt*face_length(j)*dz_face(k,j)*vnorm
            end do
         end if   ! jw
         
         ! jw
         ! jw
         ! jw
         ! jw
         ! jw
         ! jw
 			! jw
 			! jw
			! jw
      end do   ! jw
      
      ! jw
      ! jw
      ! jw
      ! jw
      ! jw
      ! jw
      ! jw
      ! jw
      
      ! jw
      if(ob_element_flag(i) > 0) then
      	do l=1,tri_or_quad(i)
      		j = facenum_at_cell(l,i)
      		if(boundary_type_of_face(j) > 0) then
      			if(top_layer_at_face(j) /= 0) then
		         	! jw
		         	! jw
		            dzT_Ai_dz = 0.0_dp 
		            
		            ! jw
		            do k = 1, top_layer_at_face(j) - bottom_layer_at_face(j) + 1
		            	! jw
		            	! jw
		               dzT_Ai_dz = dzT_Ai_dz + dz_face(top_layer_at_face(j)+1-k,j) * AinvDeltaZ1(k,j) 
		            end do
						
						! jw
		            dzT_Ai_dz = gravity * dt**2 * theta**2 * face_length(j) / delta_j(j) * dzT_Ai_dz	
		            
		            ! jw
		            ! jw
		            rhs_1(i) = rhs_1(i) + dzT_Ai_dz*eta_at_ob_new(ob_element_flag(i))
		            ! jw
		            ! jw
		         end if
      		end if
      	end do
      end if
      ! jw
   end do   ! jw
	!$omp end parallel do
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
   n = 0		! jw
   nz = 0  	! jw
   do i = 1, maxele
		n 		= n +1 						! jw
      nz 	= nz+1 						! jw
      ia(n) = nz							! jw
      ja(nz)= i      					! jw
      a1(nz) = sparsem(0,i) 			! jw
      eta_guess(n)= eta_cell_new(i) ! jw

      ! jw
      ! jw
      ! jw
      do l = 1, tri_or_quad(i)
         ii = adj_cellnum_at_cell(l,i)
         if(ii /= 0) then
				nz = nz+1 ! jw
				ja(nz) = ii
				a1(nz) = sparsem(l,i)
         end if
      end do ! jw
   end do ! jw
   
   ! jw
   ! jw
   ia(n+1) = nz+1 

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
	call pre_conj_grad(n,nz,a1,ja,ia,rhs_1,eta_guess) ! jw
	
	! jw
	! jw
	! jw
	!$omp parallel
	!$omp do private(i)
   do i = 1, maxele
		eta_cell_new(i) = eta_guess(i)
   end do
	!$omp end do

	! jw
	! jw
	if(i_sponge_layer_flag /= 0) then
		!$omp do private(i,l, rel)
		do i = 1, maxele
      	rel = 0.0_dp
			do l = 1, tri_or_quad(i)
				rel = rel + sponge_relax(nodenum_at_cell(l,i))/tri_or_quad(i)
			end do
			eta_cell_new(i) = eta_cell_new(i)*rel + etaic(i)*(1.0d0-rel)
		end do
		!$omp end do		
	end if

	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	!$omp do private(i,j, ie,sum1,sum0)
   do i = 1, maxnod
      wetdry_node(i) = 0   ! jw

      sum1  = 0.0d0 	! jw
      sum0  = 0.0d0 	! jw

      do j = 1, adj_cells_at_node(i)
         ie = adj_cellnum_at_node(j,i)
         if(top_layer_at_element(ie) == 0) then
            wetdry_node(i) = 1 ! jw
         else
            sum1 = sum1 + area(ie)*eta_cell_new(ie) 	! jw
            sum0  = sum0  + area(ie)						! jw
         end if
      end do

      if(wetdry_node(i) == 0) then
      	! jw
      	eta_node(i) = sum1/sum0 ! jw
      end if
   end do   ! jw
	!$omp end do

	! jw
	!$omp do private(i,j, nd,sum0,icount)
   do i = 1, maxnod
      if(wetdry_node(i) == 1) then ! jw
         sum0 = 0.0d0
         icount = 0
         do j = 1, adj_nodes_at_node(i)
            nd = adj_nodenum_at_node(j,i)
            if(wetdry_node(nd) == 0) then
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
	! jw

	! jw
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

	! jw
! jw
! jw
! jw
! jw
! jw
end subroutine solve_free_surface_equation
