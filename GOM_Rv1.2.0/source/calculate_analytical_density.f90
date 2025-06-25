!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine calculate_analytical_density
	use mod_global_variables
	use mod_file_definition
	implicit none

   integer :: i, j, k
   ! integer :: l, nd, icount
   ! real(dp):: sum1
	! End of local variables ==================================================!
	
	! Initialize rho_node, rho_face, and rho_cell since they are initially initialized to 0.0 in allocate_variables.f90
	rho_node = rho_o
	rho_face = rho_o
	rho_cell = rho_o
	
	! (1) calculate density at nodes ==========================================!
	!$omp parallel
	!$omp do private(i,k)
	do i = 1,maxnod
		! calculate density at wet nodes only, dry node not calculated
		if(top_layer_at_node(i) /= 0)then
			do k = bottom_layer_at_node(i), top_layer_at_node(i)
				! a simple fake density equation 
				! rho_node(k,i) = 1000.0 + 0.8*rsalt
            rho_node(k,i) = rho_o + salt_node(k,i)
            ! rho_node(k,i) = 1000.0 + salt_node(k,i) * 10.0 ! 10 times
			end do ! k=bottom_layer_at_node(i),top_layer_at_node(i)

			! extend to account for level changes
			do k = 1, bottom_layer_at_node(i)-1
			   rho_node(k,i) = rho_node(bottom_layer_at_node(i),i)
			end do
			do k = top_layer_at_node(i)+1, maxlayer
			   rho_node(k,i) = rho_node(top_layer_at_node(i),i)
			end do
		end if ! if(top_layer_at_node(i) /= 0)then
	end do   ! i = 1, maxnod
	!$omp end do nowait

	! (2) calculate density at face ===========================================!
	!$omp do private(j,k)
	do j = 1, maxface
		! dry sides will not be updated
		if(top_layer_at_face(j) /= 0)then
			do k = bottom_layer_at_face(j), top_layer_at_face(j)
				! a simple fake density equation 
				! rho_face(k,j) = 1000.0 + 0.8*rsalt
				rho_face(k,j) = (rho_o + rho_a) + salt_face(k,j)
				! rho_face(k,j) = 1000.0 + salt_face(k,j) * 10.0 ! 10 times
			end do

			! extend to account for level changes
			do k = 1, bottom_layer_at_face(j)-1
			   rho_face(k,j) = rho_face(bottom_layer_at_face(j),j)
			end do
			do k = top_layer_at_face(j)+1, maxlayer
			   rho_face(k,j) = rho_face(top_layer_at_face(j),j)
			end do
		end if ! if(top_layer_at_face(j) /= 0)then
	end do ! j=1, maxface
	!$omp end do
	
	! (3) calculate density at cell ===========================================!
	!$omp do private(i,k)
	do i = 1, maxele
		if(top_layer_at_element(i) /= 0) then
			do k = bottom_layer_at_element(i), top_layer_at_element(i) ! wet elements
				rho_cell(k,i) = (rho_o + rho_a) + salt_cell_new(k,i)
				! rho_cell(k,i) = 1000.0 + salt_cell_new(k,i) * 10.0 ! 10 times
				
! 				sum1 = 0.0
! 				icount = 0
! 				do l = 1, tri_or_quad(i)
! 					nd = nodenum_at_cell(l,i)
! 					if(top_layer_at_node(nd) /= 0) then
! 					   icount = icount + 1
! 					   sum1 = sum1 + rho_node(k,nd) ! rho_node extended
! 					end if
! 				end do
				
! 				if(icount == 0) then
! 					write(pw_run_log,*) 'There is a wet element with dry nodes at cell #:', i
! 					stop
! 				else
! 					rho_cell(k,i) = sum1 / icount
! 				end if
			end do ! k=bottom_layer_at_element(i),top_layer_at_element(i)
			
			! extend to account for level changes
			do k = 1, bottom_layer_at_element(i)-1
				rho_cell(k,i) = rho_cell(bottom_layer_at_element(i),i)
			end do
			do k = top_layer_at_element(i)+1, maxlayer
				rho_cell(k,i) = rho_cell(top_layer_at_element(i),i)
			end do
		end if ! if(top_layer_at_element(i) /= 0) then
	end do ! i=1,maxele
	!$omp end do
	!$omp end parallel
	
	
	
! 	do i=1,maxele
! 		write(333,*) rho_cell(1,i)
! 	end do
! 	stop 'I am here'

end subroutine calculate_analytical_density