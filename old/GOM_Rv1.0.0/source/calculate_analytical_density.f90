!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine calculate_analytical_density
	use mod_global_variables
	use mod_file_definition
	implicit none

   integer :: i, j, k
   ! jw
   ! jw
	! jw
	
	! jw
	rho_node = rho_o
	rho_face = rho_o
	rho_cell = rho_o
	
	! jw
	!$omp parallel
	!$omp do private(i,k)
	do i = 1,maxnod
		! jw
		if(top_layer_at_node(i) /= 0)then
			do k = bottom_layer_at_node(i), top_layer_at_node(i)
				! jw
				! jw
            rho_node(k,i) = rho_o + salt_node(k,i)
            ! jw
			end do ! jw

			! jw
			do k = 1, bottom_layer_at_node(i)-1
			   rho_node(k,i) = rho_node(bottom_layer_at_node(i),i)
			end do
			do k = top_layer_at_node(i)+1, maxlayer
			   rho_node(k,i) = rho_node(top_layer_at_node(i),i)
			end do
		end if ! jw
	end do   ! jw
	!$omp end do nowait

	! jw
	!$omp do private(j,k)
	do j = 1, maxface
		! jw
		if(top_layer_at_face(j) /= 0)then
			do k = bottom_layer_at_face(j), top_layer_at_face(j)
				! jw
				! jw
				rho_face(k,j) = (rho_o + rho_a) + salt_face(k,j)
				! jw
			end do

			! jw
			do k = 1, bottom_layer_at_face(j)-1
			   rho_face(k,j) = rho_face(bottom_layer_at_face(j),j)
			end do
			do k = top_layer_at_face(j)+1, maxlayer
			   rho_face(k,j) = rho_face(top_layer_at_face(j),j)
			end do
		end if ! jw
	end do ! jw
	!$omp end do
	
	! jw
	!$omp do private(i,k)
	do i = 1, maxele
		if(top_layer_at_element(i) /= 0) then
			do k = bottom_layer_at_element(i), top_layer_at_element(i) ! jw
				rho_cell(k,i) = (rho_o + rho_a) + salt_cell_new(k,i)
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
			end do ! jw
			
			! jw
			do k = 1, bottom_layer_at_element(i)-1
				rho_cell(k,i) = rho_cell(bottom_layer_at_element(i),i)
			end do
			do k = top_layer_at_element(i)+1, maxlayer
				rho_cell(k,i) = rho_cell(top_layer_at_element(i),i)
			end do
		end if ! jw
	end do ! jw
	!$omp end do
	!$omp end parallel
	
	
	
! jw
! jw
! jw
! jw

end subroutine calculate_analytical_density