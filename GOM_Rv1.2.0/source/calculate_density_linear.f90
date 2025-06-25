!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Solve density using temperature and salinity at node, face, and cell from a simple linear equation
!! see Wang et al., 2011 (Modeling and understanding turbulent mixing in a macrotidal salt wedge estuary, section 2.1)
!! This subroutine will calculate:
!! 		rho_node(k,i), rho_face(k,j), and rho_cell(k,i)
!! 
subroutine calculate_density_linear
   use mod_global_variables
   use mod_file_definition
   
   implicit none
   integer :: i, j, k, l, icount, nd
   real(dp):: rtemp, rsalt
   real(dp):: temp_min, temp_max, salt_min, salt_max   
   ! logical,save :: first_call = .true.
   ! End of local variables ==================================================!

	! before start, let's define min/max values:
	! Note, here, I am giving a little more room for each values to avoid termination by small error.
	temp_min = -5.0_dp	! 0.0
	temp_max = 100.0_dp	! 40.0
	salt_min = -5.0_dp	! 0.0
	salt_max = 100.0_dp	! 42.0

	! (1) calculate density at nodes ==========================================!
	!$omp parallel
	!$omp do private(i,k,rtemp,rsalt)
	do i = 1,maxnod
		! calculate density at wet nodes only, dry node not calculated
		if(top_layer_at_node(i) /= 0)then
			do k = bottom_layer_at_node(i), top_layer_at_node(i)
            rtemp = temp_node(k,i)   ! temperature at node
            rsalt = salt_node(k,i)   ! salinity at node

            if(rtemp < temp_min .or. rtemp > temp_max) then
            	write(pw_run_log,*) 'Temperature is out of range: 0.0 < temp < 40.0'
            	write(pw_run_log,*) 'it, i, k, temp(i,k) =', it, i, k, rtemp
            	stop
            end if

				if(rsalt < salt_min .or. rsalt > salt_max) then
            	write(pw_run_log,*) 'Salinity is out of range: 0.0 < salt < 42.0'
            	write(pw_run_log,*) 'it, i, k, salt(i,k) =', it, i, k, rsalt
            	stop
            end if
				
				! option for a simple linear equation ============================!
				! rho = rho_o*(1 + beta*S), where beta = 7.0e-4
				! 		see Wang et al., 2011 (Modeling and understanding turbulent mixing in a macrotidal salt wedge estuary, section 2.1)
				! rho_node(k,i) = 1000.0 + 0.8*rsalt
				rho_node(k,i) = 1000.0 + 0.7*rsalt
				! ================================================================!
				
				if(rho_node(k,i) < 980.0) then
					write(pw_run_log,*) 'Water density is too low (weird density) at node:'
					write(pw_run_log,*) 'it, i, k, temp(i,k), salt(i,k), rho_node(k,i) = ', it, i, k, rtemp, rsalt, rho_node(k,i)
					stop
				end if
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
	!$omp do private(j,k,rtemp,rsalt)
	do j = 1, maxface
		! dry sides will not be updated
		if(top_layer_at_face(j) /= 0)then
			do k = bottom_layer_at_face(j), top_layer_at_face(j)
				rtemp = temp_face(k,j)
				rsalt = salt_face(k,j)

            if(rtemp < temp_min .or. rtemp > temp_max) then
            	write(pw_run_log,*) 'Temperature is out of range: 0.0 < temp < 40.0'
            	write(pw_run_log,*) 'it, j, k, temp(j,k) =', it, j, k, rtemp
            	stop
            end if

				if(rsalt < salt_min .or. rsalt > salt_max) then
            	write(pw_run_log,*) 'Salinity is out of range: 0.0 < salt < 42.0'
            	write(pw_run_log,*) 'it, j, k, salt(j,k) =', it, j, k, rsalt
            	stop
            end if
            
				! ================================================================!
				! option for a simple linear equation
				! rho = rho_o*(1 + beta*S), where beta = 7.0e-4
				! 		see Wang et al., 2011 (Modeling and understanding turbulent mixing in a macrotidal salt wedge estuary, section 2.1)
				! jw, not yet specified, and it should not be here.
				if(i_density_flag == 1) then
					! rho_face(k,j) = 1000.0 + 0.8*rsalt
					rho_face(k,j) = 1000.0 + 0.7*rsalt
				end if
				! ================================================================!
				
				if(rho_face(k,j) < 980.0) then
					write(pw_run_log,*) 'Water density is too low (weird density) at face:'
					write(pw_run_log,*) 'it, j, k, temp(j,k), salt(j,k), rho_face(k,j) = ', it, j, k, rtemp, rsalt, rho_face(k,j)
					stop
				end if
			end do ! k=bottom_layer_at_node(i),top_layer_at_node(i)

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
	!$omp do private(i,k,icount)
	do i = 1, maxele
		if(top_layer_at_element(i) == 0) then
		 	cycle
		end if
		do k = bottom_layer_at_element(i), top_layer_at_element(i) ! wet elements
			rho_cell(k,i) = 0.0
			icount = 0
			do l = 1, tri_or_quad(i)
				nd = nodenum_at_cell(l,i)
				if(top_layer_at_node(nd) /= 0) then
				   icount = icount + 1
				   rho_cell(k,i) = rho_cell(k,i) + rho_node(k,nd) ! rho_node extended
				end if
			end do
			
			if(icount == 0) then
				write(pw_run_log,*) 'There is a wet element with dry nodes at cell #:', i
				stop
			else
				rho_cell(k,i) = rho_cell(k,i) / icount
			end if
		end do ! k=bottom_layer_at_element(i),top_layer_at_element(i)
		
		! extend to account for level changes
		do k = 1, bottom_layer_at_element(i)-1
			rho_cell(k,i) = rho_cell(bottom_layer_at_element(i),i)
		end do
		do k = top_layer_at_element(i)+1, maxlayer
			rho_cell(k,i) = rho_cell(top_layer_at_element(i),i)
		end do
	end do ! i=1,maxele
	!$omp end do
	!$omp end parallel

	! first_call=.false.
	
	! write(*,'(A, I4, 10F10.5)') 'it, rho_cell(1,i) = ', it, (rho_cell(1,i), i=46,55)
	! write(*,*) 'it, rho_cell(1,i) = ', it, (rho_cell(1,i), i=46,55)
end subroutine calculate_density_linear

