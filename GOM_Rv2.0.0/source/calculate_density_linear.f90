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
   ! jw
   ! jw

	! jw
	! jw
	temp_min = -5.0_dp	! jw
	temp_max = 100.0_dp	! jw
	salt_min = -5.0_dp	! jw
	salt_max = 100.0_dp	! jw

	! jw
	!$omp parallel
	!$omp do private(i,k,rtemp,rsalt)
	do i = 1,maxnod
		! jw
		if(top_layer_at_node(i) /= 0)then
			do k = bottom_layer_at_node(i), top_layer_at_node(i)
            rtemp = temp_node(k,i)   ! jw
            rsalt = salt_node(k,i)   ! jw

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
				
				! jw
				! jw
				! jw
				! jw
				rho_node(k,i) = 1000.0 + 0.7*rsalt
				! jw
				
				if(rho_node(k,i) < 980.0) then
					write(pw_run_log,*) 'Water density is too low (weird density) at node:'
					write(pw_run_log,*) 'it, i, k, temp(i,k), salt(i,k), rho_node(k,i) = ', it, i, k, rtemp, rsalt, rho_node(k,i)
					stop
				end if
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
	!$omp do private(j,k,rtemp,rsalt)
	do j = 1, maxface
		! jw
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
            
				! jw
				! jw
				! jw
				! jw
				! jw
				if(i_density_flag == 1) then
					! jw
					rho_face(k,j) = 1000.0 + 0.7*rsalt
				end if
				! jw
				
				if(rho_face(k,j) < 980.0) then
					write(pw_run_log,*) 'Water density is too low (weird density) at face:'
					write(pw_run_log,*) 'it, j, k, temp(j,k), salt(j,k), rho_face(k,j) = ', it, j, k, rtemp, rsalt, rho_face(k,j)
					stop
				end if
			end do ! jw

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
	!$omp do private(i,k,icount)
	do i = 1, maxele
		if(top_layer_at_element(i) == 0) then
		 	cycle
		end if
		do k = bottom_layer_at_element(i), top_layer_at_element(i) ! jw
			rho_cell(k,i) = 0.0
			icount = 0
			do l = 1, tri_or_quad(i)
				nd = nodenum_at_cell(l,i)
				if(top_layer_at_node(nd) /= 0) then
				   icount = icount + 1
				   rho_cell(k,i) = rho_cell(k,i) + rho_node(k,nd) ! jw
				end if
			end do
			
			if(icount == 0) then
				write(pw_run_log,*) 'There is a wet element with dry nodes at cell #:', i
				stop
			else
				rho_cell(k,i) = rho_cell(k,i) / icount
			end if
		end do ! jw
		
		! jw
		do k = 1, bottom_layer_at_element(i)-1
			rho_cell(k,i) = rho_cell(bottom_layer_at_element(i),i)
		end do
		do k = top_layer_at_element(i)+1, maxlayer
			rho_cell(k,i) = rho_cell(top_layer_at_element(i),i)
		end do
	end do ! jw
	!$omp end do
	!$omp end parallel

	! jw
	
	! jw
	! jw
end subroutine calculate_density_linear

