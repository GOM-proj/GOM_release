!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Solve density using temperature and salinity at node, face, and cell from pond and pickard's book
!! P310, Pond and Pickard's 2nd edition, "A.3.2 International Equation of State of Sea Water, 1980 
!! validity region: temperature: [0,40], salinity: [0:42]
!! Note: this equation is the identical equation which used in EFDC:
!! 		see more details i Ji's book p15 (Hydrodynamics and Water Quality)
!! This subroutine will calculate:
!! 		rho_node(k,i), rho_face(k,j), and rho_cell(k,i)
!! 
subroutine calculate_density_full  
   use mod_global_variables
   use mod_file_definition
   
   implicit none
   integer :: i, j, k, l, kk, icount, nd
   real(dp):: rtemp, rsalt, pres, ptmp, sbm, ratio
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
	!$omp do private(i,k,kk,rtemp,rsalt,pres,ptmp,sbm,ratio)
	do i = 1,maxnod
		! jw
		if(top_layer_at_node(i) /= 0)then
			do k = bottom_layer_at_node(i), top_layer_at_node(i)
            rtemp = temp_node(k,i)   ! jw
            rsalt = salt_node(k,i)   ! jw
				
				! jw
				! jw
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
				rho_node(k,i) = 1000.0 - 0.157406 					&
				&								+ 6.793952e-2*rtemp 		&
				&								- 9.095290e-3*rtemp**2	&
				&								+ 1.001685e-4*rtemp**3	&
				&								- 1.120083e-6*rtemp**4	&
				&								+ 6.536332e-9*rtemp**5	&
				& 					+ rsalt * (0.824493 		&
				&				 				- 4.089900e-3*rtemp		&
				&								+ 7.643800e-5*rtemp**2	&
				&								- 8.246700e-7*rtemp**3	&
				&								+ 5.387500e-9*rtemp**4)	&
				&					+ dsqrt(rsalt)**3 * (-5.72466e-3 	&
				&								+ 1.022700e-4*rtemp		&
				&								- 1.654600e-6*rtemp**2)	&
				&					+ 4.831400e-4*rsalt**2
				
				! jw
				! jw
				if(MSL-z_level(k) > 100.0) then
					! jw
					pres = 0.0 ! jw
					do kk = k, top_layer_at_node(i)
						ptmp = 9.81*rho_node(kk,i)*dz_node(kk,i) ! jw
						if(kk == k) then
							! jw
							pres = pres + 1.0e-5 * 0.5*ptmp ! jw
						else
							pres = pres + 1.0e-5 * ptmp
						end if
					end do

					! jw
					! jw
					sbm = 19652.21 + 148.4206*rtemp - 2.327105*rtemp**2 + 1.360477e-2*rtemp**3 - 5.155288e-5*rtemp**4 &
					&  + pres * (3.239908 + 1.43713e-3*rtemp + 1.16092e-4*rtemp**2 - 5.77905e-7*rtemp**3) &
					&  + pres**2 * (8.50935e-5 - 6.12293e-6*rtemp + 5.2787e-8*rtemp**2)                  &
					&  + rsalt * (54.6746 - 0.603459*rtemp + 1.09987e-2*rtemp**2 - 6.1670e-5*rtemp**3)     &
					&  + dsqrt(rsalt)**3 * (7.944e-2 + 1.6483e-2*rtemp - 5.3009e-4*rtemp**2)              &
					&  + pres * rsalt * (2.2838e-3 - 1.0981e-5*rtemp - 1.6078e-6*rtemp**2)                &
					&  + 1.91075e-4*pres*dsqrt(rsalt)**3 - 9.9348e-7*pres**2*rsalt                       &
					&  + 2.08160e-8*rtemp*pres**2*rsalt + 9.1697e-10*rtemp**2*pres**2*rsalt
					
					ratio = 1.0 / (1.0 - pres/sbm)
					rho_node(k,i) = rho_node(k,i) * ratio ! jw
				end if ! jw

				if(rho_node(k,i) < 980.0) then
					write(pw_run_log,*) 'Water density is too low (weird density) at node:'
					write(pw_run_log,*) 'i, k, temp(i,k), salt(i,k), rho_node(k,i) = ', i, k, rtemp, rsalt, rho_node(k,i)
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
	!$omp do private(j,k,kk,rtemp,rsalt,pres,ptmp,sbm,ratio)
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
				rho_face(k,j) = 1000.0 - 0.157406 					&
				&								+ 6.793952e-2*rtemp 		&
				&								- 9.095290e-3*rtemp**2	&
				&								+ 1.001685e-4*rtemp**3	&
				&								- 1.120083e-6*rtemp**4	&
				&								+ 6.536332e-9*rtemp**5	&
				& 					+ rsalt * (0.824493 		&
				&				 				- 4.089900e-3*rtemp		&
				&								+ 7.643800e-5*rtemp**2	&
				&								- 8.246700e-7*rtemp**3	&
				&								+ 5.387500e-9*rtemp**4)	&
				&					+ dsqrt(rsalt)**3 * (-5.72466e-3 	&
				&								+ 1.022700e-4*rtemp		&
				&								- 1.654600e-6*rtemp**2)	&
				&					+ 4.831400e-4*rsalt**2

				! jw
				! jw
				if(MSL-z_level(k) > 100.0) then
					! jw
					pres = 0.0 ! jw
					do kk = k, top_layer_at_face(j)
						ptmp = 9.81*rho_face(kk,j)*dz_face(kk,j) ! jw
						if(kk == k) then
							pres = pres + 1.0e-5 * 0.5*ptmp ! jw
						else
							pres = pres + 1.0e-5 * ptmp
						end if
					end do

					! jw
					! jw
					sbm = 19652.21 + 148.4206*rtemp - 2.327105*rtemp**2 + 1.360477e-2*rtemp**3 - 5.155288e-5*rtemp**4	&
					&  + pres * (3.239908 + 1.43713e-3*rtemp + 1.16092e-4*rtemp**2 - 5.77905e-7*rtemp**3) &
					&  + pres**2 * (8.50935e-5 - 6.12293e-6*rtemp + 5.2787e-8*rtemp**2)                  &
					&  + rsalt * (54.6746 - 0.603459*rtemp + 1.09987e-2*rtemp**2 - 6.1670e-5*rtemp**3)     &
					&  + dsqrt(rsalt)**3 * (7.944e-2 + 1.6483e-2*rtemp - 5.3009e-4*rtemp**2)              &
					&  + pres * rsalt * (2.2838e-3 - 1.0981e-5*rtemp - 1.6078e-6*rtemp**2)                &
					&  + 1.91075e-4*pres*dsqrt(rsalt)**3 - 9.9348e-7*pres**2*rsalt                       &
					&  + 2.08160e-8*rtemp*pres**2*rsalt + 9.1697e-10*rtemp**2*pres**2*rsalt
					
					ratio = 1.0 / (1.0 - pres/sbm)
					rho_face(k,j) = rho_face(k,j) * ratio ! jw
				end if ! jw

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
	!$omp do private(i,k,l,icount,nd)
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
end subroutine calculate_density_full  

