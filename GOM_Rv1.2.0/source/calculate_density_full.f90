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
	!$omp do private(i,k,kk,rtemp,rsalt,pres,ptmp,sbm,ratio)
	do i = 1,maxnod
		! calculate density at wet nodes only, dry node not calculated
		if(top_layer_at_node(i) /= 0)then
			do k = bottom_layer_at_node(i), top_layer_at_node(i)
            rtemp = temp_node(k,i)   ! temperature at node
            rsalt = salt_node(k,i)   ! salinity at node
				
				! Note: here, I am giving a little more room for 0.00000 value (i.e., using -5.0 not exaxt 0.0)
				! 		  This trick is applied through end of this code...
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

				! density of sea water at at one standard atmosphere pressure (when Pair = 0)
				! here : 999.842594 = (1000.0 - 0.157406)
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
				
				! account for pressure for depth > 100m
				! if( .not. first_call .and. MSL-z_level(k) > 100.0) then
				if(MSL-z_level(k) > 100.0) then
					! water pressure in [bar], p term in the book
					pres = 0.0 ! total water pressure
					do kk = k, top_layer_at_node(i)
						ptmp = 9.81*rho_node(kk,i)*dz_node(kk,i) ! [Pa], rho * g * h
						if(kk == k) then
							! note: 1.e-5 is the unit conversion factor: [Pa] to [bar]
							pres = pres + 1.0e-5 * 0.5*ptmp ! [Pa] to [bar], at the bottom layer, we only have to include a half depth
						else
							pres = pres + 1.0e-5 * ptmp
						end if
					end do

					! For the IES 80, the secant bulk modulus is given by:
					! K(s,t,p) term in the book
					sbm = 19652.21 + 148.4206*rtemp - 2.327105*rtemp**2 + 1.360477e-2*rtemp**3 - 5.155288e-5*rtemp**4 &
					&  + pres * (3.239908 + 1.43713e-3*rtemp + 1.16092e-4*rtemp**2 - 5.77905e-7*rtemp**3) &
					&  + pres**2 * (8.50935e-5 - 6.12293e-6*rtemp + 5.2787e-8*rtemp**2)                  &
					&  + rsalt * (54.6746 - 0.603459*rtemp + 1.09987e-2*rtemp**2 - 6.1670e-5*rtemp**3)     &
					&  + dsqrt(rsalt)**3 * (7.944e-2 + 1.6483e-2*rtemp - 5.3009e-4*rtemp**2)              &
					&  + pres * rsalt * (2.2838e-3 - 1.0981e-5*rtemp - 1.6078e-6*rtemp**2)                &
					&  + 1.91075e-4*pres*dsqrt(rsalt)**3 - 9.9348e-7*pres**2*rsalt                       &
					&  + 2.08160e-8*rtemp*pres**2*rsalt + 9.1697e-10*rtemp**2*pres**2*rsalt
					
					ratio = 1.0 / (1.0 - pres/sbm)
					rho_node(k,i) = rho_node(k,i) * ratio ! this is the final density equation
				end if ! .not.first_call.and.MSL-z_level(k)>100

				if(rho_node(k,i) < 980.0) then
					write(pw_run_log,*) 'Water density is too low (weird density) at node:'
					write(pw_run_log,*) 'i, k, temp(i,k), salt(i,k), rho_node(k,i) = ', i, k, rtemp, rsalt, rho_node(k,i)
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
	!$omp do private(j,k,kk,rtemp,rsalt,pres,ptmp,sbm,ratio)
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

				! density of sea water at at one standard atmosphere pressure (when Pair = 0)
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

				! account for pressure for depth > 100m
				! if( .not. first_call .and. MSL-z_level(k) > 100.0) then
				if(MSL-z_level(k) > 100.0) then
					! water pressure, p term in the book
					pres = 0.0 ! total water pressure
					do kk = k, top_layer_at_face(j)
						ptmp = 9.81*rho_face(kk,j)*dz_face(kk,j) ! rho * g * h
						if(kk == k) then
							pres = pres + 1.0e-5 * 0.5*ptmp ! at the bottom layer, we only have to include a half depth
						else
							pres = pres + 1.0e-5 * ptmp
						end if
					end do

					! For the IES 80, the secant bulk modulus is given by:
					! K(s,t,p) term in the book
					sbm = 19652.21 + 148.4206*rtemp - 2.327105*rtemp**2 + 1.360477e-2*rtemp**3 - 5.155288e-5*rtemp**4	&
					&  + pres * (3.239908 + 1.43713e-3*rtemp + 1.16092e-4*rtemp**2 - 5.77905e-7*rtemp**3) &
					&  + pres**2 * (8.50935e-5 - 6.12293e-6*rtemp + 5.2787e-8*rtemp**2)                  &
					&  + rsalt * (54.6746 - 0.603459*rtemp + 1.09987e-2*rtemp**2 - 6.1670e-5*rtemp**3)     &
					&  + dsqrt(rsalt)**3 * (7.944e-2 + 1.6483e-2*rtemp - 5.3009e-4*rtemp**2)              &
					&  + pres * rsalt * (2.2838e-3 - 1.0981e-5*rtemp - 1.6078e-6*rtemp**2)                &
					&  + 1.91075e-4*pres*dsqrt(rsalt)**3 - 9.9348e-7*pres**2*rsalt                       &
					&  + 2.08160e-8*rtemp*pres**2*rsalt + 9.1697e-10*rtemp**2*pres**2*rsalt
					
					ratio = 1.0 / (1.0 - pres/sbm)
					rho_face(k,j) = rho_face(k,j) * ratio ! this is the final density equation
				end if ! .not.first_call.and.MSL-z_level(k)>100

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
	!$omp do private(i,k,l,icount,nd)
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
end subroutine calculate_density_full  

