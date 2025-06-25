!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! 
!! This will find open boundary salinity and temperature condition for transport.
!! ===========================================================================!
subroutine find_openboundary_salt_temp_v3
	use mod_global_variables
	implicit none

	integer :: i, k, l, ii
	integer :: ie, j2, ii2
	! integer :: i3
	integer :: t_layer, b_layer

	! real(dp):: Q_jk_theta(4)
	real(dp):: Q_jk_theta
	real(dp),dimension(num_ob_cell) :: salt_at_obc, temp_at_obc
! 	real(dp),dimension(maxlayer,num_ob_cell) :: salt_at_obck, temp_at_obck
	real(dp),dimension(maxlayer) :: sum1, sum2
	integer :: icount
	
	! interpolation variables
	integer :: t2
	real(dp):: u1,u2,u3,v1,v2,v3
	! End of local variables ==================================================!

	
   ! ======================================================================!
   ! First, we have to prepare the salinity at open boundary elements,
   ! which is externally given, i.e., known value.
   ! ======================================================================!
   ! u2 is the current julian day, and it will be used for interpolation.
	! u2 = jday + (it*dt)/86400.0_dp ! time in [day] to be interpolated
	u2 = julian_day*86400.0 ! julian time [s], current time and/or time for looking
	
	salt_at_obc = ref_salt
	temp_at_obc = ref_temp
	
	! (1) at the tidal open boundary
	! write(*,*) is_salt, salt_ser_id(1), num_ob_cell
	! write(*,*) 'julian_day = ', julian_day
	! write(*,*) 'reference_diff_day = ', reference_diff_days
	do i=1,num_ob_cell
		! currently, we have four open boundary surface elevation types:
		! 		ob_eta_type(i) == -1, 1, 2
		! 		-1: radiation boundary condition
		!  	 1: harmonic tide (cosine wave)
		! 		 2: eta_ser.inp
		
		
		! Note, in this method, both salt & temp will be transported even though the option is turned off in C7_3
		! And, that is why I am doing this...
		! So, first, fill con_at_ob with initial value,
		! Then, if it is subjected to the transport, update the value with interpolated value at this time step.
		ie = ob_cell_id(i) ! open boundary element ID
		! salt_at_obc(i) = salt_cell(top_layer_at_element(ie),ie) ! let's just use surface salinity
		! temp_at_obc(i) = temp_cell(top_layer_at_element(ie),ie) ! let's just use surface temperature
		! write(*,*) i, salt_at_obc(i), temp_at_obc(i)
				
		! First, prepare the boundary concentration from the timeseries data
		! for salinity:
		if(is_salt == 1 .and. salt_ser_id(i) > 0) then 
			! if salt transport is turned on, prepare for the boundary condition.
			ii = salt_ser_id(i)
			do t2=2,salt_ser_data_num(ii)
				! Since salt_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
				! u1, and u3 should be shifted toward the simulation start year.
				! That is why I put "- reference_diff_days * 86400.0"
				u1 = salt_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! [s], lower bound time
				u3 = salt_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! [s], upper bound time
				
				! write(*,*) u1, u2, u3, reference_diff_days
				if(u2 >= u1 .and. u2 <= u3) then
					v1 = salt_ser_salt(t2-1,ii)		! [kg/m3], lower bound salt
					v3 = salt_ser_salt(t2  ,ii)		! [kg/m3], upper bound salt
					v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated salt
					exit
				end if			
			end do
		else
			v2 = ref_salt ! just use ref_salt as boundary salinity
		end if
		! salinity at ith open boundary element (at corresponding ghost cell)
		salt_at_obc(i) = spinup_function_baroclinic * v2
	
		! for temperature:
		if(is_temp == 1 .and. temp_ser_id(i) > 0) then
			! if temp transport is turned on, prepare for the boundary condition.
			ii = temp_ser_id(i)
			do t2=2,temp_ser_data_num(ii)
				! Since temp_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
				! u1, and u3 should be shifted toward the simulation start year.
				! That is why I put "- reference_diff_days * 86400.0"
				u1 = temp_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! [s], lower bound time
				u3 = temp_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! [s], upper bound time
				
				if(u2 >= u1 .and. u2 <= u3) then
					v1 = temp_ser_temp(t2-1,ii)		! [C], lower bound temp
					v3 = temp_ser_temp(t2  ,ii)		! [C], upper bound temp
					v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated temp
					exit
				end if			
			end do
		else
			v2 = ref_temp ! just use ref_temp as boundary temperature
		end if		
		! temperature at ith open boundary element (at corresponding ghost cell)
		temp_at_obc(i) = spinup_function_baroclinic * v2 ! let's just borrow spinup_function_baroclinic for temperature
			
		do k=1,maxlayer
			salt_at_obck(k,i) = salt_at_obc(i)
			temp_at_obck(k,i) = temp_at_obc(i)
		end do

		! Second, choose what to use based on flow direction
		! if inflow through boundary face, set ghost cell concentration as given boundary concentration.
		! if outflow through boundary face, set ghost cell concentration as the current boundary cell concentration (zero gradient approach).
		t_layer = top_layer_at_element(ie)
		b_layer = bottom_layer_at_element(ie)
		
		! calculate average values of neighbor element's values
		! I will use these values for outflow case
		! Note, this approach is very low level approach.
		! So, I will upgarde this approach later to a better one.
! 		icount = 0
! 		sum1 = 0.0_dp
! 		sum2 = 0.0_dp
! 		do l=1,tri_or_quad(ie)
! 			ii2 = adj_cellnum_at_cell(l,ie)
! 			if(ii2 > 0) then ! if adjacent cell exists.
! 				icount = icount + 1
! 				do k=b_layer,t_layer
! 					sum1(k) = sum1(k) + salt_cell(k,ii2)
! 					sum2(k) = sum2(k) + temp_cell(k,ii2)
! 				end do
! 			end if
! 		end do ! do l=1,tri_or_quad(ie)
! 		sum1 = sum1/max(icount,1) ! to avoid divide by 0
! 		sum2 = sum2/max(icount,1) ! to avoid divide by 0

! 		do l=1,tri_or_quad(ie)
! 			j2 = facenum_at_cell(l,ie)
! 			ii2 = adj_cellnum_at_cell(l,ie)
! 			if(boundary_type_of_face(j2) > 0) then ! if this is a boundary face, find if this face has outflow or inflow.
! 				do k=b_layer,t_layer
! 	         	Q_jk_theta = face_length(j2)*(dz_face(k,j2)) &
! 	         	&				*((1.0-theta)*un_face(k,j2)*sign_in_outflow(l,ie) + theta*un_face_new(k,j2)*sign_in_outflow(l,ie))
		         	
! 					if(Q_jk_theta < 0.0_dp) then
! 						! if inflow through this face, use given boundary concentration
! 						salt_at_obck(k,i) = salt_at_obc(i)
! 						temp_at_obck(k,i) = temp_at_obc(i)
! 					else
! 						! if outflow through this face, use neighbor cells' average for the boundary concentration.
! 						salt_at_obck(k,i) = sum1(k) ! salt_cell(k,ie)
! 						temp_at_obck(k,i) = sum2(k) ! temp_cell(k,ie)
! 					end if
! 				end do
				
				! if(i == 54) then
				! 	write(*,*) i,ie, Q_jk_theta, salt_at_obck(t_layer,i)
				! end if
! 			end if
! 		end do ! do l=1,tri_or_quad(ie)
		
		! extend salt & temp
		do k=1,b_layer-1
			salt_at_obck(k,i) = salt_at_obck(b_layer,i)
			temp_at_obck(k,i) = temp_at_obck(b_layer,i)
		end do
		do k=t_layer,maxlayer
			salt_at_obck(k,i) = salt_at_obck(t_layer,i)
			temp_at_obck(k,i) = temp_at_obck(t_layer,i)
		end do		
		
	end do ! do i=1,num_ob_cell
	! stop 'jww'
	
! 	do i=1,num_ob_cell
! 		! write(*,*) i, ob_cell_id(i), salt_at_obck(5,i), temp_at_obck(5,i)
! 		write(*,*) i, ob_cell_id(i), salt_at_obc(i), temp_at_obc(i)
! 	end do
! 	stop 'jw'
	
	! Now, I will distribute boundary cell info to the boundary nodes
! 	do i=1,num_ob_cell
! 		ie = ob_cell_id(i) ! open boundary element ID
! 		t_layer = top_layer_at_element(ie)
! 		b_layer = bottom_layer_at_element(ie)

! 		do k=b_layer,t_layer
! 			salt_at_obck2(k,1,i) = salt_at_obck(k,i)
! 			salt_at_obck2(k,2,i) = salt_at_obck(k,i)
! 			temp_at_obck2(k,1,i) = temp_at_obck(k,i)
! 			temp_at_obck2(k,2,i) = temp_at_obck(k,i)
! 		end do
		
		! extend temp & salt for below bottom & above top layers for boundary nodes
! 		do k=1,b_layer-1
! 			salt_at_obck2(k,1,i) = salt_at_obck(b_layer,i)
! 			salt_at_obck2(k,2,i) = salt_at_obck(b_layer,i)
! 			temp_at_obck2(k,1,i) = temp_at_obck(b_layer,i)
! 			temp_at_obck2(k,2,i) = temp_at_obck(b_layer,i)			
! 		end do
! 		do k=t_layer+1,maxlayer
! 			salt_at_obck2(k,1,i) = salt_at_obck(t_layer,i)
! 			salt_at_obck2(k,2,i) = salt_at_obck(t_layer,i)
! 			temp_at_obck2(k,1,i) = temp_at_obck(t_layer,i)
! 			temp_at_obck2(k,2,i) = temp_at_obck(t_layer,i)			
! 		end do		
! 	end do

	! (2) at the river boundary
	! Note:
	! 		At the river boundary, I do not specify vertical distribution.
	! 		And, I do not set boundary condition based on flow direction:
	! 			i.e., I just enforce the river boundary condition with given data.
	do i=1,num_Qb_cell
		! note: 
		! 		Q_boundary(i,1) = element id
		! 		Q_boundary(i,2) = face id
		! temp_element = Q_boundary(i,1)
		
		! for salinity
		if(is_salt == 1 .and. Q_salt_ser_id(i) > 0) then
			! if salt transport is turned on, prepare for the boundary condition.
			ii = Q_salt_ser_id(i)
			do t2=2,salt_ser_data_num(ii)
				! Since salt_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
				! u1, and u3 should be shifted toward the simulation start year.
				! That is why I put "- reference_diff_days * 86400.0"
				u1 = salt_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! [s], lower bound time
				u3 = salt_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! [s], upper bound time
				
				if(u2 >= u1 .and. u2 <= u3) then
					v1 = salt_ser_salt(t2-1,ii) 		! [kg/m3], lower bound salt
					v3 = salt_ser_salt(t2  ,ii) 		! [kg/m3], upper bound salt
					v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated salt
					exit
				end if
			end do
			! salinity at ith river boundary element (at corresponding ghost cell)
			salt_at_Qbc(i) = spinup_function_baroclinic * v2
		end if

		! for temperature
		if(is_temp == 1 .and. Q_temp_ser_id(i) > 0) then
			! if temp transport is turned on, prepare for the boundary condition.
			ii = Q_temp_ser_id(i)
			do t2=2,temp_ser_data_num(ii)
				! Since temp_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
				! u1, and u3 should be shifted toward the simulation start year.
				! That is why I put "- reference_diff_days * 86400.0"
				u1 = temp_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! [s], lower bound time
				u3 = temp_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! [s], upper bound time
				
				if(u2 >= u1 .and. u2 <= u3) then
					v1 = temp_ser_temp(t2-1,ii) 		! [C], lower bound salt
					v3 = temp_ser_temp(t2  ,ii) 		! [C], upper bound salt
					v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated temp
					exit
				end if
			end do
			! temperature at ith river boundary element (at corresponding ghost cell)
			temp_at_Qbc(i) = spinup_function_baroclinic * v2 ! let's just borrow spinup_function_baroclinic for temperature
		end if
	end do ! do i=1,num_Qb_cell
	! End of boundary concentration preparation ============================!

end subroutine find_openboundary_salt_temp_v3