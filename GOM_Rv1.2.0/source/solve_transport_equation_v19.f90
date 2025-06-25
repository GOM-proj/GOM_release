!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! solve transport equation with TVD
!! 
subroutine solve_transport_equation_v19
	use mod_global_variables
	! use mod_file_definition 
	implicit none

	integer :: i, j, k, l, ii, n1, n2
	integer :: i2, j2, ii2, i3, ie
	integer :: kk, num_vertical_layer
	integer :: t_layer, b_layer
	integer :: count1

	real(dp),dimension(0:maxlayer+1,maxele,maxtran2) :: con_cell ! including top/bottom ghost layers
	real(dp),dimension(maxlayer,maxele,maxtran2) :: con_cell_new
	real(dp),dimension(maxlayer,maxnod,maxtran2) :: con_node_new
	real(dp),dimension(maxlayer,maxface,maxtran2) :: con_face_new

	! store time series data from input:
	integer :: ith_ob, ith_Qb, ith_WR, ith_SS ! ith boundary
	real(dp),dimension(num_ob_cell,maxtran2) :: con_at_ob
	real(dp),dimension(num_Qb_cell,maxtran2) :: con_at_Qb
	real(dp),dimension(num_WR_cell,maxtran2) :: con_at_WR
	real(dp),dimension(num_SS_cell,maxtran2) :: con_at_SS
	real(dp),dimension(maxlayer,num_ob_cell,maxtran2) :: con_at_obk ! vertically distributed open boundary concentration
	real(dp),dimension(maxlayer,num_Qb_cell,maxtran2) :: con_at_Qbk ! vertically distributed river boundary concentration
	real(dp),dimension(maxlayer,num_WR_cell,maxtran2) :: con_at_WRk
	real(dp),dimension(maxlayer,num_SS_cell,maxtran2) :: con_at_SSk
	
	real(dp),dimension(maxlayer,maxtran2):: con_cell_iik_oQb ! boundary concentration at kth layer (either open_boundary/river_boundary)
	real(dp),dimension(maxtran2):: con_cell_iik ! final neighbor cell concentration at kth layer (either water neighbor/open_boundary/river_boundary)

	! advection/diffusion coefficients:	
	real(dp),dimension(maxlayer,maxface):: Q_jk_theta, D_jk_old, d_jk_theta ! here, Q_jk_theta is based on face normal coordinate.
	real(dp),dimension(0:maxlayer,maxele):: Q_ik_theta, D_ik_old, d_ik_theta
	real(dp) :: Q_jk_theta2(4) ! this includes inflow/outflow status

	! consecutive gradients:
	real(dp),dimension(maxlayer,maxface,maxtran2) :: r_jk
	real(dp),dimension(0:maxlayer,maxele,maxtran2):: r_ik
	real(dp),dimension(maxtran2):: sum1, sum1_all, sum2, sum2_all
	
	! flux limiter functions:
	real(dp),dimension(maxlayer,maxface) :: small_phi_jk
	real(dp),dimension(0:maxlayer,maxele):: small_phi_ik
	real(dp),dimension(maxlayer,maxface,maxtran2):: large_PHI_jk 
	real(dp),dimension(0:maxlayer,maxele,maxtran2):: large_PHI_ik
	real(dp),dimension(maxlayer,maxface,maxtran2):: PSI_jk
	real(dp),dimension(0:maxlayer,maxele,maxtran2):: PSI_ik
	
	! advection/diffusion terms:
	real(dp),dimension(maxtran2):: h_flux_term, v_flux_term
	real(dp),dimension(maxtran2):: h_adv, v_adv, h_diff
	real(dp),dimension(maxtran2):: sum_QC_i, sum_QC_ii, sum_hdiff, sum_hlimiter

	! tridiagonal matrix:
	real(dp),dimension(maxlayer,maxtran2):: rhs, solution
	real(dp),dimension(maxlayer):: a_lower_mat, b_diagonal_mat, c_upper_mat, gam

	! interpolation variables:
	integer :: t2
	real(dp):: u1,u2,u3,v1,v2,v3
	
	! subcycling index:
	integer :: subcycle, Nt
	real(dp):: dt2
	
	! etc:
	real(dp):: cp = 4186.0 ! specific_heat of water [J/KgC]
	real(dp):: dz_old, dz_new ! linearly interpolated dz in subcycling
	real(dp):: rtemp, rtemp1, rtemp2
	real(dp):: sum11, sum22
	! End of local variables ==================================================!
		
	! initialize before start:
	con_cell_new = 0.0_dp
	con_node_new = 0.0_dp
	con_face_new = 0.0_dp
	
	! copy active transport materials to con_cell
	! Note:
	!     "con_cell" includes top/bottom ghost layers. 
	!     I include these ghost layers to avoid array dimension mismatch when calculating vertical advection/diffusion.
	do i3=1,maxtran2
		if(tran_id(i3) == 1) then
			do i=1,maxele
				do k=1,maxlayer
					con_cell(k,i,i3) = salt_cell(k,i)
				end do
	 			! extend to top/bottom ghost cells:
				con_cell(0,i,i3) = 0.0_dp
				con_cell(maxlayer+1,i,i3) = 0.0_dp
			end do
		end if
		if(tran_id(i3) == 2) then
			if(heat_option == 1) then
				! use term-by-term calculation
				do i=1,maxele
					do k=1,maxlayer
						con_cell(k,i,i3) = temp_cell(k,i)
					end do
		 			! extend to top/bottom ghost cells:
					con_cell(0,i,i3) = 0.0_dp
					con_cell(maxlayer+1,i,i3) = 0.0_dp
					
					! include source terms here --------------------------------------!
					t_layer = top_layer_at_element(i)
					b_layer = bottom_layer_at_element(i)
					
					
! 	 				if(i==1) then
! 	 					rtemp = dt*phi_n(i)/(rho_o*cp*dz_cell(t_layer,i))
! 	 					! write(*,*) 'it_1, con_cell(t_layer,i,i3) = ', it, con_cell(t_layer,i,i3), dt, phi_n(i), dz_cell(t_layer,i), rtemp ! phi_n(i) is the problem now...
! 	 					write(*,*) 'it_1, con_cell(t_layer,i,i3) = ', it, con_cell(t_layer,i,i3), phi_n(i) ! phi_n(i) is the problem now...
! 	 				end if
					
! jw					con_cell(t_layer,i,i3) = con_cell(t_layer,i,i3) + dt*phi_n(i)/(rho_o*cp*dz_cell(t_layer,i))
					
					! add sediment back radiation to all water columns
					! con_cell(b_layer,i,i3) = con_cell(b_layer,i,i3) + dt*phi_sw(i)/(rho_o*cp*dz_cell(t_layer,i))
! jw					do k=1,maxlayer
! jw					 	con_cell(k,i,i3) = con_cell(k,i,i3) + dt*phi_sw(i)/(rho_o*cp*dz_cell(t_layer,i))
! jw					end do
					
					
					do k=b_layer,t_layer
						con_cell(k,i,i3) = con_cell(k,i,i3) + dt*phi_sz(k,i)/(rho_o*cp*dz_cell(k,i))
					end do
					
					
					! if(i==1) then
					! 	write(*,*) 'it_2, con_cell(t_layer,i,i3) = ', it, con_cell(t_layer,i,i3)
					! end if				
				end do
			else if(heat_option == 2) then ! not yet included... Note... this part should be deleted later...
				! use equilibrium temperature
				do i=1,maxele
					do k=1,maxlayer
						con_cell(k,i,i3) = temp_cell(k,i)
					end do
		 			! extend to top/bottom ghost cells:
					con_cell(0,i,i3) = 0.0_dp
					con_cell(maxlayer+1,i,i3) = 0.0_dp
					
					! include source terms here --------------------------------------!
					t_layer = top_layer_at_element(i)
					b_layer = bottom_layer_at_element(i)
					
					! here I have to include equilibrium temperature:
					con_cell(t_layer,i,i3) = con_cell(t_layer,i,i3) ! <here> I have to include equilibrium temperature ! dt*phi_n(i)/(rho_o*cp*dz_cell(t_layer,i))

					! add sediment back radiation to all water columns
					do k=1,maxlayer
						con_cell(k,i,i3) = con_cell(k,i,i3) + dt*phi_sw(i)/(rho_o*cp*dz_cell(t_layer,i))
					end do
	
					! if(i==1) then
					! 	write(*,*) 'it_2, con_cell(t_layer,i,i3) = ', it, con_cell(t_layer,i,i3)
					! end if				
				end do				
			end if
			
			! check...
			! rtemp1 = dt*phi_n(1)/(rho_o*cp*dz_cell(t_layer,1))
			! write(*,*) 'jw1'
			! t_layer = top_layer_at_element(1)
			! write(*,*) it, dt, phi_n(1),rho_o,cp,t_layer,dz_cell(t_layer,1), dz_cell_new(t_layer,1), rtemp1
		end if
	end do
	

	Nt = trans_sub_iter ! number of sub-cycling
	
	! if trans_sub_iter is given with "0", reset this number to "1" to do "1" iteration.
	! Now, following part will be executed at least "one" time.
	! so, from now, I don't need to keep the original non-subcycling version.
	if(Nt == 0) then
		Nt = 1
	end if

	! jw = 0.0_dp
	
	
	dt2 = dt/Nt ! new dt for sub-cycling
	do subcycle=0,(Nt-1) ! subcycle starts
		! Note:
		! 		I simply use Q for this time step at the outer time step not this sub timestep.
		! 		i.e., I update Q in "compute_boundary_velocity.f90" not in this sub timestep.
		! 		But, I update concentration at this each sub timestep.
		! 		This is for the speedup.
		
		
		! Before start, let's initialize some variables ========================!
		! initialize open boundary and river boundary concentrations:
		con_at_ob = 0.0_dp
		con_at_Qb = 0.0_dp
		con_at_WR = 0.0_dp
		con_at_SS = 0.0_dp
		con_at_obk = 0.0_dp
		con_at_Qbk = 0.0_dp
		con_at_WRk = 0.0_dp
		con_at_SSk = 0.0_dp

		! initialize r_jk & r_ik:
  		r_ik = 1.0_dp
 		r_jk = 1.0_dp

		! setup coefficient matrix, which will have a tridiagonal system on each horizontal grid element
		! no boundary conditions are involved yet.
	   a_lower_mat = 0.0_dp
	   b_diagonal_mat = 0.0_dp
	   c_upper_mat = 0.0_dp
	   gam = 0.0_dp   
		rhs = 0.0_dp
		solution = 0.0_dp

	   ! ======================================================================!
	   ! First, we have to prepare the transport materials (e.g., salt, temp, ...) at open boundary elements,
	   ! which is externally given, i.e., known value.
	   ! ======================================================================!
	   ! u2 is the current julian day, and it will be used for interpolation.
		! u2 = jday + (it*dt)/86400.0_dp ! time in [day] to be interpolated
		! Note: julian_day is the current time for outer simulation.
		! thus, we have to go back to previous time, and add subcycling time from there.
		u2 = (julian_day*86400.0 - dt) + (subcycle+1)*dt2 ! julian time [s], current time and/or time for looking
		
		! (1) at the tidal open boundary
		do i=1,num_ob_cell
			! currently, we have four open boundary surface elevation types:
			! 		ob_eta_type(i) == -1, 1, 2
			! 		-1: radiation boundary condition
			!  	 1: harmonic tide (cosine wave)
			! 		 2: eta_ser.inp
			
			! First, prepare the boundary concentration from the timeseries data,
			! and interpolate it to the current time frame.
			do i3=1,maxtran2
				! for salinity:
				if(tran_id(i3) == 1 .and. salt_ser_id(i) > 0) then
					ii = salt_ser_id(i)
					do t2=2,salt_ser_data_num(ii)
						! Since salt_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
						! u1, and u3 should be shifted toward the simulation start year.
						! That is why I put "- reference_diff_days * 86400.0"
						u1 = salt_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! [s], lower bound time
						u3 = salt_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! [s], upper bound time
						
						if(u2 >= u1 .and. u2 <= u3) then
							v1 = salt_ser_salt(t2-1,ii)		! [kg/m3], lower bound salt
							v3 = salt_ser_salt(t2  ,ii)		! [kg/m3], upper bound salt
							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated salt
							exit
						end if			
					end do
					! salinity at ith open boundary element (at corresponding ghost cell)
					con_at_ob(i,i3) = spinup_function_baroclinic * v2
				end if
			
				! for temperature:
				if(tran_id(i3) == 2 .and. temp_ser_id(i) > 0) then
					ii = temp_ser_id(i)
					do t2=2,temp_ser_data_num(ii)
						! Since temp_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
						! u1, and u3 should be shifted toward the simulation start year.
						! That is why I put "- reference_diff_days * 86400.0"
						u1 = temp_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! [s], lower bound time
						u3 = temp_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! [s], upper bound time
						
						if(u2 >= u1 .and. u2 <= u3) then
							v1 = temp_ser_temp(t2-1,ii)		! [C], lower bound salt
							v3 = temp_ser_temp(t2  ,ii)		! [C], upper bound salt
							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated temp
							exit
						end if			
					end do
					! temperature at ith open boundary element (at corresponding ghost cell)
					con_at_ob(i,i3) = spinup_function_baroclinic * v2 ! let's just borrow spinup_function_baroclinic for temperature
				end if
			end do ! do i3=1,maxtran2, end of time interpolation

			! Second, choose what to use based on flow direction, and distribute them in vertical columns. 
			! if inflow through boundary face, set ghost cell concentration as given boundary concentration.
			! if outflow through boundary face, set ghost cell concentration as the current boundary cell concentration (zero gradient approach).
			! Then, distribute them in vertical, i.e.:
			! 		"con_cell" 	-> con_at_obk	if outflow or zero flow from boundary cell
			! 		"con_at_ob" -> con_at_obk	if inflow to boundary cell
			i2 = ob_cell_id(i) ! element ID (i.e., element number)
			! j2 = ob_face_id(i) ! face ID, I will find this in the next "do loop"
			t_layer = top_layer_at_element(i2)
			b_layer = bottom_layer_at_element(i2)

			do l=1,tri_or_quad(i2)
				j2 = facenum_at_cell(l,i2)
				if(boundary_type_of_face(j2) > 0) then ! if this is a boundary face, find if this face has outflow or inflow.
					do k=b_layer,t_layer
						! Note:
						!      Here, "dz_face" has been linearly interpolated between n and n+1 time step
						!      i.e., linearly interpolated at current sub-time-step.
						!      i.e., "dz_old" at current sum_time_step.
						dz_old = (dz_face(k,j2) + (dz_face_new(k,j2) - dz_face(k,j2))*(subcycle)/Nt)
						Q_jk_theta2(l) = face_length(j2)* dz_old &
						&				*((1.0-theta)*un_face(k,j2)*sign_in_outflow(l,i2) + theta*un_face_new(k,j2)*sign_in_outflow(l,i2))
						
						if(Q_jk_theta2(l) <= 0.0_dp) then ! if inflow or no flow through this face, use given boundary concentration
							do i3=1,maxtran2
								con_at_obk(k,i,i3) = con_at_ob(i,i3) ! distribute in vertical
							end do
						else ! if outflow through this face, use self concentration for the boundary concentration (i.e., zero gradient condtion).
							do i3=1,maxtran2
								con_at_obk(k,i,i3) = con_cell(k,i2,i3) ! distribute in vertical
								! con_at_obk(k,i,i3) = con_at_ob(i,i3) ! distribute in vertical
							end do
						end if
					end do
				end if
			end do ! do l=1,tri_or_quad(i2)
		end do ! do i=1,num_ob_cell

		! (2) at the river boundary
		! Note:
		! 		At the river boundary, I do not set boundary condition based on flow direction:
		! 			i.e., I just enforce the river boundary condition with given data.
		do i=1,num_Qb_cell
			! note: 
			! 		Q_boundary(i,1) = element id
			! 		Q_boundary(i,2) = face id
			! temp_element = Q_boundary(i,1)
			
			! First, prepare the boundary concentration from the timeseries data
			do i3=1,maxtran2
				! for salinity
				if(tran_id(i3) == 1 .and. Q_salt_ser_id(i) > 0) then
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
					con_at_Qb(i,i3) = spinup_function_baroclinic * v2
					
					! write(*,*) 'con_at_Qb(i,i3) = ', i, i3, con_at_Qb(i,i3)
				end if
				
				! for temperature
				if(tran_id(i3) == 2 .and. Q_temp_ser_id(i) > 0) then
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
					con_at_Qb(i,i3) = spinup_function_baroclinic * v2 ! let's just borrow spinup_function_baroclinic for temperature
					! write(*,*) 'con_at_Qb = ', i, i3, con_at_Qb(i,i3)
					! stop 'jw'
				end if
			end do ! do i3=1,maxtran2
			
			! Second, distribute them in vertical. 
			! Note: I do not set boundary condition based on flow direction for river boundary:
			! 		i.e., I just enforce the river boundary condition with given data.
			i2 = Q_boundary(i,1) ! element id
			
			t_layer = top_layer_at_element(i2)
			b_layer = bottom_layer_at_element(i2)
			
			do k=b_layer,t_layer
				do i3=1,maxtran2
					con_at_Qbk(k,i,i3) = con_at_Qb(i,i3)
					
					! write(*,*) 'con_at_Qbk(k,i,i3) = ', k,i,i3, con_at_Qb(i,i3)
				end do
			end do
		end do ! do i=1,num_Qb_cell

		! (3) at the Withdraw/Return boundary
		! Note:
		! 		It is similar to the river boundary, but I have to set boundary condition based on flow direction:
		do i=1,num_WR_cell
			! note: 
			! 		WR_boundary(i,1) = element id
			! 		WR_boundary(i,2) = face id
			! temp_element = WR_boundary(i,1)
			
			! First, prepare the boundary concentration from the timeseries data
			do i3=1,maxtran2
				! for salinity
				if(tran_id(i3) == 1 .and. WR_salt_ser_id(i) > 0) then
					ii = WR_salt_ser_id(i)
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
					! salinity at ith withdraw/return boundary element (at corresponding ghost cell)
					con_at_WR(i,i3) = spinup_function_baroclinic * v2
					! write(*,*) 'con_at_WR(i,i3) = ', i, i3, con_at_WR(i,i3)
				end if
				
				! for temperature
				if(tran_id(i3) == 2 .and. WR_temp_ser_id(i) > 0) then
					ii = WR_temp_ser_id(i)
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
					! temperature at ith withdraw/return boundary element (at corresponding ghost cell)
					con_at_WR(i,i3) = spinup_function_baroclinic * v2 ! let's just borrow spinup_function_baroclinic for temperature
					! write(*,*) 'con_at_WR(i,i3) = ', i, i3, con_at_WR(i,i3)
					! stop 'jw'
				end if
			end do ! do i3=1,maxtran2

			! Second, choose what to use based on flow direction.
			! if inflow through boundary face, set ghost cell concentration as given boundary concentration.
			! if outflow through boundary face, set ghost cell concentration as the current boundary cell concentration (zero gradient approach).
			i2 = WR_boundary(i,1)! element ID (i.e., element number)
			j = WR_boundary(i,2) ! face number of the WR boundary

			if(WR_layer(i) == 999) then
				! anytime surface layer
				k = top_layer_at_face(j)
			else if(WR_layer(i) == 0) then
				! anytime bottom layer
				k = bottom_layer_at_face(j)
			else
				k = WR_layer(i)
			end if
			
			if(WRu_boundary(i) >= 0.0) then ! outflow, i.e., withdrawl flow
				! if outflow through this face, use self concentration for the boundary concentration (i.e., zero gradient condtion).
				do i3=1,maxtran2
					con_at_WRk(k,i,i3) = con_cell(k,i2,i3)
					! write(*,*) 'con_at_WRk(k,i,i3) = ', k, i, i3, con_at_WRk(k,i,i3)
				end do
			else ! inflow, i.e., return flow
				! if inflow or no flow through this face, use given boundary concentration
				do i3=1,maxtran2
					con_at_WRk(k,i,i3) = con_at_WR(i,i3)
					! write(*,*) 'con_at_WRk(k,i,i3) = ', k, i, i3, con_at_WRk(k,i,i3)
				end do
			end if
		end do ! do i=1,num_WR_cell
		

		! (4) at the Source/Sink boundary
 		do i=1,num_SS_cell			
 			! First, prepare the boundary concentration from the timeseries data
 			do i3=1,maxtran2
 				! for salinity
 				if(tran_id(i3) == 1 .and. SS_salt_ser_id(i) > 0) then
 					ii = SS_salt_ser_id(i)
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
 					! salinity at ith source/sink boundary element (at corresponding ghost cell)
 					con_at_SS(i,i3) = spinup_function_baroclinic * v2
 				end if
 				
 				! for temperature
 				if(tran_id(i3) == 2 .and. SS_temp_ser_id(i) > 0) then
 					ii = SS_temp_ser_id(i)
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
 					! temperature at ith source/sink boundary element (at corresponding ghost cell)
 					con_at_SS(i,i3) = spinup_function_baroclinic * v2
 				end if
 			end do ! do i3=1,maxtran2
 		end do ! do i=1,num_SS_cell
		! End of boundary concentration preparation ============================!
		
		
		
		! Let's pre-calculate flux limiter function related variables ==========!
		! 		(1) Q_jk_theta, D_jk_old, d_jk_theta, small_phi_jk
		! 		(2) Q_ik_theta, D_ik_old, d_ik_theta, small_phi_ik
		! 		(3) r_jk & r_ik
		! 		(4) large_PHI_ik & PSI_ik
		! 		(5) large_PHI_jk & PSI_jk
		! Note: upwind scheme does not need:
		! 		small_phi_jk, small_phi_ik, r_jk, r_ik, large_PHI_ik & PSI_ik, large PHI_jk & PSI_jk
		! 		However, I will calculate some parts even for the upwind method since it is faster than using "if" statement.
		! omp starts here...
 		!$omp parallel
		! (1) calculate: Q_jk_theta, D_jk_old, d_jk_theta, small_phi_jk:
		!$omp do private(j, k, dz_old, b_layer, t_layer)
			do j=1,maxface
				b_layer = bottom_layer_at_face(j)
				t_layer = top_layer_at_face(j)
				
				do k=b_layer,t_layer
				! do k=1,maxlayer
					! Note:
					!     Here, I do not include "sign_in_outflow".
					!     So, this flow at the face is regarding to the face normal direction only.
					!     i.e., it does not tell the inflow/outflow through this face to this element,
					!     but it only tells that the amount of flow through this face.
					!     The inflow/outflow will be calculated whenever it requires later.
					!     The following comment equation is the original equation:
					! Q_jk_theta(k,j) = face_length(j)*(dz_face(k,j) + (dz_face_new(k,j) - dz_face(k,j))*(subcycle)/Nt) &
					! &						*((1.0-theta)*un_face(k,j)*sign_in_outflow(l,i) + theta*un_face_new(k,j)*sign_in_outflow(l,i))
					dz_old = dz_face(k,j) + (dz_face_new(k,j) - dz_face(k,j))*(subcycle)/Nt
					Q_jk_theta(k,j) = face_length(j)*dz_old &
					&						*((1.0-theta)*un_face(k,j) + theta*un_face_new(k,j))
					D_jk_old(k,j) = face_length(j) * dz_old * Kh(k,j) &
					&					/delta_j(j)
					d_jk_theta(k,j) = max(0.0, D_jk_old(k,j) - 0.5*abs(Q_jk_theta(k,j)))
					
					! to avoid "divide by 0"
					if(Q_jk_theta(k,j) == 0.0_dp) then
						small_phi_jk(k,j) = 1.0_dp
					else
						small_phi_jk(k,j) = min(1.0, 2.0*D_jk_old(k,j)/abs(Q_jk_theta(k,j)))
					end if
				end do ! do k=1,maxlayer
			end do ! do j=1,maxface
		!$omp end do
		
		! (2) calculate: Q_ik_theta, D_ik_old, d_ik_theta, small_phi_ik:
		!$omp do private(i, k, dz_old, b_layer, t_layer)
			do i=1,maxele
				b_layer = bottom_layer_at_element(i)
				t_layer = top_layer_at_element(i)
				! do k=0,maxlayer
				! do k=b_layer-1,t_layer
				do k=b_layer,t_layer-1 ! mid-levels only (i.e., excluding at m-1 & M levels)
					! Note:
					!      The value of the "wn_cell" at below bottom layers are 0.0 (see solve_velocities.f90),
					!      thus Q_ik_theta at below bottom layers also will be 0.0.
					Q_ik_theta(k,i) = area(i)*((1.0-theta)*wn_cell(  k,i) + theta*wn_cell_new(  k,i))
					
					! Note: D_ik_old(0,i) will be 0.0 since Kv(0,i) = 0.0 (see allocate_variables.f90 & read_main_inp.f90)
					dz_old = dzhalf_cell(k,i) + (dzhalf_cell_new(k,i) - dzhalf_cell(k,i))*(subcycle)/Nt
					D_ik_old(k,i) = area(i)*Kv(k,i)/dz_old
					d_ik_theta(k,i) = max(0.0, D_ik_old(k,i) - 0.5*abs(Q_ik_theta(k,i)))
						
					! to avoid "divide by 0"
					if(Q_ik_theta(k,i) == 0.0_dp) then
						small_phi_ik(k,i) = 1.0_dp
					else
						small_phi_ik(k,i) = min(1.0, 2.0*D_ik_old(k,i)/abs(Q_ik_theta(k,i)))
					end if
				end do ! do k=0,maxlayer				
			end do ! do i=1,maxele
		!$omp end do
		!$omp end parallel
		
		! (3) calculate: r_jk & r_ik:
		! 		calculate gradient terms only of flux_limiter option is used...
		if(h_flux_limiter > 0 .or. v_flux_limiter > 0) then
			!$omp parallel
	 		!$omp do private(i, j, k, l, ii, i3,											&
	 		!$omp & 			  b_layer, t_layer, num_vertical_layer,					&
	 		!$omp &			  ith_ob, ith_Qb,	ith_WR,										&
	 		!$omp &			  con_cell_iik, con_cell_iik_oQb,							&
	 		!$omp & 			  sum1, sum1_all, sum2, sum2_all,							&
	 		!$omp &			  Q_jk_theta2,														&
	 		!$omp &			  rtemp, rtemp1, rtemp2)
			do i=1,maxele
				t_layer = top_layer_at_element(i)
				b_layer = bottom_layer_at_element(i)
				
				if(t_layer == 0) then
					! if this is a dry element, do nothing, and move to the next element.
					cycle
				end if
				num_vertical_layer = t_layer - b_layer + 1 ! number of active vertical layers at this element
				
				! First, if this element belongs to tidal open boundary or river boundary,
				! copy boundary concentration value to the concentration at the neighbor element 
				! Note, that I am setting vertically identical open boundary ghost value for river boundary, i.e., no vertical difference.
				! But, I am using vertically different open boundary ghost value for tidal boundary (just for outflow case).
				ith_ob = ob_element_flag(i) ! This is the ith open boundary
				ith_Qb = Qb_element_flag(i) ! This is the ith river boundary
				ith_WR = WR_element_flag(i) ! This is the ith withdrawl/return boundary
				
				if(ith_ob > 0) then ! if this element belongs to tidal open boundary
		 			do k=b_layer,t_layer
						do i3=1,maxtran2
							con_cell_iik_oQb(k,i3) = con_at_obk(k,ith_ob,i3) ! salt, temp, ...
						end do
					end do
				end if
				if(ith_Qb > 0) then ! if this element belongs to river boundary
					do k=b_layer,t_layer
						do i3=1,maxtran2
							con_cell_iik_oQb(k,i3) = con_at_Qbk(k,ith_Qb,i3) ! salt, temp, ...
							! write(*,*) con_cell_iik_oQb(k,i3)
						end do
					end do
				end if
				if(ith_WR > 0) then ! if this element belongs to withdraw/return boundary
					! i2 = WR_boundary(i,1)! element ID (i.e., element number)
					j = WR_boundary(i,2) ! face number of the WR boundary
					if(WR_layer(ith_WR) == 999) then
						! anytime surface layer
						k = top_layer_at_face(j)
					else if(WR_layer(ith_WR) == 0) then
						! anytime bottom layer
						k = bottom_layer_at_face(j)
					else
						! specific vertical layer
						k = WR_layer(ith_WR)
					end if
					
					do i3=1,maxtran2
						con_cell_iik_oQb(k,i3) = con_at_WRk(k,ith_WR,i3)
						! write(*,*) 'con_cell_iik_oQb(k,i3) = ', k,i3,con_cell_iik_oQb(k,i3)
					end do
				end if

				! calculate r_jk ----------------------------------------------------!
				do k=b_layer,t_layer
					! calculate numerator and denominator terms in the "r" factor.
					! this should be done for inflow faces in this element.
					sum1 		= 0.0_dp ! sum for horizontal flux limiter term for leaving faces, in r_jk term, numerator term
					sum2 		= 0.0_dp ! sum for horizontal flux limiter term for leaving faces, in r_jk term, denominator term
					sum1_all = 0.0_dp ! sum1 + vertical term
					sum2_all = 0.0_dp ! sum2 + vertical term
					
					! first, include vertical terms:
					! To successfully kill unwanted terms in the original full equation, I use the following approach.
					! There were other approaches I used, but I think this is dirty but the most straightforward approach.
					if(num_vertical_layer == 1) then
						! if this is a 2D element, there will be no vertical terms.
					else
						! if this is a 3D element, there will be vertical terms.
						if(k == b_layer) then
							do i3=1,maxtran2
								rtemp2 = abs(0.5*(Q_ik_theta(k  ,i) + abs(Q_ik_theta(k  ,i))))
								sum1_all(i3) = rtemp2 * (con_cell(k,i,i3) - con_cell(k+1,i,i3))
								sum2_all(i3) = rtemp2
							end do
						else if(k == t_layer) then
							do i3=1,maxtran2
								rtemp1 = abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))
								sum1_all(i3) = rtemp1 * (con_cell(k,i,i3) - con_cell(k-1,i,i3))
								sum2_all(i3) = rtemp1
							end do						
						else
							do i3=1,maxtran2
								rtemp1 = abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))
								rtemp2 = abs(0.5*(Q_ik_theta(k  ,i) + abs(Q_ik_theta(k  ,i))))
								sum1_all(i3) = rtemp1 * (con_cell(k,i,i3) - con_cell(k-1,i,i3)) &
								&				 + rtemp2 * (con_cell(k,i,i3) - con_cell(k+1,i,i3))
								sum2_all(i3) = rtemp1 + rtemp2
							end do
						end if
					end if
					
					! second, include horizontal terms:
					! do this for inflow faces in this element.
					do l=1,tri_or_quad(i)
						j = facenum_at_cell(l,i)
						ii = adj_cellnum_at_cell(l,i)
						if(boundary_type_of_face(j) /= -1 .or. isflowside3(j) > 0 .or. isflowside4(j) > 0) then ! this face is not a land boundary face (i.e., neighbour is water cell or open boundary) or water boundary face
							! Note: 
							! 		boundary_type_of_face(j):
							! 			land = -1; 
							! 			water =  0; 
							! 			open boundary = 1 ~ n
							! so, if this face is not a land boundary face, we have flow and can calculate gradient
							! so, following will be calculated only if this face has neighbor water cell or this is one of the boundary faces
							! Note, river boundary face is also set to "land" boundary in boundary_type_of_face,
							! and that is why I need additional criteria, isflowside3 to include river boundary face
							
							! Note: 
							!      Previously, Q_jk_theta(k,j) was calculated without "sign_in_outflow",
							!      thus, it represented the flow based on the face normal direction.
							!      Now, I change it to inflow/outflow to/from this element
							! Q_jk_theta(k,j) = face_length(j)*(dz_face(k,j) + (dz_face_new(k,j) - dz_face(k,j))*(subcycle)/Nt) &
							! &						*((1.0-theta)*un_face(k,j) + theta*un_face_new(k,j))
							Q_jk_theta2(l) = Q_jk_theta(k,j)*sign_in_outflow(l,i)
								
							if(Q_jk_theta2(l) < 0.0_dp) then ! inflow through this face
								! find neighbor cell's concentration either from a true neighbor cell or a ghost boundary cell
								! so, here "con_cell_iik(i3)" is this element's (in vertical also) neighbor cell's concentration 
								! either the neighbor cell is boundary ghost cell or just normal water cell.
								if(ii > 0) then	! if adjacent element exists (water cell exists)
									do i3=1,maxtran2
										con_cell_iik(i3) = con_cell(k,ii,i3)
									end do
								else ! if adjacent element doesn't exist (i.e., it is a boundary face)
									do i3=1,maxtran2
										con_cell_iik(i3) = con_cell_iik_oQb(k,i3) ! i.e., either con_at_obk(k,ibnd) or con_at_Qbk(k,ibnd)
									end do
								end if
								
								do i3=1,maxtran2 ! for both salt & temp
									sum1(i3) = sum1(i3) + abs(Q_jk_theta2(l))*(con_cell(k,i,i3) - con_cell_iik(i3))
									sum2(i3) = sum2(i3) + abs(Q_jk_theta2(l))
								end do
							end if ! if(Q_jk_theta2 < 0.0_dp) then ! inflow through this face
						end if ! if(boundary_type_of_face(j) /= -1 .or. isflowside3(j) > 0) then ! note: land = -1; water =  0; open boundary = 1 ~ n
					end do ! do l=1,tri_or_quad(i)
					
					! third, combine vertical terms and horizontal terms for the r_jk term:
					do i3=1,maxtran2
						sum1_all(i3) = sum1_all(i3) + sum1(i3)
						sum2_all(i3) = sum2_all(i3) + sum2(i3)
					end do
					
					! finally, calculate r_jk at only outflow faces in this element
					do l=1,tri_or_quad(i)
						j = facenum_at_cell(l,i)
						ii = adj_cellnum_at_cell(l,i)
						
						if(boundary_type_of_face(j) /= -1 .or. isflowside3(j) > 0 .or. isflowside4(j) > 0) then ! note: land = -1; water =  0;open boundary = 1 ~ n
							! outflow through this face, 
							! Note: Q_jk_theta2(l) is already calculated at the previous "do loop", 
							! and that is why I don't need to calculate it again here.
							if(Q_jk_theta2(l) >= 0.0_dp) then ! to include zero flow faces also...
								if(ii > 0) then	! if adjacent element exists (water cell exists)
									do i3=1,maxtran2
										con_cell_iik(i3) = con_cell(k,ii,i3)
									end do
								else ! if adjacent element doesn't exist (i.e., it is a boundary face)
									do i3=1,maxtran2
										con_cell_iik(i3) = con_cell_iik_oQb(k,i3) ! i.e., either con_at_obk(k,ibnd) or con_at_Qbk(k,ibnd)
									end do
								end if
								
								! now, we are ready to solve r_jk for each transport material:
								do i3=1,maxtran2
									rtemp = sum2_all(i3)*(con_cell_iik(i3) - con_cell(k,i,i3))
									if(rtemp == 0.0_dp) then ! this is to avoid "divide by 0"
										r_jk(k,j,i3) = 1.0 ! check this again
									else
										r_jk(k,j,i3) = sum1_all(i3)/rtemp
									end if
								end do
							end if
						end if
					end do ! do l=1,tri_or_quad(i)
				end do ! do k=b_layer,t_layer
				! end of r_jk calculation -------------------------------------------!
				
				! calculate r_ik ----------------------------------------------------!
				! I will calculate r_ik excluding at r_ik(m-1,,) & r_ik(M,,)
				! r_ik(0:m-1,,) & r_ik(M:maxlayer,,) are set to 1.0 from the initialization process.
				if(Q_ik_theta(t_layer,i) >= 0.0_dp) then ! w > 0
					do k=b_layer,t_layer-1
						do i3=1,maxtran2
							rtemp1 = sum1(i3) + abs(Q_ik_theta(k-1,i)) * (con_cell(k,i,i3)-con_cell(k-1,i,i3))
							rtemp2 = sum2(i3) + abs(Q_ik_theta(k-1,i))
							rtemp  = (con_cell(k+1,i,i3)-con_cell(k,i,i3))*rtemp2
							if(rtemp == 0.0) then ! this is to avoid "divide by 0"
								r_ik(k,i,i3) = 1.0_dp
							else
								r_ik(k,i,i3) = rtemp1/rtemp
							end if						
						end do
					end do
					do i3=1,maxtran2
						r_ik(t_layer  ,i,i3) = 1.0_dp
						r_ik(b_layer-1,i,i3) = 1.0_dp
					end do
				else ! w < 0
					do k=b_layer+1,t_layer
						do i3=1,maxtran2
							rtemp1 = sum1(i3) + abs(Q_ik_theta(k,i)) * (con_cell(k,i,i3)-con_cell(k+1,i,i3))
							rtemp2 = sum2(i3) + abs(Q_ik_theta(k,i))
							rtemp  = (con_cell(k-1,i,i3)-con_cell(k,i,i3))*rtemp2
							if(rtemp == 0.0_dp) then ! this is to avoid "divide by 0"
								r_ik(k-1,i,i3) = 1.0_dp
							else
								r_ik(k-1,i,i3) = rtemp1/rtemp
							end if
						end do
					end do
					do i3=1,maxtran2
						r_ik(t_layer  ,i,i3) = 1.0_dp
						r_ik(b_layer-1,i,i3) = 1.0_dp
					end do
				end if
				! end of r_ik calculation -------------------------------------------!
			end do ! do i=1,maxele
			!$omp end do
			
			! (4) calculate: large_PHI_ik & PSI_ik:
			!$omp do private(i,k,i3,b_layer,t_layer)
				do i=1,maxele
					b_layer = bottom_layer_at_element(i)
					t_layer = top_layer_at_element(i)
					
					do k=b_layer,t_layer-1 ! mid-levels only (i.e., excluding at m-1 & M levels)
						do i3=1,maxtran2
							if(h_flux_limiter == 0) then
								PSI_ik(k,i,i3) = 0.0_dp
							else if(h_flux_limiter == 1) then ! Minmod
								large_PHI_ik(k,i,i3) = max(small_phi_ik(k,i), &
								&									min(1.0,r_ik(k,i,i3)))
								PSI_ik(k,i,i3) = large_PHI_ik(k,i,i3) - small_phi_ik(k,i)
							else if(h_flux_limiter == 2) then ! van Leer
								large_PHI_ik(k,i,i3) = max(small_phi_ik(k,i), &
								&									(r_ik(k,i,i3) + abs(r_ik(k,i,i3)))/(1.0+abs(r_ik(k,i,i3))))
								PSI_ik(k,i,i3) = large_PHI_ik(k,i,i3) - small_phi_ik(k,i)
							else if(h_flux_limiter == 3) then ! Superbee
								large_PHI_ik(k,i,i3) = max(small_phi_ik(k,i), &
								&									min(1.0,2.0*r_ik(k,i,i3)), &
								&									min(2.0,r_ik(k,i,i3)))
								PSI_ik(k,i,i3) = large_PHI_ik(k,i,i3) - small_phi_ik(k,i)
							end if
						end do
					end do
				end do
			!$omp end do
			
			! (5) calculate large_PHI_jk & PSI_jk:
			!$omp do private(j,k,i3,b_layer,t_layer)
				do j=1,maxface
					b_layer = bottom_layer_at_face(j)
					t_layer = top_layer_at_face(j)
					
					do k=b_layer,t_layer
						do i3=1,maxtran2
							if(h_flux_limiter == 0) then
								PSI_jk(k,j,i3) = 0.0_dp
							else if(h_flux_limiter == 1) then ! Minmod
								large_PHI_jk(k,j,i3) = max(small_phi_jk(k,j), &
								&									min(1.0,r_jk(k,j,i3)))
								PSI_jk(k,j,i3) = large_PHI_jk(k,j,i3) - small_phi_jk(k,j)
							else if(h_flux_limiter == 2) then ! van Leer
								large_PHI_jk(k,j,i3) = max(small_phi_jk(k,j), &
								&									(r_jk(k,j,i3) + abs(r_jk(k,j,i3)))/(1.0+abs(r_jk(k,j,i3))))
								PSI_jk(k,j,i3) = large_PHI_jk(k,j,i3) - small_phi_jk(k,j)
							else if(h_flux_limiter == 3) then ! Superbee
								large_PHI_jk(k,j,i3) = max(small_phi_jk(k,j), &
								&									min(1.0,2.0*r_jk(k,j,i3)), &
								&									min(2.0,r_jk(k,j,i3)))
								PSI_jk(k,j,i3) = large_PHI_jk(k,j,i3) - small_phi_jk(k,j)
							end if
						end do
					end do
				end do
			!$omp end do
			!$omp end parallel
			! End of pre-calculate flux limiter function related variables ======!
		end if ! if(h_flux_limiter > 0 .or. v_flux_liminter > 0) then


		! ======================================================================!
		! main transport calculation start =====================================!
		!$omp parallel
		!$omp do private(i,j,k,l,ii,kk,i3, 						&
		!$omp & 			b_layer,t_layer,num_vertical_layer, &
		!$omp &			dz_old, dz_new,							&
		!$omp &			Q_jk_theta2, 								&
		!$omp & 			ith_ob, ith_Qb, ith_WR, ith_SS,		&
		!$omp &			con_cell_iik,							 	&
		!$omp &			con_cell_iik_oQb,							&
		!$omp &			sum_QC_i,	sum_QC_ii,					&
		!$omp &			sum_hdiff, sum_hlimiter, 				&
		!$omp &			h_adv, h_diff, h_flux_term, v_adv, v_flux_term, &
		!$omp &			a_lower_mat,b_diagonal_mat,c_upper_mat,rhs,solution,gam)
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			! Skip transport equation for dry element
			! If it is a dry element, do not calculate transport equation.
			if(t_layer == 0) then
				cycle ! do nothing and go to the next element
			end if
			
			! For all other wet element solve transport equation ================!
			! Now, we have to solve transport equation for all other elements at each vertical layer.
			num_vertical_layer = t_layer - b_layer + 1 ! number of active vertical layers at this element
			
			! First, setup interior elements of matrix A: =======================!
			! 		note for numbering for vertical direction
			!     matrix number is reversed for top and bottom layer
			!     top layer is always 1 and bottom layer is equal to actual surface layer number
			!     i.e. bottom layer = actual vertical layer number, 1 = actual top layer
			!     this is for just convinience for solving tridiagonal matrix
			
			! -------------------------------------------------------------------!
			! setup diagonal elements of matrix A: a, b, & c
			! A = [b c 0 0 0 0 		-> surface layer
			!      a b c 0 0 0 
			!      0 a b c 0 0
			!      - - - - - -
			!      0 0 0 a b c
			!      0 0 0 0 a b]		-> bottom layer
			! -------------------------------------------------------------------!
			! Note, we are constructing this matrix at vertical layers not vertical levels.
			
			! (1) setup matrix A for mid-layers: ================================!
			! From here, I will use "kk" for matrix ordering, and "k" for corresponding variable's layer ordering (i.e, upside-down ordering).
			! Note: 
			! 		"A" matrix is identical for all transport material.
			do k = b_layer+1, t_layer-1 ! "k" is the index of variable system
				! First, setup A matrix at the middle of the water column (mid-layers)
				! If it is a 2D (i.e., one layer) or if it is a 3D but with only two layers, 
				! this part will not be performed since it only has top & bottom layers (not mid-layers)
				kk = t_layer + 1 - k ! "kk" is the index of matrix
				! i.e., from here "kk" is the matrix index in Ax = b system
				!                 "k" is the actual variable index
				
				! calculate left-hand-side matrix, A: ----------------------------!
				! note: this part is different from non-subcycling
				a_lower_mat(kk) = -dt2*d_ik_theta(k  ,i)
				c_upper_mat(kk) = -dt2*d_ik_theta(k-1,i)
				! Note:
				!      Here, dz_new is calculated for the current subcycle timestep using (subcycle+1).
				dz_new = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle+1)/Nt
				b_diagonal_mat(kk) = - a_lower_mat(kk) &
				&							+ area(i)*dz_new &
				&							- c_upper_mat(kk)
			end do ! do kk = bottom_layer_at_element(i)+1, top_layer_at_element(i)-1 ; end mid-layer
			
			
			! (2) setup the top and bottom layer's elements of matrix A =========!
			if(t_layer == b_layer) then ! if it is 2D (i.e., one layer)
				! We only has one layer and the layer number is "one"
				! to keep "k" index in the following equations, I just define "k" to 1.
				! Otherwise, I have to change "k" to "1" in following equations.
				! For the one layer system, matrix index is equal to variable index
				kk = 1 ! this is the matrix index
				k = t_layer ! this is the variable index
				
				! if it is 2D, we only need the diagonal matrix term, b_diagonal_mat, at the left-hand-side matrix
				! i.e., no vertical diffusion terms (between top and current layer, and between bottom and current layer) are required
				! a_lower_mat(kk) = 0.0_dp
				! note: this part is different from non-subcycling
				dz_new = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle+1)/Nt
				b_diagonal_mat(kk) = area(i)*dz_new
				! c_upper_mat(kk) = 0.0_dp
			else ! if it is with only 2 layers (i.e., top and bottom layer only) or full 3D
				! (1) at the top(surface) layer ----------------------------------!
				! Only b and c terms are required
				! i.e., only vertical advection and vertical diffusion between current layer and the bottom layer
				! No vertical advection and diffusion between current layer and the surface (since no above layer exists)
				kk = 1
				k = t_layer ! k should be 2; we are in the layer #2
				
				! no vertical diffusion between current layer and above layer (i.e., no more above layer exists)
				! only diffusion between current layer and bottom layer exists
				! note: this part is different from non-subcycling
				! a_lower_mat(kk) = 0.0_dp
				c_upper_mat(kk) = -dt2*d_ik_theta(k-1,i)
				dz_new = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle+1)/Nt
				b_diagonal_mat(kk) = area(i)*dz_new &
				&						  -c_upper_mat(kk)
				
				! (2) at the bottom layer ----------------------------------------!
				! Only a and b terms are requried
				kk = num_vertical_layer ! now, we are in the bottom row in a matrix
				k = b_layer
				
				! now turn of vertical advection and diffusion between current layer and 0th layer (we do not have lower layer)
	        	! no vertical diffusion between current layer and bottom layer 
	        	! only diffusion between current layer and top layer exists
	      	! note: this part is different from non-subcycling
	         a_lower_mat(kk) = -dt2*d_ik_theta(k,i)
	         dz_new = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle+1)/Nt
	         b_diagonal_mat(kk) = -a_lower_mat(kk) &
	         &						   + area(i)*dz_new
	         ! c_upper_mat(kk) = 0.0_dp
	      end if ! if(top_layer_at_element(i) == bottom_layer_at_element(i)) then ! if it is 2D (i.e., one layer)      
	      ! End of A matrix ===================================================!

			! Now, setup right-hand side vector "b" =============================!
			! First, include horizontal terms: ==================================!
			! (1) horizontal advection, (2) horizontal diffusion, and (3) horizontal flux limiter term
			! 		Note that: at all layers, calculation of horizontal terms are identical, and
			! 		this routine should search through each side (face) in the current element.
			! Then, I will include vertical terms in "b" 			
			! Note: 
			! 		"b" vectors are different for all transport material.

			! First, if this element belongs to tidal open boundary or river boundary,
			! copy boundary concentration value to the concentration at the neighbor element 
			ith_ob = ob_element_flag(i) ! This is the ith (i.e., ibnd'th) open boundary
			ith_Qb = Qb_element_flag(i) ! This is the ith river boundary
			ith_WR = WR_element_flag(i) ! This is the ith withdraw/return boundary
			if(ith_ob > 0) then ! if this element belongs to tidal open boundary					
	 			do k=b_layer,t_layer
					do i3=1,maxtran2
						con_cell_iik_oQb(k,i3) = con_at_obk(k,ith_ob,i3) ! salt, temp, ...
					end do
				end do
			end if
			if(ith_Qb > 0) then ! if this element belongs to river boundary
				do k=b_layer,t_layer
					do i3=1,maxtran2
						con_cell_iik_oQb(k,i3) = con_at_Qbk(k,ith_Qb,i3) ! salt, temp, ...
						! write(*,*) 'con_cell_iik_oQb(k,i3) = ', i, k, i3, con_cell_iik_oQb(k,i3)
					end do
				end do
			end if
			if(ith_WR > 0) then ! if this element belongs to withdraw/return boundary
				! i2 = WR_boundary(i,1)! element ID (i.e., element number)
				j = WR_boundary(i,2) ! face number of the WR boundary
				if(WR_layer(ith_WR) == 999) then
					! anytime surface layer
					k = top_layer_at_face(j)
				else if(WR_layer(ith_WR) == 0) then
					! anytime bottom layer
					k = bottom_layer_at_face(j)
				else
					! specific vertical layer
					k = WR_layer(ith_WR)
				end if
					
				do i3=1,maxtran2
					con_cell_iik_oQb(k,i3) = con_at_WRk(k,ith_WR,i3)
					! write(*,*) 'con_cell_iik_oQb(k,i3) = ', k, i3, con_cell_iik_oQb(k,i3)
				end do				
			end if

			do k=b_layer,t_layer
				kk = t_layer + 1 - k ! "kk" is the index of matrix
				
	         sum_QC_i = 0.0_dp ! sum of QC for leaving face (i.e., QC_out), for horizontal advection
	         sum_QC_ii = 0.0_dp ! sum of QC for entering face(i.e., QC_in), for horizontal advection
	         sum_hdiff = 0.0_dp ! sum for horizontal diffusion
	         sum_hlimiter = 0.0_dp ! sum for horizontal adv/diff for flux limiter
				do l=1,tri_or_quad(i) ! calculate at each horizontal face
					j = facenum_at_cell(l,i)
					ii = adj_cellnum_at_cell(l,i) ! ii = neighbor cell number which shareing this side with the current cell

	        		if(ii > 0) then	! if adjacent element exists (water cell exists)
	        			do i3=1,maxtran2
	        				con_cell_iik(i3) = con_cell(k,ii,i3)
	        			end do
	        		else ! if adjacent element doesn't exist (i.e., it is a boundary face)
	        			do i3=1,maxtran2
	        				con_cell_iik(i3) = con_cell_iik_oQb(k,i3)
	        			end do
	        		end if

					! calculate horizontal advection and diffusion, only if adjacent cell, which shares this side, exists.
					! otherwise, i.e., if this side locates at the land boundary, horizontal advection and diffusion will be "0.0",
					! and if this side locates at the open boundary, i.e., this element is the boundary element,
					! the concentration is already defined with the boundary value.
					! if(boundary_type_of_face(j) == 0) then
					if(boundary_type_of_face(j) /= -1  .or. isflowside3(j) > 0 .or. isflowside4(j) > 0) then ! if this side is not a land face (i.e., face to water or open boundary) or it is the river boundary face
						! for horizontal advection through j'th side
			      	! note: this part is different from non-subcycling
			      	Q_jk_theta2(l) = Q_jk_theta(k,j)*sign_in_outflow(l,i)

	 					! Calculate horizontal advection at this face --------------!
	 					! if(sign_in_outflow(l,i)*un_face(k,j) >= 0.0) then ! this is the original approach, but I used the next approach since:
	 					! what if un_face and un_face_new has differrent direction?; i.e., when flow direction is changed.
	 					! That is why I think using Q_jk_theta2 (as shown above) might be better.
	 					! note: if Q_jk_theta2 == 0.0, there will be no advection, so both "sum_QC_i" and "sum_QC_ii" will not be changed from the previous value.
						if(Q_jk_theta2(l) > 0.0_dp) then
							! if flow leaves through this face from the current element, i,
							! then concentration moves from the current cell to the neighbor element.
							do i3=1,maxtran2
								sum_QC_i(i3) = sum_QC_i(i3) +  abs(Q_jk_theta2(l))*con_cell(k,i,i3) ! QC_out; advection for leaving face
							end do
						else if(Q_jk_theta2(l) < 0.0_dp) then ! inflow through this face
							! if flow enters through this face from the neighbor element, ii, to the current element, i,
							! then concentration moves from the neighbor cell to the current element.
							do i3=1,maxtran2
								sum_QC_ii(i3) = sum_QC_ii(i3) + abs(Q_jk_theta2(l))*con_cell_iik(i3) ! QC_in; advection for entering face
							end do
						end if ! if(Q_jk_theta2(l) > 0.0) then
						
						! Now, let's calculate horizontal diffusion through j'th side (this side)
	 					do i3=1,maxtran2
	 						! this tells horizontal diffusion occurs only if concentration difference exists between two cells.
	 						sum_hdiff(i3) = sum_hdiff(i3) + d_jk_theta(k,j)*(con_cell_iik(i3) - con_cell(k,i,i3))
	 					end do

						! calculate flux limiter term at each horizontal side ------!
						do i3=1,maxtran2
		 					sum_hlimiter(i3) = sum_hlimiter(i3) + &
		 					&						 PSI_jk(k,j,i3)*abs(Q_jk_theta2(l))*(con_cell(k,i,i3) - con_cell_iik(i3))
		 				end do
					end if ! if(boundary_type_of_face(j) /= -1) then
				end do ! do l=1,tri_or_quad(i)
	         
	         ! now, include horizontal terms in vector "b":
	         do i3=1,maxtran2
		         ! note that this part is different from non-subcycling
		         ! 		dt -> dt2
		         ! (1) horizontal advection term through all sides
		         h_adv(i3) = dt2*(sum_QC_ii(i3) - sum_QC_i(i3)) ! (QC_in - QC_out)
		         
		         ! (2) horizontal diffusion term through all sides
		         h_diff(i3) = dt2*sum_hdiff(i3)
		         
		         ! (3) horizontal flux limiter term
		         h_flux_term(i3) = 0.5*dt2*sum_hlimiter(i3)
	         
	         
		         ! update right-hand side vector including:
		         ! Note: this part is different from no-subcycling
		         ! 	 (0) explicit term
		         ! + (1) horizontal advection
		         ! + (2) horizontal diffusion
		         ! + (3) horizontal flux limiter
		         ! Note: vertical diffusion term is in the left hand side, i.e., in Matrix A
		         !       and vertical advection terms will be included later

	         	! do not include flux limiter correction
	         	! actually, I do not need to do this since PSI_jk is already 0.0 if "h_flux_limiter == 0",
	         	! thus sum_hlimiter = 0 & h_flux_term = 0 for "h_flux_limiter == 0.
	         	! however, I just use this approach to make sure and to clearly show the equation.
	         	! note: this part is different from non-subcycling
	         	dz_old = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle)/Nt
	         	rhs(kk,i3) = area(i)*dz_old*con_cell(k,i,i3) &
	         	&				+ h_adv(i3) + h_diff(i3)

		         if(h_flux_limiter > 0) then
		         	rhs(kk,i3) = rhs(kk,i3) + h_flux_term(i3)
		         end if
		      end do ! do i3=1,maxtran2
		   end do ! do k=b_layer,t_layer
		   ! End of including horizontal terms in vector "b" ===================!
			
			
			! Second, let's include vertical advection terms in "b" =============!
			! Note: vertical diffusion term is already included in the left-hand side matrix.
			!       so, I just need to include vertical advection and vertical flux limter terms in the righ-hand side vector.
			! Note: at the top and bottom boundary layers, the vertical advection term is different from mid-layers.
			if(num_vertical_layer == 1) then
				! if this water column has only one vertical layer (i.e., 2D), we don't need to include vertical advection & vertical flux terms.
				! so, do nothing (i.e., do not add vertical advection terms in the "rhs").
			else
				do k=b_layer,t_layer
					kk = t_layer + 1 - k ! "kk" is the index of matrix
					
					if(k == b_layer) then ! at the bottom layer
						do i3=1,maxtran2
							! if we have positive vertical velocity, we will have outflux to upper layer
							! if we have negative vertical velocity, we will have influx from upper layer
							! note: this part is different from non-subcycling:
							! 		dt -> dt2
							! this is the original approach when source/sink was not included
							v_adv(i3) = dt2*( &
							&	&
							&	-abs(0.5*(Q_ik_theta(k,i) + abs(Q_ik_theta(k,i))))*con_cell(k  ,i,i3) &
							&	+abs(0.5*(Q_ik_theta(k,i) - abs(Q_ik_theta(k,i))))*con_cell(k+1,i,i3) &
							&	)
														
							rhs(kk,i3) = rhs(kk,i3) + v_adv(i3)
							
							if(v_flux_limiter > 0) then
								v_flux_term(i3) = 0.5*dt2*( &
								& &
								& -PSI_ik(k,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k,i) + abs(Q_ik_theta(k,i))))*(con_cell(k+1,i,i3)-con_cell(k  ,i,i3)) &
								& +PSI_ik(k,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k,i) - abs(Q_ik_theta(k,i))))*(con_cell(k  ,i,i3)-con_cell(k+1,i,i3)) &
								&  )
								
								rhs(kk,i3) = rhs(kk,i3) + v_flux_term(i3)
							end if
						end do
					else if(k == t_layer) then ! at the top layer
						do i3=1,maxtran2
							! if we have positive vertical velocity, we will have outflux to upper layer
							! if we have negative vertical velocity, we will have influx from upper layer
							! note: this part is different from non-subcycling:
							! 		dt -> dt2								
							v_adv(i3) = dt2*( &
							&	 abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))*con_cell(k-1,i,i3) &
							&	&
							&	&
							&	-abs(0.5*(Q_ik_theta(k-1,i) - abs(Q_ik_theta(k-1,i))))*con_cell(k  ,i,i3) )
														
							rhs(kk,i3) = rhs(kk,i3) + v_adv(i3)
							
							if(v_flux_limiter > 0) then
								v_flux_term(i3) = 0.5*dt2*( &
								&  PSI_ik(k-1,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))*(con_cell(k  ,i,i3)-con_cell(k-1,i,i3)) &
								&  &
								&  &
								& -PSI_ik(k-1,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k-1,i) - abs(Q_ik_theta(k-1,i))))*(con_cell(k-1,i,i3)-con_cell(k  ,i,i3)) )
								
								rhs(kk,i3) = rhs(kk,i3) + v_flux_term(i3)
							end if
						end do
					else ! at the mid-layers
						do i3=1,maxtran2
							! if we have positive vertical velocity, we will have outflux to upper layer
							! if we have negative vertical velocity, we will have influx from upper layer
							! note: this part is different from non-subcycling:
							! 		dt -> dt2								
							v_adv(i3) = dt2*( &
							&	 abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))*con_cell(k-1,i,i3) &
							&	-abs(0.5*(Q_ik_theta(k  ,i) + abs(Q_ik_theta(k  ,i))))*con_cell(k  ,i,i3) &
							&	+abs(0.5*(Q_ik_theta(k  ,i) - abs(Q_ik_theta(k  ,i))))*con_cell(k+1,i,i3) &
							&	-abs(0.5*(Q_ik_theta(k-1,i) - abs(Q_ik_theta(k-1,i))))*con_cell(k  ,i,i3) )
														
							rhs(kk,i3) = rhs(kk,i3) + v_adv(i3)
							
							if(v_flux_limiter > 0) then
								v_flux_term(i3) = 0.5*dt2*( &
								&	PSI_ik(k-1,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))*(con_cell(k  ,i,i3)-con_cell(k-1,i,i3)) &
								& -PSI_ik(k  ,i,i3)*&
								&		abs(0.5*(Q_ik_theta(  k,i) + abs(Q_ik_theta(k  ,i))))*(con_cell(k+1,i,i3)-con_cell(k  ,i,i3)) &
								& +PSI_ik(k  ,i,i3)*&
								&		abs(0.5*(Q_ik_theta(  k,i) - abs(Q_ik_theta(k  ,i))))*(con_cell(k  ,i,i3)-con_cell(k+1,i,i3)) &
								& -PSI_ik(k-1,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k-1,i) - abs(Q_ik_theta(k-1,i))))*(con_cell(k-1,i,i3)-con_cell(k  ,i,i3)) )
								
								rhs(kk,i3) = rhs(kk,i3) + v_flux_term(i3)	
							end if							
						end do
					end if ! if(k == b_layer) then ! at the bottom layer
				end do ! do k=b_layer,t_layer
			end if ! if(num_vertical_layer == 1) then
			! End of including vertical advection & limiter terms in the rhs ====!
			
	      			
			! Now, get new solutions in the current element (water column) ======!
	      ! So, call Thomas Algorithm
	      ! Note, I only have one vector for the right hand side vector
			! gam & solution are output results, we don't need gam, which is the pivot value for tridiagonal matrix solver
			! solution(:,1) = is our solution (i.e., concentration at the center of each vertical layer in this element
	      ! call tridiagonal_solver(maxlayer+1, num_vertical_layer, 1, a_lower_mat, b_diagonal_mat, c_upper_mat, rhs, solution, gam)
	      call tridiagonal_solver(maxlayer, num_vertical_layer, maxtran2, a_lower_mat, b_diagonal_mat, c_upper_mat, rhs, solution, gam)
	      
	      do k = b_layer,t_layer
	      	kk = t_layer - k + 1 ! this is the corresponding matrix index
	      	
	      	! There could be (slight) negative solution due to (roundoff) numerical errors, and that negative values cause an error.
	      	! Thus, to suppress the error, let's set the minimum value to 0.0 and maximum value to 50.0, both for salt & temp
	      	do i3=1,maxtran2
	      		if(tran_id(i3) == 1 .or. tran_id(i3) == 2) then ! both for salt & temp
			      	if(solution(kk,i3) >= 0.0_dp .and. solution(kk,i3) <= 50.0_dp) then
			      		con_cell_new(k,i,i3) = solution(kk,i3)
			      	else if(solution(kk,i3) < 0.0_dp) then
			      		con_cell_new(k,i,i3) = 0.0_dp
			      	else if(solution(kk,i3) > 50.0_dp) then
			      		con_cell_new(k,i,i3) = 50.0_dp
			      	end if
			      end if
		      end do
	      end do
	      
	      ! let's include source/sink here... 
	      ! this is explicit inclusion test...
	      if(SS_element_flag(i) > 0) then
	      	ith_SS = SS_element_flag(i)
	      	k = SS_layer(ith_SS)
	      	if(Q_add_SS(ith_SS) >= 0) then
	      		! source
		      	do i3=1,maxtran2
			      	if(k == 999) then
			      		! at surface
			      		con_cell_new(t_layer,i,i3) = con_cell_new(t_layer,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_at_SS(ith_SS,i3)*dt2/(area(i)*dz_cell_new(t_layer,i))
			      	else if(k == 0) then
			      		! at bottom
			      		con_cell_new(b_layer,i,i3) = con_cell_new(b_layer,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_at_SS(ith_SS,i3)*dt2/(area(i)*dz_cell_new(b_layer,i))
			      	else
			      		! at any layer
			      		con_cell_new(k,i,i3) = con_cell_new(k,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_at_SS(ith_SS,i3)*dt2/(area(i)*dz_cell_new(k,i))	
			      	end if
			      end do
			   else
			   	! sink
			   	do i3=1,maxtran2
			      	if(k == 999) then
			      		! at surface
			      		con_cell_new(t_layer,i,i3) = con_cell_new(t_layer,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_cell(t_layer,i,i3)*dt2/(area(i)*dz_cell_new(t_layer,i))
			      	else if(k == 0) then
			      		! at bottom
			      		con_cell_new(b_layer,i,i3) = con_cell_new(b_layer,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_cell(b_layer,i,i3)*dt2/(area(i)*dz_cell_new(b_layer,i))
			      	else
			      		! at any layer
			      		con_cell_new(k,i,i3) = con_cell_new(k,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_cell(k,i,i3)*dt2/(area(i)*dz_cell_new(k,i))	
			      	end if
			      end do
		      end if
	      end if
	   end do ! do i=1,maxele
		!$omp end do
		! End of transport equation calculation for each water column ==========!
	   
	   ! note: this part is different from non-subcycling
	   ! update old variable to new variable for the next subcycle
 	   !$omp do private(i,k,i3)
		do i=1,maxele
			do k=bottom_layer_at_element(i),top_layer_at_element(i)
		   	do i3=1,maxtran2
		   		con_cell(k,i,i3) = con_cell_new(k,i,i3)
		   	end do
		   end do
	 		
	 		! extend to top/bottom ghost cells:
	 		do k=0,bottom_layer_at_element(i)-1
	 			do i3=1,maxtran2
					con_cell(k,i,i3) = 0.0_dp
				end do
			end do
			do k=top_layer_at_element(i)+1,maxlayer+1
	 			do i3=1,maxtran2
					con_cell(k,i,i3) = 0.0_dp
				end do				
			end do
		end do
 		!$omp end do
 	   !$omp end parallel
 	   ! End of OMP parallel region
 	   
 	   ! Enforce open boundary condition to the boundary cells:
 	   ! I am not sure if this part is needed or not.
! 		do i=1,num_ob_cell
! 			! currently, we have four open boundary surface elevation types:
! 			! 		ob_eta_type(i) == -1, 1, 2, 3
! 			! 		-1: radiation boundary condition
! 			!  	 1: harmonic tide (cosine wave)
! 			! 		 2: eta_ser.inp
! 			ie = ob_cell_id(i) ! open boundary element ID
! 			! First, prepare the boundary concentration from the timeseries data
! 			do i3=1,maxtran2
! 				! for salinity:
! 				if(tran_id(i3) == 1 .and. salt_ser_id(i) > 0) then
! 					! salinity at ith open boundary element (at corresponding ghost cell)
! 					do k=1,maxlayer
! 						con_cell_new(k,ie,i3) = con_at_obk(k,i,i3)
! 						con_cell(k,ie,i3) = con_at_obk(k,i,i3)
! 					end do
! 				end if
! 			end do
! 		end do
		! End of transport equation at cell center =============================!
	end do ! do subcycle=0,(Nt-1) ! subcycle starts
	
	
	! Now, calculate concentrations at node and face ==========================!
	! (1) calculate concentration at node 
	!$omp parallel
	!$omp do private(i,k,l,ii,i3,sum11,sum22,count1)
	do i=1,maxnod
		! note, here I am doing from top to bottom, and there is reason
		! previously, there was an unintended node value (i.e., 0.0) could be asigned when all adjacent cells are located at the next vertical layer than this node.
		! To prevent this error, I updated this part as follows:
		! If all adjacent cells are located above this node's vertical level, 
		! just use the top node value for the current nodal value.
		do k=top_layer_at_node(i),bottom_layer_at_node(i),-1
			do i3=1,maxtran2
				sum11 = 0.0_dp ! total volume
				sum22 = 0.0_dp ! total mass
				count1 = 0 ! count adjacent cells which are in the same vertical layer with this node
				do l = 1,adj_cells_at_node(i)
					ii = adj_cellnum_at_node(l,i)
					
					! if(dz_cell_new(k,ii) > 0.0_dp) then
					if(k>=bottom_layer_at_element(ii) .and. k<=top_layer_at_element(ii)) then ! this is better criteria than above...
						sum11 = sum11 + area(ii)*dz_cell_new(k,ii) ! total volume
						sum22 = sum22 + area(ii)*dz_cell_new(k,ii) * con_cell_new(k,ii,i3) ! total mass
					else
						count1 = count1 + 1
					end if
				end do

				if(sum11 == 0.0_dp) then
					con_node_new(k,i,i3) = 0.0_dp
				else
					con_node_new(k,i,i3) = sum22/sum11 ! [mass]/[m^3], volumetric weighting average
				end if
				
				! if all adjacent cells are located above this node's vertical level,
				! just use the top node value for the current nodal value.
				if(count1 == adj_cells_at_node(i)) then
					if(k /= top_layer_at_node(i)) then
						con_node_new(k,i,i3) = con_node_new(k+1,i,i3)
					end if
				end if				
			end do	
		end do
		
! 		do i3=1,maxtran2
! 			if(con_node_new(top_layer_at_node(i),i,i3) > 0.0_dp) then
! 				do k=top_layer_at_node(i)-1,bottom_layer_at_node(i),-1
! 					if(con_node_new(k,i,i3) == 0.0_dp) then
! 						con_node_new(k,i,i3) = con_node_new(k+1,i,i3)
! 					end if
! 				end do			
! 			end if
! 		end do
	end do
	!$omp end do

	! (2) calculate concentration at face
	!$omp do private(j,n1,n2,k,i3)
	do j=1,maxface
		n1 = nodenum_at_face(1,j)
		n2 = nodenum_at_face(2,j)
		do k=bottom_layer_at_face(j),top_layer_at_face(j)
			do i3=1,maxtran2
				con_face_new(k,j,i3) = (con_node_new(k,n1,i3) + con_node_new(k,n2,i3)) * 0.5
			end do
		end do
	end do
	!$omp end do
	!$omp end parallel

	! update "con" to corresponding variable ==================================!
	do i3=1,maxtran2
		if(tran_id(i3) == 1) then ! salinity
			salt_cell_new = con_cell_new(:,:,1)
			salt_node = con_node_new(:,:,1)
			salt_face = con_face_new(:,:,1)
		end if
		if(tran_id(i3) == 2) then ! temperature
			temp_cell_new = con_cell_new(:,:,2)
			temp_node = con_node_new(:,:,2)
			temp_face = con_face_new(:,:,2)
		end if
	end do
	! End of concentration calculation at node and face =======================!

! 	do i=1,maxele
! 		write(*,'(A,I3,*(F5.2))') 'i=', i, (salt_cell_new(k,i), k=1,10)
! 	end do
! 	write(*,*)
	! stop 'I am here'
end subroutine solve_transport_equation_v19