!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will calculate eta at the ghost cell at boundary elements at both previous and current time step
!! This explicitly given boundary elevations will be used at:
!! 		solve_momentum_equation.f90
!! 		solve_free_surface_equation.f90
!! 
subroutine compute_boundary_eta
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, j, l
	integer :: nd, i2, ie, temp_element

	! interpolation variables
	integer :: t2
	real(dp):: u1,u3,v1,v3
	real(dp):: u2_old, u2_new, v2_old, v2_new
	real(dp):: sum0
	integer :: icount
	! End of local variables ==================================================!
	
   ! =========================================================================!
   ! First, we have to prepare the surface elevation at open boundary elements,
   ! which is externally given, i.e., known value, for both previous and current time step.
   ! =========================================================================!
   ! u2 is the julian day for the previous time step, and it will be used for interpolation.
	u2_old = julian_day*86400.0_dp - dt ! julian time [s], for the previous time and/or time for looking (since we are in the next time level, so we have subtract dt)
	u2_new = julian_day*86400.0_dp 		! julian time [s], for the current time

   do i=1,num_ob_cell
      temp_element = ob_cell_id(i)
		
		! Note: we have three open boundary surface elevation types:
		! 		ob_eta_type(i) == -1, 1, 2
		! 		-1: radiation boundary condition
		!  	 1: harmonic tide (cosine wave)
		! 				The sine wave can be achieved with setting [tidal_phase_shift] = 90.0 in harmonic_ser.inp
		! 		 2: eta_ser.inp
		
		! (1) for the radiation boundary condition
		if(ob_eta_type(i) == -1) then
			sum0 = 0.0_dp
			icount = 0
			do l = 1, tri_or_quad(temp_element)
				nd = nodenum_at_cell(l,temp_element)
				do i2 = 1,adj_cells_at_node(nd)
					ie = adj_cellnum_at_node(i2,nd)
					
					! ob_element_flag(i): 0 		(non open boundary element)
					! 						  : 1 ~ n 	(open boundary element)
					if(ob_element_flag(ie) == 0) then 
						icount = icount + 1
						sum0 = sum0 + eta_cell(ie)
					end if
				end do
			end do 
         
			! check correctness of open boundary 
			if(icount == 0) then
				write(pw_run_log,*)'isolated obe cannot have ob condition', temp_element
				stop
			else
				! eta_at_ob_old(temp_element) = sum0/icount
				eta_at_ob_old(i) = sum0/icount
			end if   
		
		! (2) from harmonic tide (tidal constituents), cosine wave
      else if(ob_eta_type(i) == 1) then
      	! use traditional cosine wave
			eta_at_ob_old(i) = 0.0_dp
			eta_at_ob_new(i) = 0.0_dp
			do j = 1, no_tidal_constituent(harmonic_ser_id(i))
				! Note: There will be a slight error associated with the following calculation:
				! For example, at the first time step, "eta_at_ob_old" should be equal to 0.0, but this calculation is not exactly 0.0.
				! thus, it will create a slight error when calculating barotropic gradient term from the explicit term calculation in solve_momentum_equation.f90
				! However, it should be o.k. anyway since this is a minor numerical calculation error...
				eta_at_ob_old(i) = eta_at_ob_old(i) 						&
				&  + spinup_function_tide * tidal_amplitude(j,harmonic_ser_id(i)) &
				&	* tidal_nodal_factor(j,harmonic_ser_id(i))   						&
				&  * cos(2.0_dp*pi/(tidal_period(j,harmonic_ser_id(i))*3600.0_dp)*(elapsed_time-dt)	&
				&	+ equilibrium_argument(j,harmonic_ser_id(i))*deg2rad 				&
				&	- tidal_phase(j,harmonic_ser_id(i))*deg2rad - tidal_phase_shift(harmonic_ser_id(i))*deg2rad)
				
 				eta_at_ob_new(i) = eta_at_ob_new(i) 						&
 				&  + spinup_function_tide * tidal_amplitude(j,harmonic_ser_id(i)) &
 				&	* tidal_nodal_factor(j,harmonic_ser_id(i))   						&
 				&  * cos(2.0_dp*pi/(tidal_period(j,harmonic_ser_id(i))*3600.0_dp)*elapsed_time		&
 				&	+ equilibrium_argument(j,harmonic_ser_id(i))*deg2rad 				&
 				&	- tidal_phase(j,harmonic_ser_id(i))*deg2rad - tidal_phase_shift(harmonic_ser_id(i))*deg2rad) 				
			end do
		
		! (3) from time history of the elevation, eta_ser.inp			
      else if(ob_eta_type(i) == 2) then
      	! write(*,*) eta_ser_id(i), eta_ser_data_num(eta_ser_id(i))
      	
      	! for previous time step
			do t2=2,eta_ser_data_num(eta_ser_id(i)) ! 2 ~ 121
				! Since eta_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
				! u1, and u3 should be shifted toward the simulation start year.
				! That is why I put "- reference_diff_days * 86400.0"
				u1 = eta_ser_time(t2-1,eta_ser_id(i)) - reference_diff_days * 86400.0_dp  		! [s], lower bound time
				u3 = eta_ser_time(t2  ,eta_ser_id(i)) - reference_diff_days * 86400.0_dp		! [s], upper bound time
				
				if(u2_old >= u1 .and. u2_old <= u3) then
					v1 = eta_ser_eta(t2-1,eta_ser_id(i))	! [m], lower bound eta
					v3 = eta_ser_eta(t2,  eta_ser_id(i))	! [m], upper bound eta
					v2_old = (u2_old-u1)*(v3-v1)/(u3-u1)+v1	! interpolated eta
					exit
				end if			
			end do
      	eta_at_ob_old(i) = spinup_function_tide * v2_old
      	
      	! for current (new) time step
			do t2=2,eta_ser_data_num(eta_ser_id(i)) ! 2 ~ 121
				! Since eta_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
				! u1, and u3 should be shifted toward the simulation start year.
				! That is why I put "- reference_diff_days * 86400.0"
				u1 = eta_ser_time(t2-1,eta_ser_id(i)) - reference_diff_days * 86400.0_dp  		! [s], lower bound time
				u3 = eta_ser_time(t2  ,eta_ser_id(i)) - reference_diff_days * 86400.0_dp		! [s], upper bound time
				
				if(u2_new >= u1 .and. u2_new <= u3) then
					v1 = eta_ser_eta(t2-1,eta_ser_id(i))	! [m], lower bound eta
					v3 = eta_ser_eta(t2,  eta_ser_id(i))	! [m], upper bound eta
					v2_new = (u2_new-u1)*(v3-v1)/(u3-u1)+v1	! interpolated eta
					exit
				end if			
			end do
      	eta_at_ob_new(i) = spinup_function_tide * v2_new
      end if   !   if(ob_eta_type(k) == -1, 1, 2)
      
		! Add storm surge elevation from holland storm surge model
		! jw, it is not yet included...
		if(holland_flag == 1) then
			eta_at_ob_old(i) = eta_at_ob_old(i) + eta_from_Holland_at_ob(i) ! not yet corrected
			eta_at_ob_new(i) = eta_at_ob_new(i) + eta_from_Holland_at_ob(i) ! not yet corrected
		end if
   end do   ! i=1,num_ob_cell   
   ! End of boundary elevation preparation ===================================!
end subroutine compute_boundary_eta