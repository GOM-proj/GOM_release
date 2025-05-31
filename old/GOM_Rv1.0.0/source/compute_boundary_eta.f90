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

	! jw
	integer :: t2
	real(dp):: u1,u3,v1,v3
	real(dp):: u2_old, u2_new, v2_old, v2_new
	real(dp):: sum0
	integer :: icount
	! jw
	
   ! jw
   ! jw
   ! jw
   ! jw
   ! jw
	u2_old = julian_day*86400.0 - dt ! jw
	u2_new = julian_day*86400.0 		! jw

   do i=1,num_ob_cell
      temp_element = ob_cell_id(i)
		
		! jw
		! jw
		! jw
		! jw
		! jw
		! jw
		
		! jw
		if(ob_eta_type(i) == -1) then
			sum0 = 0.0_dp
			icount = 0
			do l = 1, tri_or_quad(temp_element)
				nd = nodenum_at_cell(l,temp_element)
				do i2 = 1,adj_cells_at_node(nd)
					ie = adj_cellnum_at_node(i2,nd)
					
					! jw
					! jw
					if(ob_element_flag(ie) == 0) then 
						icount = icount + 1
						sum0 = sum0 + eta_cell(ie)
					end if
				end do
			end do 
         
			! jw
			if(icount == 0) then
				write(pw_run_log,*)'isolated obe cannot have ob condition', temp_element
				stop
			else
				! jw
				eta_at_ob_old(i) = sum0/icount
			end if   
		
		! jw
      else if(ob_eta_type(i) == 1) then
      	! jw
      	
      	! jw
			do t2=2,eta_ser_data_num(eta_ser_id(i)) ! jw
				! jw
				! jw
				! jw
				u1 = eta_ser_time(t2-1,eta_ser_id(i)) - reference_diff_days * 86400.0  		! jw
				u3 = eta_ser_time(t2  ,eta_ser_id(i)) - reference_diff_days * 86400.0		! jw
				
				if(u2_old >= u1 .and. u2_old <= u3) then
					v1 = eta_ser_eta(t2-1,eta_ser_id(i))	! jw
					v3 = eta_ser_eta(t2,  eta_ser_id(i))	! jw
					v2_old = (u2_old-u1)*(v3-v1)/(u3-u1)+v1	! jw
					exit
				end if			
			end do
      	eta_at_ob_old(i) = spinup_function_tide * v2_old
      	
      	! jw
			do t2=2,eta_ser_data_num(eta_ser_id(i)) ! jw
				! jw
				! jw
				! jw
				u1 = eta_ser_time(t2-1,eta_ser_id(i)) - reference_diff_days * 86400.0  		! jw
				u3 = eta_ser_time(t2  ,eta_ser_id(i)) - reference_diff_days * 86400.0		! jw
				
				if(u2_new >= u1 .and. u2_new <= u3) then
					v1 = eta_ser_eta(t2-1,eta_ser_id(i))	! jw
					v3 = eta_ser_eta(t2,  eta_ser_id(i))	! jw
					v2_new = (u2_new-u1)*(v3-v1)/(u3-u1)+v1	! jw
					exit
				end if			
			end do
      	eta_at_ob_new(i) = spinup_function_tide * v2_new
      	
      	
		! jw
      else if(ob_eta_type(i) == 2) then
      	! jw
			eta_at_ob_old(i) = 0.0_dp
			eta_at_ob_new(i) = 0.0_dp
			do j = 1, no_tidal_constituent(harmonic_ser_id(i))
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
      end if   ! jw
      
		! jw
		! jw
		if(holland_flag == 1) then
			eta_at_ob_old(i) = eta_at_ob_old(i) + eta_from_Holland_at_ob(i) ! jw
			eta_at_ob_new(i) = eta_at_ob_new(i) + eta_from_Holland_at_ob(i) ! jw
		end if
   end do   ! jw
   ! jw
end subroutine compute_boundary_eta