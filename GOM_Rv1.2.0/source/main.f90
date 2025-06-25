!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
program main
! 	use omp_lib
   use mod_global_variables
   use mod_file_definition
	implicit none
	
	integer :: i, j, k
   ! End of local variables ==================================================!
! 	call omp_set_nested(.true.)
	
	! write a welcome messange on run.log and the screen ======================!
	call welcome_message
! 	write(*,*) 'jw01'
	! Scan & Go step ==========================================================!
	! Find the maximum array size for allocatable variables & allocate global variables
	! Most allocatable global variables should be allocated before "Time marching start!"
	! "scan_gom" includes "set_geometry_1.f90", and some variables are allocated in "set_geometry_1.f90"
	! Most of other variables will be allocated in "allocate_variables.f90"
	call scan_gom 
! 	write(*,*) 'jw02'	
	
	! find "maxface", allocate some variables, and set geometry
	call set_geometry_1
! 	write(*,*) 'jw03'	
	
	! allocate rest global variables
	call allocate_variables
! 	write(*,*) 'jw04'	
	
	! This includes "call set_geometry_2"
	call read_main_inp
! 	write(*,*) 'jw05'	
	! If space and time interpolation is activated, 
	! create a regular grid or find searching neighbour nodes at each node before start simulation for speeding up IDW process.
	! The second method (find_neighbour_nodes) is more stable and faster.
	if(hurricane_flag == 1 .and. hurricane_interp_method == 2) then
		! call create_regular_grid
		call find_neighbour_nodes ! this is more stable
	end if
! 	write(*,*) 'jw06'	
	! Now you are ready to go =================================================!
	
	! prepare the GOM for the simulation ======================================!
	call prepare_gom
! 	write(*,*) 'jw07'	

!=============================================================================!
	! Time marching start!
!=============================================================================!
	do it = 1, ndt
		elapsed_time      = it*dt	! [sec]
      elapsed_time_hr   = elapsed_time/3600.0_dp
      julian_day        = julian_day + dt/86400.0_dp 				! elapsed julian day from the given start year; see C3 in main.inp
      ! or
      ! julian_day			= jday + (it*dt)/86400.0
      current_jday_1900 = current_jday_1900  + dt/86400.0_dp	! julian day from year of 1900. I think this value can be deleted (jw).

		! convert elapsed time into local time =================================!
		! write(*,*) year, month,day,hour,minute,julian_day
      call julian_to_localtime(year, month, day, hour, minute, julian_day)
      
		! screen output of elapaed time of simulation
	   ! Note: here the [julian_day] is not actuall julian_day but [julian_day - 1].
	   ! i.e., it starts from "0.0" day not from "1.0" day. 
	   ! e.g., actual julian day at January 1st is 1.0 day, but in GOM it is 0.0 day (in elapsed julian day concept); 
	   ! but, now I will show the correct Julian day on the screen (i.e., julian_day + 1)
		if(ishow == 1) then
	      if(mod(it,ishow_frequency) == 0) then
	         write(*,'(A4,I10, A20,F12.5, A20,F12.5, A13,F12.5, A13,1X,I2.2,A1,I2.2,A1,I4,1X,I2.2,A1,I2.2)') 	&
	         &     	' it=', it, 																						&
	         &			', elapsed_time (hr)=', elapsed_time_hr,                        					&
	         &      	', elapsed_time(day)=', elapsed_time/86400.0_dp,                         		&
	         &     	', julian day=', julian_day +1 ,                           							&
	         &     	', local time=', month, '/', day, '/',year, hour, ':', minute 
	      end if
	   end if
	   ! write(*,*) julian_day, month, day, hour, minute
		
		! before start, define spinup function for boundary elevation forcing, 
		! wind and pressure forcing, and tidal potential forcing ===============!
		call setup_spinup		
! 		write(*,*) 'jw1'

		!=======================================================================!
		! First,
		! Boundary conditions at the bottom and surface should be calculated 
		! before solving the momentum equation.
		!=======================================================================!
		! Update bottom boundary condition =====================================!
		! calculate bottom friction --------------------------------------------!
		! This will calculate Cd (drag coefficient) and Gamma_B (bottom stress coefficient, Cd * sqrt(u^2 * v^2))
      call calculate_bottom_friction ! omp done
! 		write(*,*) 'jw2'
! 		stop 'jw2'
      
      ! Update surface boundary condition ====================================!
		! include wind and air pressure ----------------------------------------!
		if(wind_flag == 1 .or. airp_flag == 1) then
			if(num_windp_ser > 0) then
      		call get_windp_at_face ! omp done
      	end if
      end if
! 		write(*,*) 'jw3'
      
		! include storm surge model --------------------------------------------!
		! Let's skip Holland storm surge now...
! 	   if(holland_flag == 1) then
! 	      call holland_storm_surge
! 	   end if
! 		write(*,*) 'jw4'		

		! include hurricane data (for storm surge simulation) ------------------!
		if(hurricane_flag == 1) then
		   if(julian_day >= hurricane_start_jday .and. julian_day <= hurricane_end_jday) then
	         if(hurricane_interp_method == 1) then
	         	! linear interpolation version
	            call read_hurricane_ser_1 ! omp done
	         elseif(hurricane_interp_method == 2) then
	         	! space & time interpolation version: Inverse Distance Weighting (IDW) Interpolation
	            call read_hurricane_ser_2 ! omp done
	         end if
	      end if
      end if
! 		write(*,*) 'jw5'
				      
		! Calculate wind stress from given wind speed at each face -------------!
		if(wind_flag == 0) then	! no wind applied to the model
			wind_stress_normal = 0.0_dp
			wind_stress_tangnt = 0.0_dp
		end if
		
		if(wind_flag == 1 .or. hurricane_flag == 1) then
			! This will calculate [wind_stress_normal] & [wind_stress_tangnt] at each face
      	call calculate_wind_stress ! omp done
		end if
		
		! Analytical wind stress will be used.
		if(wind_flag == 2) then
			call ana_windp ! omp done
		end if
! 		write(*,*) 'jw6'		

		!=======================================================================!
		! Second,
		! Boundary velocities (from open boundary) should be included 
		! before solving the momentum equation.
		!=======================================================================!
		! compute temporary velocitiy at boundary
		! True river or tidal river boundary condition
     	call compute_boundary_velocity
! 		write(*,*) 'jw7'		
		
		!=======================================================================!
		! Third,
		! Now it is ready to solve the momentum equation.
		!=======================================================================!
		! These are the main part of the GOM ===================================!
		! (0) Before start the main part, let's calculate water surface elevation at the boundary ghost cells from given data:
		call compute_boundary_eta
! 		write(*,*) 'jw8'		
		
		! (1) solve nonlinear advection
		! This will calculate face_normal & face_tangential velocities at (n+1, which is the current) time step, from the nonlinear advection term
		! solve_nonlinear_advection
		if(advection_flag == 0) then
			! To include tendency term automatically even though turning off nonlinear advection term,
			! we include the tendency term in the nonlinear advection equation.
			! That is why we are doing this.
			do j=1,maxface
				do k=1,maxlayer
 					un_ELM(k,j) = un_face(k,j)
 					vn_ELM(k,j) = vn_face(k,j)
				end do
			end do
		else
			! When advection is on, it will include tendency term + nonlinear advection term.
			call solve_nonlinear_advection ! omp done
		end if
! 		write(*,*) 'jw9'		
		
		! write diagnostic file
		if(dia_advection == 1) then
			write(pw_dia_advection,*) 'it = ', it, ', elapsed_time = ', elapsed_time
			write(pw_dia_advection,*) 'un_ELM(1,j), vn_ELM(1,j)'
			do j=1,maxface
				write(pw_dia_advection,'(A3, I5, 2f10.4)') 'j=', j, un_ELM(1,j), vn_ELM(1,j)
			end do
		end if
		
		! (2) calculate matrix and vectors in the momentum equation
      call solve_momentum_equation ! omp done
! 		write(*,*) 'jw10'
! 		stop 'jw10'

		!=======================================================================!
		! Fourth,
		! Now solve the free surface equation & calculate velocities
		!=======================================================================!
		! (3) solve wave equation using the calculated matrix: free surface wave equation
		call solve_free_surface_equation ! jw's jcg version, omp done
! 		write(*,*) 'jw11'
! 		stop 'jw11'

		! Solve momentum Eq(44) & Continuity Eq(47) to calculate velocities at face:
		! 		u & v: horizontal face center
		! 		w: vertical face center
		call solve_velocities ! omp done
! 		write(*,*) 'jw12'				

		call update_vertical_layers_1	! omp done
! 		write(*,*) 'jw13'		

		call calculate_velocity_at_node ! omp done
! 		write(*,*) 'jw14'

		! Solve transport equations ============================================!
		if(transport_flag == 1) then
			! if heat transport is activated, calculate heat exchange first, then solve transport equation:
			if(is_temp == 1) then
				if(heat_option == 1) then
					call heat_term_by_term
				else if(heat_option == 2) then ! not yet included
					call heat_equilibrium_temperature
				end if
			end if
			
			if(transport_solver == 1) then
				! use TVD approach
				if(maxval(is_tran) > 0) then ! double check if transport is activated.
					! call solve_transport_equation_v18
					call solve_transport_equation_v19
				else
					write(pw_run_log,*) 'Error: solve_transport_equation.f90 - transport material is not activated.'
					write(*,*) 'Error: solve_transport_equation.f90 - transport material is not activated.'
					stop
				end if
			else if(transport_solver == 2) then
				! use ELM approach
				if(maxval(is_tran) > 0) then ! double check if transport is activated.
					! if(mod(it,trans_ext_iter) == 0) then
						call solve_transport_equation_ELM_v8
					! end if
				else
					write(pw_run_log,*) 'Error: solve_transport_equation.f90 - transport material is not activated.'
					write(*,*) 'Error: solve_transport_equation.f90 - transport material is not activated.'
					stop
				end if
			end if
		end if
! 		write(*,*) 'jw15'
		! write(*,*) un_face_new(1,61), un_face_new(1,64), un_face_new(1,66)
		! write(*,*) un_face_new(1,68), un_face_new(1,71), un_face_new(1,73)
		
		
		! After solving transport equation, we have to update density at face, node, and cell.
		! If baroclinic_flag is off (i.e., 0), baroclinic velocity will not be included in the momentum equation.
		! That means that eventhough salinity trasport is activated, baroclinic velocity either can or cannot be included by setting the "baroclinic_flag".
		if(baroclinic_flag == 1) then
			if(ana_density == 0) then
				! write(*,*) 'jw15_1'
				call calculate_density_full ! omp done
				! write(*,*) 'jw15_2'
			else if(ana_density == 1) then
				call calculate_density_linear
			else if(ana_density == 2) then ! omp done
				call calculate_analytical_density ! omp done
			end if
		end if
! 		write(*,*) 'jw16'
		
		! Update old variables to the new variables for the next time step =====!
		call update_variables ! omp done
! 		write(*,*) 'jw17'

		! write output files ===================================================!		
		call write_output_files ! omp done
		if(IS2D_switch == 1) then
			if(IS2D_flood_map > 0) then
				call update_flood_map ! omp done
			end if
		end if
! 		write(*,*) 'jw18'		
		
		! write restart.out ====================================================!
		if(restart_out == -1) then
			! write restart file once at the end of simulation
			if(it==ndt) then
				call write_restart
			end if
		else if(restart_out == 1) then
			! write restart file every restart_freq times
			if(mod(it,restart_freq) == 0) then
				call write_restart
			end if
		end if
! 		write(*,*) 'jw19'		
		
		! check if termination required based on min/max values ================!
		if(mod(it,terminate_check_freq) == 0) then
			call model_termination_check
		end if
		
		! stop 'I am here'
   end do ! it=1,ndt   
!=== End of time marching ====================================================!
	
	! write flood map after the end of the simulation =========================!
	if(IS2D_switch > 0) then
		if(IS2D_flood_map > 0) then
			!$omp parallel sections
			!$omp section
			if(IS2D_format == 1 .or. IS2D_format == 3) then
				call write_flood_map_tec
			end if
			!$omp section
			if(IS2D_format == 2 .or. IS2D_format == 3) then
				call write_flood_map_vtk
			end if
			!$omp end parallel sections
		end if
	end if

	! === call closing message ================================================!
	call closing_message
	stop

end program main
