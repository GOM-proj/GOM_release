!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
program main
 	use omp_lib
   use mod_global_variables
   use mod_file_definition
	implicit none
	
	integer :: i, j, k
   ! jw
 	call omp_set_nested(.true.)
	
	! jw
	call welcome_message
! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	call scan_gom 
! jw
	
	! jw
	call set_geometry_1
! jw
	
	! jw
	call allocate_variables
! jw
	
	! jw
	call read_main_inp
! jw
	! jw
	! jw
	! jw
	if(hurricane_flag == 1 .and. hurricane_interp_method == 2) then
		! jw
		call find_neighbour_nodes ! jw
	end if
! jw
	! jw
	
	! jw
	call prepare_gom
! jw

!=============================================================================!
	! jw
!=============================================================================!
	do it = 1, ndt
		elapsed_time      = it*dt	! jw
      elapsed_time_hr   = elapsed_time/3600.0_dp
      julian_day        = julian_day + dt/86400.0_dp 				! jw
      ! jw
      ! jw
      current_jday_1900 = current_jday_1900  + dt/86400.0_dp	! jw

		! jw
		! jw
      call julian_to_localtime(year, month, day, hour, minute, julian_day)
      
		! jw
	   ! jw
	   ! jw
	   ! jw
	   ! jw
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
	   ! jw
		
		! jw
		! jw
		call setup_spinup		
! jw

		!=======================================================================!
		! jw
		! jw
		! jw
		!=======================================================================!
		! jw
		! jw
		! jw
      call calculate_bottom_friction ! jw
! jw
      
      ! jw
		! jw
		if(wind_flag == 1 .or. airp_flag == 1) then
			if(num_windp_ser > 0) then
      		call get_windp_at_face ! jw
      	end if
      end if
! jw
      
		! jw
		! jw
! jw
! jw
! jw
! jw

		! jw
		if(hurricane_flag == 1) then
		   if(julian_day >= hurricane_start_jday .and. julian_day <= hurricane_end_jday) then
	         if(hurricane_interp_method == 1) then
	         	! jw
	            call read_hurricane_ser_1 ! jw
	         elseif(hurricane_interp_method == 2) then
	         	! jw
	            call read_hurricane_ser_2 ! jw
	         end if
	      end if
      end if
! jw
				      
		! jw
		if(wind_flag == 0) then	! jw
			wind_stress_normal = 0.0_dp
			wind_stress_tangnt = 0.0_dp
		end if
		
		if(wind_flag == 1 .or. hurricane_flag == 1) then
			! jw
      	call calculate_wind_stress ! jw
		end if
		
		! jw
		if(wind_flag == 2) then
			call ana_windp ! jw
		end if
! jw

		!=======================================================================!
		! jw
		! jw
		! jw
		!=======================================================================!
		! jw
		! jw
     	call compute_boundary_velocity
! jw
		
		!=======================================================================!
		! jw
		! jw
		!=======================================================================!
		! jw
		! jw
		call compute_boundary_eta
! jw
		
		! jw
		! jw
		! jw
		if(advection_flag == 0) then
			! jw
			! jw
			! jw
			do j=1,maxface
				do k=1,maxlayer
					un_ELM(k,j) = un_face(k,j)
					vn_ELM(k,j) = vn_face(k,j)
				end do
			end do
		else
			! jw
			call solve_nonlinear_advection ! jw
		end if
! jw
		
		! jw
		if(dia_advection == 1) then
			write(pw_dia_advection,*) 'it = ', it, ', elapsed_time = ', elapsed_time
			write(pw_dia_advection,*) 'un_ELM(1,j), vn_ELM(1,j)'
			do j=1,maxface
				write(pw_dia_advection,'(A3, I5, 2f10.4)') 'j=', j, un_ELM(1,j), vn_ELM(1,j)
			end do
		end if
		
		! jw
      call solve_momentum_equation ! jw
! jw

		!=======================================================================!
		! jw
		! jw
		!=======================================================================!
		! jw
		call solve_free_surface_equation ! jw
! jw

		! jw
		! jw
		! jw
		call solve_velocities ! jw
! jw

		call update_vertical_layers_1	! jw
! jw

		call calculate_velocity_at_node ! jw
! jw

		! jw
		if(transport_flag == 1) then
			! jw
			if(is_temp == 1) then
				if(heat_option == 1) then
					call heat_term_by_term
				else if(heat_option == 2) then ! jw
					call heat_equilibrium_temperature
				end if
			end if
			
			if(transport_solver == 1) then
				! jw
				if(maxval(is_tran) > 0) then ! jw
					! jw
					call solve_transport_equation_v19
				else
					write(pw_run_log,*) 'Error: solve_transport_equation.f90 - transport material is not activated.'
					write(*,*) 'Error: solve_transport_equation.f90 - transport material is not activated.'
					stop
				end if
			else if(transport_solver == 2) then
				! jw
				if(maxval(is_tran) > 0) then ! jw
					! jw
						call solve_transport_equation_ELM_v8
					! jw
				else
					write(pw_run_log,*) 'Error: solve_transport_equation.f90 - transport material is not activated.'
					write(*,*) 'Error: solve_transport_equation.f90 - transport material is not activated.'
					stop
				end if
			end if
		end if
! jw
		! jw
		! jw
		
		
		! jw
		! jw
		! jw
		if(baroclinic_flag == 1) then
			if(ana_density == 0) then
				! jw
				call calculate_density_full ! jw
				! jw
			else if(ana_density == 1) then
				call calculate_density_linear
			else if(ana_density == 2) then ! jw
				call calculate_analytical_density ! jw
			end if
		end if
! jw
		
		! jw
		call update_variables ! jw
! jw

		! jw
		call write_output_files ! jw
		if(IS2D_switch == 1) then
			if(IS2D_flood_map > 0) then
				call update_flood_map ! jw
			end if
		end if
! jw
		
		! jw
		if(restart_out == -1) then
			! jw
			if(it==ndt) then
				call write_restart
			end if
		else if(restart_out == 1) then
			! jw
			if(mod(it,restart_freq) == 0) then
				call write_restart
			end if
		end if
! jw
		
		! jw
		if(mod(it,terminate_check_freq) == 0) then
			call model_termination_check
		end if
		
		! jw
   end do ! jw
!=== End of time marching ====================================================!
	
	! jw
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

	! jw
	call closing_message
	stop

end program main
