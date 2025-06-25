!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Prepare the GOM for the simulation
!! 
subroutine prepare_gom
	use mod_global_variables
	use mod_file_definition
	implicit none

   integer :: i
	character(len=40) :: format_string
	! End of local variables ==================================================!

	! set it to 0 for the beginning point =====================================!
	it = 0

	format_string = '(I4.4)'	! an integer of width 4 with zeros at the left. For example: 0001 ~ 9999
	
	! convert elapsed time to julian day and local time =======================!
	! if(restart_in == 0)then ! if it is cold start
	write(*,'(A35,F10.4)') 'Simulation start julian day = ', jday+1
	write(*,'(A35,I2.2,A1,I2.2,A1,I4,I3,A1,I2.2)') 				&
	&	'Simulation start local time = ', 							&
	&	start_month, '/', start_day, '/', start_year, start_hour, ':', start_minute 
	write(*,*) ! put one empty line

	julian_day = jday
	year       = start_year
	month      = start_month
	day        = start_day
	hour       = start_hour
	minute     = start_minute      
	! end if
   	! send back to the original value since these values will be calculated in main.f90 again
   	! current_jday_1900 = current_jday_1900 - dt/86400.0
	
	! time stepping ===========================================================!   
	! start time outeration iteration 
   if(restart_in /= 0) then
      ! it_hot_start = it_hot_start - 1
   else
      current_jday_1900 = jday_1900
   end if	
	
   write(pw_run_log,'(A35,I10)') 	'Time stepping begins at  [it] = ', 1 
   write(pw_run_log,'(A35,I10)') 	'Time stepping   ends at [ndt] = ', ndt
   write(pw_run_log,'(A35,F10.4)') 	'                         [dt] = ', dt
   write(pw_run_log,'(A35,F10.4)') 	'Total simulation time   [day] = ', ndt*dt/86400.0

   write(*,'(A35,I10)') 	'Time stepping begins at  [it] = ', 1
   write(*,'(A35,I10)')	 	'Time stepping   ends at [ndt] = ', ndt
   write(*,'(A35,F10.4)') 	'                         [dt] = ', dt
   write(*,'(A35,F10.4)') 	'Total simulation time   [day] = ', ndt*dt/86400.0
   write(*,*)
   write(*,*) "!==============================================================================!"
   
   elapsed_time      = it*dt
   elapsed_time_hr   = elapsed_time/3600.0
   ! current_jday_1900 = current_jday_1900 + dt/86400.0
   
   ! Note: here the [julian_day] is not actuall julian_day but [julian_day - 1].
   ! i.e., it starts from "0.0" day not from "1.0" day. 
   ! e.g., actual julian day at January 1st is 1.0 day, but in GOM it is 0.0 day (in elapsed julian day concept).
   ! but, now I will show the correct Julian day on the screen (i.e., julian_day + 1)
	write(*,'(A4,I10, A20,F12.5, A20,F12.5, A13,F12.5, A13,1X,I2.2,A1,I2.2,A1,I4,1X,I2.2,A1,I2.2)') 	&
	&     	' it=', it, 																						&
	&			', elapsed_time (hr)=', elapsed_time_hr,                        					&
	&      	', elapsed_time(day)=', elapsed_time/86400.0, 	                        		&
	&     	', julian day=', julian_day +1 ,                           							&
	&     	', local time=', month, '/', day, '/',year, hour, ':', minute 
	! =========================================================================!
	
		
	! open diagonostic files ==================================================!
	if(dia_advection == 1) then
		open(pw_dia_advection,file=id_dia_advection,form='formatted',status='replace')
		write(pw_dia_advection,'(A)') 'solve_nonlinear_advection.f90'
	end if
	
	if(dia_momentum == 1) then
		open(pw_dia_momentum,file=id_dia_momentum,form='formatted',status='replace')
		write(pw_dia_momentum,'(A)') 'solve_momentum_equation.f90'
	end if
	
	if(dia_freesurface == 1) then
		open(pw_dia_freesurface,file=id_dia_freesurface,form='formatted',status='replace')
		write(pw_dia_freesurface,'(A)') 'solve_free_surface_equation.f90'
	end if
	
	if(dia_eta_at_ob == 1) then
		open(pw_dia_eta_at_ob,file=id_dia_eta_at_ob,form='formatted',status='replace')
		write(pw_dia_eta_at_ob,'(A)') 'solve_free_surface_equation.f90'
	end if
	
	if(dia_bottom_friction == 1) then
		open(pw_dia_bottom_friction,file=id_dia_bottom_friction,form='formatted',status='replace')
		write(pw_dia_bottom_friction,'(A)') 'calculate_bottom_friction.f90'
	end if
	
	if(dia_face_velocity == 1) then
		! uv at each face
		open(pw_dia_face_velocity_uv,file=id_dia_face_velocity_uv,form='formatted',status='replace')
		write(pw_dia_face_velocity_uv,'(A)') 'calculate_horizontal_velocities.f90'
		
		! w at each element
		open(pw_dia_face_velocity_w,file=id_dia_face_velocity_w,form='formatted',status='replace')
		write(pw_dia_face_velocity_w,'(A)') 'calculate_vertical_velocities.f90'
	end if
	
	if(dia_node_velocity == 1) then
		open(pw_dia_node_velocity,file=id_dia_node_velocity,form='formatted',status='replace')
		write(pw_dia_node_velocity,'(A)') 'calculate_velocity_at_node.f90'
	end if
	! =========================================================================!
! 	write(*,*) 'jw001'
	! compute initial vertical layers =========================================!
   call update_vertical_layers_0	! jw, check this
! 	write(*,*) 'jw002'
	
	! initialize heat budget model ============================================!
	! jw, i_transport_model_flag has never been set.
!   if(i_transport_model_flag == 1) then
!      call initialize_transport_model
!   end if
! 	write(*,*) 'jw003'
	! Calculation of coriolis factor according to Coriolis options ============!
   call calculate_coriolis_factor
! 	write(*,*) 'jw004'
	
	! initialization for cold start alone =====================================!
   if(restart_in == 0) then ! if it is cold start
      call setup_cold_start
      call calculate_velocity_at_node
   else if(restart_in == 1) then
   	! velocities at node (e.g., u_node, v_node, and w_node) will be read from restart.inp
   	call setup_hot_start
   end if
!  write(*,*) 'jw005'
	! compute initial vertical layers =========================================!
   ! call update_vertical_layers_0	! jw, check this

	! compute velocities at nodes =============================================!
   ! call calculate_velocity_at_node
	
	! calculate initial density ===============================================!
	! This is required to solve baroclinic gradient equation in the momentum equation.
	! If one of "is_salt" or "is_temp" is activated but "baroclinic_flag" is off, only transport equation will be solved without density flow.
	! If one of "is_salt" or "is_temp" is activated and "baroclinic_flag" is on, transport equation will be solved and also the density flow will be calulated.
	if(baroclinic_flag == 1) then
		if(ana_density == 0) then
			call calculate_density_full
		else if(ana_density == 1) then
			call calculate_density_linear
		else if(ana_density == 2) then
			call calculate_analytical_density
		end if
	end if
! 	write(*,*) 'jw006'
	
	! set initial elevation field (hot start) for relaxation in sponge layer
	! jw, i_sponge_layer_flag has never been set.
   if(i_sponge_layer_flag /= 0) then
      do i = 1, maxele
         etaic(i) = eta_cell(i)
      end do
   end if
! 	write(*,*) 'jw007'
	
	! calculate vertically averaged u,v,salt,temp,rho for initial value write =!
	call calculate_vertically_averaged_values

	! Before starting simulation, write grid checking files: check_grid2D.dat, check_grid2DO.dat, check_grid_3D.dat
	!$omp parallel sections
	!$omp section
	if(check_grid_2D == 1 .or. check_grid_2DO == 1 .or. check_grid_3D == 1) then
		call write_grid_checking_files
	end if
! 	write(*,*) 'jw008'
	
	! === open output files and write initial conditions ======================!
   ! prepare timeseries out files: ~_tser.dat
   !$omp section
   if(tser_station_num > 0) then
	   call open_tser_out_files
	end if
! 	write(*,*) 'jw009' 
	
	! open 2D output files: for Tecplot & vtk
	!$omp section
	if(IS2D_switch == 1) then
		call open_2D_out_files
		
		if(IS2D_flood_map == 1) then
			! initialize flood_id
			! If the node is initially set to land, set "flood_id" to "1"
			do i=1,maxnod
				if(h_node(i) <= 0.0_dp) then ! node is dry, if h_node <= 0
					flood_id(i) = 1
				end if
			end do
		end if
	end if
! 	write(*,*) 'jw0010'
	
	! open 3D output files: for Tecplot & vtk
	!$omp section
	if(IS3D_surf_switch == 1 .or. IS3D_full_switch == 1) then
		call open_3D_out_files
	end if
! 	write(*,*) 'jw0011'
	
	! open 2D dump output files: for dump
	!$omp section
	if(IS2D_dump_switch == 1) then
		call open_dump2D
	end if
! 	write(*,*) 'jw0012'
	
	! open 3D dump output files: for dump
	!$omp section
	if(IS3D_dump_switch == 1) then
		call open_dump3D
	end if
! 	write(*,*) 'jw0013'

	!$omp section
	if(netcdf_switch == 1) then
		! initialize netcdf_File_num to "0" for the later use
		netcdf_File_num = 0
		
		! total number of time series (data) in each netcdf dump file, this will be identical for all netcdf files,
		! and that is why I am calculating one time.
		netcdf_maxtime = INT(netcdf_File_Freq/netcdf_frequency)
		
		! write one-time grid information
		call write_netcdf_grid
	end if
	!$omp end parallel sections
end subroutine prepare_gom