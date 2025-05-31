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
	! jw

	! jw
	it = 0

	format_string = '(I4.4)'	! jw
	
	! jw
	! jw
	write(*,'(A35,F10.4)') 'Simulation start julian day = ', jday+1
	write(*,'(A35,I2.2,A1,I2.2,A1,I4,I3,A1,I2.2)') 				&
	&	'Simulation start local time = ', 							&
	&	start_month, '/', start_day, '/', start_year, start_hour, ':', start_minute 
	write(*,*) ! jw

	julian_day = jday
	year       = start_year
	month      = start_month
	day        = start_day
	hour       = start_hour
	minute     = start_minute      
	! jw
   	! jw
   	! jw
	
	! jw
	! jw
   if(restart_in /= 0) then
      ! jw
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
   ! jw
   
   ! jw
   ! jw
   ! jw
   ! jw
	write(*,'(A4,I10, A20,F12.5, A20,F12.5, A13,F12.5, A13,1X,I2.2,A1,I2.2,A1,I4,1X,I2.2,A1,I2.2)') 	&
	&     	' it=', it, 																						&
	&			', elapsed_time (hr)=', elapsed_time_hr,                        					&
	&      	', elapsed_time(day)=', elapsed_time/86400.0, 	                        		&
	&     	', julian day=', julian_day +1 ,                           							&
	&     	', local time=', month, '/', day, '/',year, hour, ':', minute 
	! jw
	
		
	! jw
	if(dia_advection == 1) then
		open(pw_dia_advection,file=id_dia_advection,form='formatted',status='replace')
		write(pw_dia_advection,'(A)') 'solve_nonlinear_advection.f90'
	end if
	
	if(dia_momentum == 1) then
		open(pw_dia_momentum,file=id_dia_momentum,form='formatted',status='replace')
		write(pw_dia_advection,'(A)') 'solve_momentum_equation.f90'
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
		! jw
		open(pw_dia_face_velocity_uv,file=id_dia_face_velocity_uv,form='formatted',status='replace')
		write(pw_dia_face_velocity_uv,'(A)') 'calculate_horizontal_velocities.f90'
		
		! jw
		open(pw_dia_face_velocity_w,file=id_dia_face_velocity_w,form='formatted',status='replace')
		write(pw_dia_face_velocity_w,'(A)') 'calculate_vertical_velocities.f90'
	end if
	
	if(dia_node_velocity == 1) then
		open(pw_dia_node_velocity,file=id_dia_node_velocity,form='formatted',status='replace')
		write(pw_dia_node_velocity,'(A)') 'calculate_velocity_at_node.f90'
	end if
	! jw
! jw
	! jw
   call update_vertical_layers_0	! jw
! jw
	
	! jw
	! jw
! jw
! jw
! jw
! jw
	! jw
   call calculate_coriolis_factor
! jw
	
	! jw
   if(restart_in == 0) then ! jw
      call setup_cold_start
      call calculate_velocity_at_node
   else if(restart_in == 1) then
   	! jw
   	call setup_hot_start
   end if
! jw
	! jw
   ! jw

	! jw
   ! jw
	
	! jw
	! jw
	! jw
	! jw
	if(baroclinic_flag == 1) then
		if(ana_density == 0) then
			call calculate_density_full
		else if(ana_density == 1) then
			call calculate_density_linear
		else if(ana_density == 2) then
			call calculate_analytical_density
		end if
	end if
! jw
	
	! jw
	! jw
   if(i_sponge_layer_flag /= 0) then
      do i = 1, maxele
         etaic(i) = eta_cell(i)
      end do
   end if
! jw
	
	! jw
	call calculate_vertically_averaged_values

	! jw
	!$omp parallel sections
	!$omp section
	if(check_grid_2D == 1 .or. check_grid_2DO == 1 .or. check_grid_3D == 1) then
		call write_grid_checking_files
	end if
! jw
	
	! jw
   ! jw
   !$omp section
   if(tser_station_num > 0) then
	   call open_tser_out_files
	end if
! jw
	
	! jw
	!$omp section
	if(IS2D_switch == 1) then
		call open_2D_out_files
		
		if(IS2D_flood_map == 1) then
			! jw
			! jw
			do i=1,maxnod
				if(h_node(i) <= 0.0_dp) then ! jw
					flood_id(i) = 1
				end if
			end do
		end if
	end if
! jw
	
	! jw
	!$omp section
	if(IS3D_surf_switch == 1 .or. IS3D_full_switch == 1) then
		call open_3D_out_files
	end if
! jw
	
	! jw
	!$omp section
	if(IS2D_dump_switch == 1) then
		call open_dump2D
	end if
! jw
	
	! jw
	!$omp section
	if(IS3D_dump_switch == 1) then
		call open_dump3D
	end if
! jw
	!$omp end parallel sections
end subroutine prepare_gom