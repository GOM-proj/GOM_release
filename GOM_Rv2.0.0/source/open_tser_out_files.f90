!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! write the time series data
!! water elevation, total waer depth, velocities, salinity, temperature, and air pressure
!! Note: Header lines are different for Tecplot (***.dat) and Paraview (*.txt)
!! 
subroutine open_tser_out_files
   use mod_global_variables
   use mod_file_definition
   implicit none
   
   integer :: pw_tser
   character(len=200) :: id_tser
   integer :: tser_list
   ! jw

	! jw
	! jw
	if(tser_time == 1) then	! jw
		tser_time_conv = 86400.0
	else if(tser_time == 2) then ! jw
		tser_time_conv = 1440.0
	else if(tser_time == 3) then ! jw
		tser_time_conv = 24.0
	else if(tser_time == 4) then ! jw
		tser_time_conv = 1.0
	end if

	! jw
	! jw
 	if(tser_eta == 1) then
 		tser_list = 1
 		if(tser_format == 1 .or. tser_format == 3) then
	 		pw_tser = pw_eta_tser_from_msl_tec
	 		id_tser = id_eta_tser_from_msl_tec
	 		call open_tser_outfiles_tec(pw_tser,id_tser,tser_list)
	 	end if
 		if(tser_format == 2 .or. tser_format == 3) then
	 		pw_tser = pw_eta_tser_from_msl_vtk
	 		id_tser = id_eta_tser_from_msl_vtk
	 		call open_tser_outfiles_vtk(pw_tser,id_tser,tser_list)
	 	end if	 		
 	end if

 	if(tser_H == 1) then
 		tser_list = 2
 		if(tser_format == 1 .or. tser_format == 3) then
 			pw_tser = pw_H_tser_tec
 			id_tser = id_H_tser_tec
 			call open_tser_outfiles_tec(pw_tser,id_tser,tser_list)
 		end if
 		if(tser_format == 2 .or. tser_format == 3) then
 			pw_tser = pw_H_tser_vtk
 			id_tser = id_H_tser_vtk
 			call open_tser_outfiles_vtk(pw_tser,id_tser,tser_list)
 		end if 		
 	end if

	if(tser_u == 1) then
		tser_list = 3
		if(tser_format == 1 .or. tser_format == 3) then
			pw_tser = pw_u_tser_tec
			id_tser = id_u_tser_tec
			call open_tser_outfiles_tec(pw_tser,id_tser,tser_list)
		end if
		if(tser_format == 2 .or. tser_format == 3) then
			pw_tser = pw_u_tser_vtk
			id_tser = id_u_tser_vtk
			call open_tser_outfiles_vtk(pw_tser,id_tser,tser_list)
		end if		
	end if	

	if(tser_v == 1) then
		tser_list = 4
		if(tser_format == 1 .or. tser_format == 3) then
			pw_tser = pw_v_tser_tec
			id_tser = id_v_tser_tec
			call open_tser_outfiles_tec(pw_tser,id_tser,tser_list)
		end if
		if(tser_format == 2 .or. tser_format == 3) then
			pw_tser = pw_v_tser_vtk
			id_tser = id_v_tser_vtk
			call open_tser_outfiles_vtk(pw_tser,id_tser,tser_list)
		end if		
	end if	

	if(tser_salt == 1) then
		tser_list = 5
		if(tser_format == 1 .or. tser_format == 3) then
			pw_tser = pw_salt_tser_tec
			id_tser = id_salt_tser_tec
			call open_tser_outfiles_tec(pw_tser,id_tser,tser_list)
		end if
		if(tser_format == 2 .or. tser_format == 3) then
			pw_tser = pw_salt_tser_vtk
			id_tser = id_salt_tser_vtk
			call open_tser_outfiles_vtk(pw_tser,id_tser,tser_list)
		end if		
	end if	

	if(tser_temp == 1) then ! jw
		tser_list = 6
		if(tser_format == 1 .or. tser_format == 3) then
			pw_tser = pw_temp_tser_tec
			id_tser = id_temp_tser_tec
			call open_tser_outfiles_tec(pw_tser,id_tser,tser_list)
		end if
		if(tser_format == 2 .or. tser_format == 3) then
			pw_tser = pw_temp_tser_vtk
			id_tser = id_temp_tser_vtk
			call open_tser_outfiles_vtk(pw_tser,id_tser,tser_list)
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
! jw
! jw
! jw
! jw
! jw
	
end subroutine open_tser_out_files
