!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine open_2D_out_files
	use mod_global_variables
	use mod_file_definition	
	implicit none
	
	character(len=40) :: format_string
	character(len= 4) :: File_num_buff
	! End of local variables --------------------------------------------------!
	
	format_string = '(I4.4)'	! an integer of width 4 with zeros at the left. For example: 0001 ~ 9999
	
	! define time conversion factor
	if(IS2D_time == 1) then	! second
		IS2D_time_conv = 86400.0
	else if(IS2D_time == 2) then ! minute
		IS2D_time_conv = 1440.0
	else if(IS2D_time == 3) then ! hour
		IS2D_time_conv = 24.0
	else if(IS2D_time == 4) then ! day
		IS2D_time_conv = 1.0
	end if
	
	! set the first file number and zone number
	IS2D_File_num = 1 ! for tecplot
	zone_num_2D = 1	! for tecplot
	IS2D_vtk_num = 0 	! for vtk file format
		
	!=== 2D contour plot files ================================================!
	! (1) flood map plot ------------------------------------------------------!
	if(IS2D_flood_map == 1) then
		! tecplot file:
		if(IS2D_format == 1 .or. IS2D_format == 3) then
			if(IS2D_binary == 0) then
				open(pw_flood_map_tec, file = trim(id_flood_map_tec), form = 'formatted', status = 'replace')
			else if(IS2D_binary == 1) then 
				! not yet included				
			end if
		end if
		
		! vtk file:
		if(IS2D_format == 2 .or. IS2D_format == 3) then
			if(IS2D_binary == 0) then
				open(pw_flood_map_vtk, file = trim(id_flood_map_vtk), form = 'formatted', status = 'replace')
			else if(IS2D_binary == 1) then 
				! not yet included
			end if
		end if
	end if
	
	! (2) main 2D plots -------------------------------------------------------!
	! tecplot file:
	if(IS2D_format == 1 .or. IS2D_format == 3) then
		write(File_num_buff,format_string) IS2D_File_num		
		IS2D_File_name = trim(id_tec2D)//trim(File_num_buff)//'.dat'
		
		if(IS2D_binary == 0) then
			open(pw_tec2D, file = trim(IS2D_File_name), form = 'formatted', status = 'replace')	! tec2D_0000.dat
			call write_tec2D_head
		else if(IS2D_binary == 1) then 
			! not yet included
			! open(pw_tec2D_binary, file = id_tec2D_binary, form = 'formatted', status = 'replace')	! 
		end if
	end if	
	
	! vtk file:
	! don't require this part for vtk file format
! 	if(IS2D_format == 2 .or. IS2D_format == 3) then
		! write(File_num_buff,format_string) IS2D_File_num		
		! IS2D_File_name = trim(id_vtk2D)//trim(File_num_buff)//'.vtk'
		
		! open(pw_vtk2D, file = trim(IS2D_File_name), form = 'formatted', status = 'replace')	! vtk2D_0000.vtk
! 	end if
	
end subroutine open_2D_out_files
