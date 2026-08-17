!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine open_2D_out_files
	use mod_global_variables
	use mod_file_definition	
	implicit none
	
	character(len=40) :: format_string
	character(len= 4) :: File_num_buff
	! jw
	
	format_string = '(I4.4)'	! jw
	
	! jw
	if(IS2D_time == 1) then	! jw
		IS2D_time_conv = 86400.0
	else if(IS2D_time == 2) then ! jw
		IS2D_time_conv = 1440.0
	else if(IS2D_time == 3) then ! jw
		IS2D_time_conv = 24.0
	else if(IS2D_time == 4) then ! jw
		IS2D_time_conv = 1.0
	end if
	
	! jw
	IS2D_File_num = 1 ! jw
	zone_num_2D = 1	! jw
	IS2D_vtk_num = 0 	! jw
		
	!=== 2D contour plot files ================================================!
	! jw
	if(IS2D_flood_map == 1) then
		! jw
		if(IS2D_format == 1 .or. IS2D_format == 3) then
			if(IS2D_binary == 0) then
				open(pw_flood_map_tec, file = trim(id_flood_map_tec), form = 'formatted', status = 'replace')
			else if(IS2D_binary == 1) then 
				! jw
			end if
		end if
		
		! jw
		if(IS2D_format == 2 .or. IS2D_format == 3) then
			if(IS2D_binary == 0) then
				open(pw_flood_map_vtk, file = trim(id_flood_map_vtk), form = 'formatted', status = 'replace')
			else if(IS2D_binary == 1) then 
				! jw
			end if
		end if
	end if
	
	! jw
	! jw
	if(IS2D_format == 1 .or. IS2D_format == 3) then
		write(File_num_buff,format_string) IS2D_File_num		
		IS2D_File_name = trim(id_tec2D)//trim(File_num_buff)//'.dat'
		
		if(IS2D_binary == 0) then
			open(pw_tec2D, file = trim(IS2D_File_name), form = 'formatted', status = 'replace')	! jw
			call write_tec2D_head
		else if(IS2D_binary == 1) then 
			! jw
			! jw
		end if
	end if	
	
	! jw
	! jw
! jw
		! jw
		! jw
		
		! jw
! jw
	
end subroutine open_2D_out_files
