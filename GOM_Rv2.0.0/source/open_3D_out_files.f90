!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine open_3D_out_files
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	character(len=40) :: format_string	
	character(len= 4) :: File_num_surf_buff, File_num_full_buff
	! jw

	format_string = '(I4.4)'	! jw

	! jw
	if(IS3D_time == 1) then	! jw
		IS3D_time_conv = 86400.0
	else if(IS3D_time == 2) then ! jw
		IS3D_time_conv = 1440.0
	else if(IS3D_time == 3) then ! jw
		IS3D_time_conv = 24.0
	else if(IS3D_time == 4) then ! jw
		IS3D_time_conv = 1.0
	end if
	
	! jw
	IS3D_File_num_surf = 1
	IS3D_File_num_full = 1
	zone_num_3D_surf = 1
	zone_num_3D_full = 1		
	IS3D_vtk_num = 0

	! jw
	! jw
	! jw
	!=== 3D surface contour plot files ========================================!
	if(IS3D_surf_switch == 1) then
		! jw
		if(IS3D_format == 1 .or. IS3D_format == 3) then
			if(IS3D_binary == 0) then
				write(File_num_surf_buff,format_string) IS3D_File_num_surf
				IS3D_File_name_surf = trim(id_tec3D_surf)//trim(File_num_surf_buff)//'.dat'
				
				open(pw_tec3D_surf, file = trim(IS3D_File_name_surf), form = 'formatted', status = 'replace')	! jw
				call write_tec3D_surf_head
			else if(IS3D_binary == 1) then
				! jw
			end if
		end if
		
		! jw
		! jw
		! jw
		! jw
	end if

	!=== 3D full contour plot files ===========================================!
	if(IS3D_full_switch == 1) then
		! jw
		if(IS3D_format == 1 .or. IS3D_format == 3) then
			if(IS3D_binary == 0) then
				write(File_num_full_buff,format_string) IS3D_File_num_full
				IS3D_File_name_full = trim(id_tec3D_full)//trim(File_num_full_buff)//'.dat'

				open(pw_tec3D_full, file = trim(IS3D_File_name_full), form = 'formatted', status = 'replace')	! jw
				if(IS3D_grid_format == 1) then
					call write_tec3D_full_head_zto_sigma
				else if(IS3D_grid_format == 2) then
					call write_tec3D_full_head_z
				end if
			else if(IS3D_binary == 1) then
				! jw
			end if
		end if		

		! jw
		! jw
		! jw
		! jw
	end if
end subroutine open_3D_out_files