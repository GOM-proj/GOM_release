!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_output_files
   use mod_global_variables
   use mod_file_definition   
	implicit none

	integer :: IS2D_buff = 0, IS2D_buff2 = 0, data_in_2D
 	integer :: IS3D_buff = 0, IS3D_buff2 = 0, data_in_3D
 	integer :: IS2D_dump_buff = 0, IS2D_dump_buff2 = 0, data_in_2D_dump
 	integer :: IS3D_dump_buff = 0, IS3D_dump_buff2 = 0, data_in_3D_dump
	character(len=40) :: format_string
	character(len= 4) :: File_num_buff1, File_num_buff2, File_num_buff3
	character(len= 4) :: File_num_surf_buff, File_num_full_buff

	format_string = '(I4.4)'	! jw
	! jw
	
	!$omp parallel	
	!$omp sections

	! jw
	!$omp section
	if(tser_station_num > 0) then
      if(mod(it, tser_frequency) == 0) then
      	! jw
        	call write_tser_out_files
      end if
   end if
      
	! jw
	!$omp section
	if(IS2D_switch == 1) then
		data_in_2D = INT(IS2D_File_freq/IS2D_frequency) ! jw
		if(it >= IS2D_start .and. it <= IS2D_end) then
			IS2D_buff = IS2D_buff + 1 ! jw
			if(mod(IS2D_buff,IS2D_frequency) == 0) then	! jw
				IS2D_buff2 = IS2D_buff2 + 1 ! jw
				
				! jw
				! jw
				if(IS2D_format == 1) then
					! jw
					if(IS2D_buff2 /=1 .and. mod(IS2D_buff2,data_in_2D) == 1) then ! jw
						! jw
						if(IS2D_binary == 0) then
							close(pw_tec2D)	! jw
							
							! jw
							IS2D_File_num = IS2D_File_num + 1
							zone_num_2D = 1	! jw
														
							write(File_num_buff1,format_string) IS2D_File_num		
							IS2D_File_name = trim(id_tec2D)//trim(File_num_buff1)//'.dat'
									
							open(pw_tec2D, file = trim(IS2D_File_name), form = 'formatted', status = 'replace')	! jw
							
							call write_tec2D_head
						else if(IS2D_binary == 1) then
							! jw
						end if
					else
						! jw
						if(IS2D_binary == 0) then
							zone_num_2D = zone_num_2D + 1
							call write_tec2D_body
						else if(IS2D_binary == 1) then
							! jw
						end if
					end if
				end if

				! jw
				! jw
				if(IS2D_format == 2) then
					IS2D_vtk_num = IS2D_vtk_num + 1
					if(IS2D_binary == 0) then
						call write_vtk2D
					else if(IS2D_binary == 1) then 
						! jw
					end if
				end if
				
				! jw
				! jw
				if(IS2D_format == 3) then
					!$omp parallel	sections
					!$omp section ! jw
					if(IS2D_buff2 /=1 .and. mod(IS2D_buff2,data_in_2D) == 1) then ! jw
						! jw
						if(IS2D_binary == 0) then
							close(pw_tec2D)	! jw
							
							! jw
							IS2D_File_num = IS2D_File_num + 1
							zone_num_2D = 1	! jw
														
							write(File_num_buff1,format_string) IS2D_File_num		
							IS2D_File_name = trim(id_tec2D)//trim(File_num_buff1)//'.dat'
									
							open(pw_tec2D, file = trim(IS2D_File_name), form = 'formatted', status = 'replace')	! jw
							
							call write_tec2D_head
						else if(IS2D_binary == 1) then
							! jw
						end if
					else
						! jw
						if(IS2D_binary == 0) then
							zone_num_2D = zone_num_2D + 1
							call write_tec2D_body
						else if(IS2D_binary == 1) then
							! jw
						end if
					end if

					!$omp section ! jw
					IS2D_vtk_num = IS2D_vtk_num + 1
					if(IS2D_binary == 0) then
						call write_vtk2D
					else if(IS2D_binary == 1) then 
						! jw
					end if
					!$omp end parallel sections
				end if
			end if
		end if
	end if
	
	
   ! jw
   !$omp section
	if(IS3D_surf_switch == 1 .or. IS3D_full_switch == 1) then
		data_in_3D = INT(IS3D_File_freq/IS3D_frequency) ! jw
		if(it >= IS3D_start .and. it <= IS3D_end) then
			IS3D_buff = IS3D_buff + 1
			if(mod(IS3D_buff,IS3D_frequency) == 0) then	! jw
				IS3D_buff2 = IS3D_buff2 + 1
				
				! jw
				! jw
				if(IS3D_format == 1) then
					! jw
					if(IS3D_buff2 /=1 .and. mod(IS3D_buff2,data_in_3D) == 1) then ! jw
						if(IS3D_surf_switch == 1) then
							if(IS3D_binary == 0) then
								close(pw_tec3D_surf)
								
								IS3D_File_num_surf = IS3D_File_num_surf + 1
								zone_num_3D_surf = 1	! jw
								
								write(File_num_surf_buff,format_string) IS3D_File_num_surf
								IS3D_File_name_surf = trim(id_tec3D_surf)//trim(File_num_surf_buff)//'.dat'
								
								open(pw_tec3D_surf, file = trim(IS3D_File_name_surf), form = 'formatted', status = 'replace')	! jw
								call write_tec3D_surf_head		
							else if(IS3D_binary == 1) then
								! jw
							end if
						end if
						
						if(IS3D_full_switch == 1) then
							if(IS3D_binary == 0) then
								close(pw_tec3D_full)
								
								IS3D_File_num_full = IS3D_File_num_full + 1
								zone_num_3D_full = 1	! jw
								
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
					else	! jw
						if(IS3D_surf_switch == 1) then
							if(IS3D_binary == 0) then
								zone_num_3D_surf = zone_num_3D_surf + 1
								call write_tec3D_surf_body
							else if(IS3D_binary == 1) then
								! jw
							end if
						end if
						if(IS3D_full_switch == 1) then
							if(IS3D_binary == 0) then
								zone_num_3D_full = zone_num_3D_full + 1
								if(IS3D_grid_format == 1) then
									call write_tec3D_full_body_zto_sigma
								else if(IS3D_grid_format == 2) then
									call write_tec3D_full_body_z
								end if
							else if(IS3D_binary == 1) then
								! jw
							end if							
						end if					
					end if
				end if ! jw

				! jw
				! jw
				if(IS3D_format == 2) then
					IS3D_vtk_num = IS3D_vtk_num + 1
					if(IS3D_binary == 0) then
						call write_vtk3D
					else if(IS3D_binary == 1) then
						! jw
					end if
				end if
				
				! jw
				if(IS3D_format == 3) then
					!$omp parallel	sections
					!$omp section ! jw
					! jw
					if(IS3D_buff2 /=1 .and. mod(IS3D_buff2,data_in_3D) == 1) then ! jw
						if(IS3D_surf_switch == 1) then
							if(IS3D_binary == 0) then
								close(pw_tec3D_surf)
								
								IS3D_File_num_surf = IS3D_File_num_surf + 1
								zone_num_3D_surf = 1	! jw
								
								write(File_num_surf_buff,format_string) IS3D_File_num_surf
								IS3D_File_name_surf = trim(id_tec3D_surf)//trim(File_num_surf_buff)//'.dat'
								
								open(pw_tec3D_surf, file = trim(IS3D_File_name_surf), form = 'formatted', status = 'replace')	! jw
								call write_tec3D_surf_head		
							else if(IS3D_binary == 1) then
								! jw
							end if
						end if
						
						if(IS3D_full_switch == 1) then
							if(IS3D_binary == 0) then
								close(pw_tec3D_full)
								
								IS3D_File_num_full = IS3D_File_num_full + 1
								zone_num_3D_full = 1	! jw
								
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
					else	! jw
						if(IS3D_surf_switch == 1) then
							if(IS3D_binary == 0) then
								zone_num_3D_surf = zone_num_3D_surf + 1
								call write_tec3D_surf_body
							else if(IS3D_binary == 1) then
								! jw
							end if
						end if
						if(IS3D_full_switch == 1) then
							if(IS3D_binary == 0) then
								zone_num_3D_full = zone_num_3D_full + 1
								if(IS3D_grid_format == 1) then
									call write_tec3D_full_body_zto_sigma
								else if(IS3D_grid_format == 2) then
									call write_tec3D_full_body_z
								end if
							else if(IS3D_binary == 1) then
								! jw
							end if							
						end if					
					end if
					
					!$omp section ! jw
					IS3D_vtk_num = IS3D_vtk_num + 1
					if(IS3D_binary == 0) then
						call write_vtk3D
					else if(IS3D_binary == 1) then
						! jw
					end if
					!$omp end parallel sections				
				end if ! jw
			end if ! jw
		end if
	end if
	
	! jw
	!$omp section
	if(IS2D_dump_switch == 1) then
		data_in_2D_dump = INT(IS2D_dump_File_freq/IS2D_dump_frequency) ! jw
		if(it >= IS2D_dump_start .and. it <= IS2D_dump_end) then
			IS2D_dump_buff = IS2D_dump_buff + 1 ! jw
			if(mod(IS2D_dump_buff,IS2D_dump_frequency) == 0) then	! jw
				IS2D_dump_buff2 = IS2D_dump_buff2 + 1
				if(IS2D_dump_buff2 /=1 .and. mod(IS2D_dump_buff2,data_in_2D_dump) == 1) then ! jw
					! jw
					close(pw_dump2D)
					
					! jw
					! jw
					IS2D_dump_File_num = IS2D_dump_File_num + 1
					write(File_num_buff2,format_string) IS2D_dump_File_num
					
					IS2D_dump_File_name = trim(id_dump2D)//trim(File_num_buff2)//'.dat' ! jw
					if(IS2D_dump_binary == 0) then
						! jw
						open(pw_dump2D, file = trim(IS2D_dump_File_name), form = 'formatted', status = 'replace')	! jw
					else if(IS2D_dump_binary == 1) then
						! jw
						open(pw_dump2D, file = trim(IS2D_dump_File_name), form = 'unformatted', status = 'replace')	! jw
					end if
				end if
				
				call write_dump2D
			end if
		end if
	end if
	
	! jw
	!$omp section
	if(IS3D_dump_switch == 1) then
		data_in_3D_dump = INT(IS3D_dump_File_freq/IS3D_dump_frequency) ! jw
		if(it >= IS3D_dump_start .and. it <= IS3D_dump_end) then
			IS3D_dump_buff = IS3D_dump_buff + 1 ! jw
			if(mod(IS3D_dump_buff,IS3D_dump_frequency) == 0) then	! jw
				IS3D_dump_buff2 = IS3D_dump_buff2 + 1
				if(IS3D_dump_buff2 /=1 .and. mod(IS3D_dump_buff2,data_in_3D_dump) == 1) then ! jw
					! jw
					close(pw_dump3D)
					
					! jw
					! jw
					IS3D_dump_File_num = IS3D_dump_File_num + 1
					write(File_num_buff3,format_string) IS3D_dump_File_num
					
					IS3D_dump_File_name = trim(id_dump3D)//trim(File_num_buff3)//'.dat' ! jw
					if(IS3D_dump_binary == 0) then
						! jw
						open(pw_dump3D, file = trim(IS3D_dump_File_name), form = 'formatted', status = 'replace')	! jw
					else if(IS3D_dump_binary == 1) then
						! jw
						open(pw_dump3D, file = trim(IS3D_dump_File_name), form = 'unformatted', status = 'replace')	! jw
					end if
				end if
				
				call write_dump3D
			end if
		end if
	end if
		
	
	!$omp end sections
	!$omp end parallel
end subroutine write_output_files