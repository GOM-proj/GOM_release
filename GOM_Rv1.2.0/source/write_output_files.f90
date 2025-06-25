!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_output_files
   use mod_global_variables
   use mod_file_definition   
   use mod_netcdf_utils
   use netcdf
	implicit none

	integer :: data_in_2D
 	integer :: data_in_3D
 	integer :: data_in_2D_dump
 	integer :: data_in_3D_dump
	character(len=40) :: format_string
	character(len= 4) :: File_num_buff1, File_num_buff2, File_num_buff3
	character(len= 4) :: File_num_surf_buff, File_num_full_buff
	
	! for netcdf file:
	character(len= 4) :: File_num_buff4
	integer :: netcdf_status
	! End of local variables ==================================================!

	format_string = '(I4.4)'	! an integer of width 4 with zeros at the left. For example: 0001 ~ 9999
	
	!$omp parallel
	!$omp sections

	! write time series data ==================================================!
	!$omp section
	if(tser_station_num > 0) then
      if(mod(it, tser_frequency) == 0) then
      	! this subroutine has both tecplot & vtk writing formats in it:
        	call write_tser_out_files
      end if
   end if
      
	! write 2D output files ===================================================!
	!$omp section
	if(IS2D_switch == 1) then
		data_in_2D = INT(IS2D_File_freq/IS2D_frequency) ! total number of data in 2D file
		if(it >= IS2D_start .and. it <= IS2D_end) then
			IS2D_buff = IS2D_buff + 1 ! this is a new it for output file writing
			if(mod(IS2D_buff,IS2D_frequency) == 0) then	! write 2D file every IS2D_frequency					
				IS2D_buff2 = IS2D_buff2 + 1 ! writing data count from the beginning
				
				! write Tecplot files --------------------------------------------!
				! if(IS2D_format == 1 .or. IS2D_format == 3) then
				if(IS2D_format == 1) then
					! if(IS2D_buff2 /=1 .and. mod(IS2D_buff2,IS2D_File_freq) == 1) then	! this is an old method
					if(IS2D_buff2 /=1 .and. mod(IS2D_buff2,data_in_2D) == 1) then ! this is a new method
						! close previous tec2D_***.dat and create a new (next) tec2D_***.dat
						if(IS2D_binary == 0) then
							close(pw_tec2D)	! close previous tec2D_***.dat
							
							! open a new tec2D_***.dat
							IS2D_File_num = IS2D_File_num + 1
							zone_num_2D = 1	! reset zone number to 1
														
							write(File_num_buff1,format_string) IS2D_File_num		
							IS2D_File_name = trim(id_tec2D)//trim(File_num_buff1)//'.dat'
									
							open(pw_tec2D, file = trim(IS2D_File_name), form = 'formatted', status = 'replace')	! tec2D_***.out	
							
							call write_tec2D_head
						else if(IS2D_binary == 1) then
							! not yet activated
						end if
					else
						! append data to the previous tectplot file
						if(IS2D_binary == 0) then
							zone_num_2D = zone_num_2D + 1
							call write_tec2D_body
						else if(IS2D_binary == 1) then
							! not yet activated
						end if
					end if
				end if

				! write VTK files ------------------------------------------------!
				! if(IS2D_format == 2 .or. IS2D_format == 3) then
				if(IS2D_format == 2) then
					IS2D_vtk_num = IS2D_vtk_num + 1
					if(IS2D_binary == 0) then
						call write_vtk2D
					else if(IS2D_binary == 1) then 
						! not yet included
					end if
				end if
				
				! write Tecplot/VTK files in parallel ----------------------------!
				! This looks ugly, but it is faster with this nesting approach...
				if(IS2D_format == 3) then
					!$omp parallel	sections
					!$omp section ! write tecplot file
					if(IS2D_buff2 /=1 .and. mod(IS2D_buff2,data_in_2D) == 1) then ! this is a new method
						! close previous tec2D_***.dat and create a new (next) tec2D_***.dat
						if(IS2D_binary == 0) then
							close(pw_tec2D)	! close previous tec2D_***.dat
							
							! open a new tec2D_***.dat
							IS2D_File_num = IS2D_File_num + 1
							zone_num_2D = 1	! reset zone number to 1
														
							write(File_num_buff1,format_string) IS2D_File_num		
							IS2D_File_name = trim(id_tec2D)//trim(File_num_buff1)//'.dat'
									
							open(pw_tec2D, file = trim(IS2D_File_name), form = 'formatted', status = 'replace')	! tec2D_***.out	
							
							call write_tec2D_head
						else if(IS2D_binary == 1) then
							! not yet activated
						end if
					else
						! append data to the previous tectplot file
						if(IS2D_binary == 0) then
							zone_num_2D = zone_num_2D + 1
							call write_tec2D_body
						else if(IS2D_binary == 1) then
							! not yet activated
						end if
					end if

					!$omp section ! write vtk file
					IS2D_vtk_num = IS2D_vtk_num + 1
					if(IS2D_binary == 0) then
						call write_vtk2D
					else if(IS2D_binary == 1) then 
						! not yet included
					end if
					!$omp end parallel sections
				end if
			end if
		end if
	end if
	
	
   ! write 3D tecplot files ==================================================!
   !$omp section
	if(IS3D_surf_switch == 1 .or. IS3D_full_switch == 1) then
		data_in_3D = INT(IS3D_File_freq/IS3D_frequency) ! total number of data in 3D file
		if(it >= IS3D_start .and. it <= IS3D_end) then
			IS3D_buff = IS3D_buff + 1
			if(mod(IS3D_buff,IS3D_frequency) == 0) then	! write 3D file every IS3D_frequency
				IS3D_buff2 = IS3D_buff2 + 1
				
				! write Tecplot files -----------------------------------------!
				! if(IS3D_format == 1 .or. IS3D_format == 3) then
				if(IS3D_format == 1) then
					! if(IS3D_buff2 /= 1 .and. mod(IS3D_buff2,IS3D_File_freq) == 1) then	! close previous file and open the next file
					if(IS3D_buff2 /=1 .and. mod(IS3D_buff2,data_in_3D) == 1) then ! this is a new method
						if(IS3D_surf_switch == 1) then
							if(IS3D_binary == 0) then
								close(pw_tec3D_surf)
								
								IS3D_File_num_surf = IS3D_File_num_surf + 1
								zone_num_3D_surf = 1	! reset zone number to 1
								
								write(File_num_surf_buff,format_string) IS3D_File_num_surf
								IS3D_File_name_surf = trim(id_tec3D_surf)//trim(File_num_surf_buff)//'.dat'
								
								open(pw_tec3D_surf, file = trim(IS3D_File_name_surf), form = 'formatted', status = 'replace')	! tec3D_surf_***.out
								call write_tec3D_surf_head		
							else if(IS3D_binary == 1) then
								! not yet activated
							end if
						end if
						
						if(IS3D_full_switch == 1) then
							if(IS3D_binary == 0) then
								close(pw_tec3D_full)
								
								IS3D_File_num_full = IS3D_File_num_full + 1
								zone_num_3D_full = 1	! reset zone number to 1
								
								write(File_num_full_buff,format_string) IS3D_File_num_full
								IS3D_File_name_full = trim(id_tec3D_full)//trim(File_num_full_buff)//'.dat'
								
								open(pw_tec3D_full, file = trim(IS3D_File_name_full), form = 'formatted', status = 'replace')	! tec3D_full_***.out
								if(IS3D_grid_format == 1) then
									call write_tec3D_full_head_zto_sigma
								else if(IS3D_grid_format == 2) then
									call write_tec3D_full_head_z
								end if
							else if(IS3D_binary == 1) then
								! not yet activated
							end if
						end if
					else	! append to the previous tectplot file
						if(IS3D_surf_switch == 1) then
							if(IS3D_binary == 0) then
								zone_num_3D_surf = zone_num_3D_surf + 1
								call write_tec3D_surf_body
							else if(IS3D_binary == 1) then
								! not yet activated
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
								! not yet activated
							end if							
						end if					
					end if
				end if ! if(IS3D_format == 1) then

				! write VTK files ---------------------------------------------!
				! if(IS3D_format == 2 .or. IS3D_format == 3) then
				if(IS3D_format == 2) then
					IS3D_vtk_num = IS3D_vtk_num + 1
					if(IS3D_binary == 0) then
						call write_vtk3D
					else if(IS3D_binary == 1) then
						! not yet included
					end if
				end if
				
				! write both tecplot/vtk files -----------------------------------!
				if(IS3D_format == 3) then
					!$omp parallel	sections
					!$omp section ! write tecplot file
					! if(IS3D_buff2 /= 1 .and. mod(IS3D_buff2,IS3D_File_freq) == 1) then	! close previous file and open the next file
					if(IS3D_buff2 /=1 .and. mod(IS3D_buff2,data_in_3D) == 1) then ! this is a new method
						if(IS3D_surf_switch == 1) then
							if(IS3D_binary == 0) then
								close(pw_tec3D_surf)
								
								IS3D_File_num_surf = IS3D_File_num_surf + 1
								zone_num_3D_surf = 1	! reset zone number to 1
								
								write(File_num_surf_buff,format_string) IS3D_File_num_surf
								IS3D_File_name_surf = trim(id_tec3D_surf)//trim(File_num_surf_buff)//'.dat'
								
								open(pw_tec3D_surf, file = trim(IS3D_File_name_surf), form = 'formatted', status = 'replace')	! tec3D_surf_***.out
								call write_tec3D_surf_head		
							else if(IS3D_binary == 1) then
								! not yet activated
							end if
						end if
						
						if(IS3D_full_switch == 1) then
							if(IS3D_binary == 0) then
								close(pw_tec3D_full)
								
								IS3D_File_num_full = IS3D_File_num_full + 1
								zone_num_3D_full = 1	! reset zone number to 1
								
								write(File_num_full_buff,format_string) IS3D_File_num_full
								IS3D_File_name_full = trim(id_tec3D_full)//trim(File_num_full_buff)//'.dat'
								
								open(pw_tec3D_full, file = trim(IS3D_File_name_full), form = 'formatted', status = 'replace')	! tec3D_full_***.out
								if(IS3D_grid_format == 1) then
									call write_tec3D_full_head_zto_sigma
								else if(IS3D_grid_format == 2) then
									call write_tec3D_full_head_z
								end if
							else if(IS3D_binary == 1) then
								! not yet activated
							end if
						end if
					else	! append to the previous tectplot file
						if(IS3D_surf_switch == 1) then
							if(IS3D_binary == 0) then
								zone_num_3D_surf = zone_num_3D_surf + 1
								call write_tec3D_surf_body
							else if(IS3D_binary == 1) then
								! not yet activated
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
								! not yet activated
							end if							
						end if					
					end if
					
					!$omp section ! write vtk file				
					IS3D_vtk_num = IS3D_vtk_num + 1
					if(IS3D_binary == 0) then
						call write_vtk3D
					else if(IS3D_binary == 1) then
						! not yet included
					end if
					!$omp end parallel sections				
				end if ! if(IS3D_format == 3) then
			end if ! if(mod(IS3D_buff,IS3D_frequency) == 0) then	! write 3D file every IS3D_frequency
		end if
	end if
	
	! write 2D dump files =====================================================!
	!$omp section
	if(IS2D_dump_switch == 1) then
		data_in_2D_dump = INT(IS2D_dump_File_freq/IS2D_dump_frequency) ! total number of data in each 2D dump file
		if(it >= IS2D_dump_start .and. it <= IS2D_dump_end) then
			IS2D_dump_buff = IS2D_dump_buff + 1 ! this is a new it for output file writing
			if(mod(IS2D_dump_buff,IS2D_dump_frequency) == 0) then	! write 2D file every IS2D_dump_frequency
				IS2D_dump_buff2 = IS2D_dump_buff2 + 1
				if(IS2D_dump_buff2 /=1 .and. mod(IS2D_dump_buff2,data_in_2D_dump) == 1) then ! this is a new method
					! close previous dump2D_***.dat and create a new (next) dump2D_***.dat
					close(pw_dump2D)
					
					! open a new dump2D_***.dat
					! This is the identical routine as in open_dump2D.f90, but explicitly show here.
					IS2D_dump_File_num = IS2D_dump_File_num + 1
					write(File_num_buff2,format_string) IS2D_dump_File_num
					
					IS2D_dump_File_name = trim(id_dump2D)//trim(File_num_buff2)//'.dat' ! dump2D_****.dat
					if(IS2D_dump_binary == 0) then
						! ASCII format
						open(pw_dump2D, file = trim(IS2D_dump_File_name), form = 'formatted', status = 'replace')	! dump2D_***.dat
					else if(IS2D_dump_binary == 1) then
						! Binary format
						open(pw_dump2D, file = trim(IS2D_dump_File_name), form = 'unformatted', status = 'replace')	! dump2D_***.dat
					end if
				end if
				
				call write_dump2D
			end if
		end if
	end if
	
	! write 3D dump files =====================================================!
	!$omp section
	if(IS3D_dump_switch == 1) then
		data_in_3D_dump = INT(IS3D_dump_File_freq/IS3D_dump_frequency) ! total number of data in each 3D dump file
		if(it >= IS3D_dump_start .and. it <= IS3D_dump_end) then
			IS3D_dump_buff = IS3D_dump_buff + 1 ! this is a new it for output file writing
			if(mod(IS3D_dump_buff,IS3D_dump_frequency) == 0) then	! write 3D file every IS3D_dump_frequency
				IS3D_dump_buff2 = IS3D_dump_buff2 + 1
				if(IS3D_dump_buff2 /=1 .and. mod(IS3D_dump_buff2,data_in_3D_dump) == 1) then ! this is a new method
					! close previous dump3D_***.dat and create a new (next) dump3D_***.dat
					close(pw_dump3D)
					
					! open a new dump3D_***.dat
					! This is the identical routine as in open_dump3D.f90, but explicitly show here.
					IS3D_dump_File_num = IS3D_dump_File_num + 1
					write(File_num_buff3,format_string) IS3D_dump_File_num
					
					IS3D_dump_File_name = trim(id_dump3D)//trim(File_num_buff3)//'.dat' ! dump3D_****.dat
					if(IS3D_dump_binary == 0) then
						! ASCII format
						open(pw_dump3D, file = trim(IS3D_dump_File_name), form = 'formatted', status = 'replace')	! dump3D_***.dat
					else if(IS3D_dump_binary == 1) then
						! Binary format
						open(pw_dump3D, file = trim(IS3D_dump_File_name), form = 'unformatted', status = 'replace')	! dump3D_***.dat
					end if
				end if
				
				call write_dump3D
			end if
		end if
	end if
		
	! write netcdf output file ================================================!
	!$omp section
	if(netcdf_switch == 1) then
		! This is calculated in prepare_gom.f90
		! netcdf_maxtime = INT(netcdf_File_Freq/netcdf_frequency)
		if(it >= netcdf_start .and. it <= netcdf_end) then
			netcdf_buff = netcdf_buff + 1 ! this is a new it for output file writing
			
			if(mod(netcdf_buff,netcdf_frequency) == 0) then ! write netcdf data every netcdf_frequency
				netcdf_buff2 = netcdf_buff2 + 1 ! total number of netcdf writing count from the beginning

				! at every new netcdf file generation time, we have to close the previous file and write header.
				! but, at the first time, we do not have a file to close, and that is why I include one more "if" statement
				if(mod(netcdf_buff2,netcdf_maxtime) == 1) then
					! close previous netcdf_out_***.nc and create a new (next) netcdf_out_***.nc
					
					! close previous file only after the 1st netcdf file writing
					if(netcdf_buff2 > 1) then
						call netcdf_check(nf90_close(ncid))
					end if
										
					! open a new netcdf_out_***.nc
					netcdf_File_num = Netcdf_File_num + 1
					write(File_num_buff4,format_string) netcdf_File_num
					netcdf_File_name = trim(id_netcdf_out)//trim(File_num_buff4)//'.nc'

					netcdf_status = nf90_create(netcdf_File_name, nf90_clobber, ncid) ! this will create: ./output/netcdf_out_***.nc
					call netcdf_check(netcdf_status)
					
					! write header (dimension & attribute) information in the new netcdf file
					call write_netcdf_head
					
					! reset netcdf_buff3 to 0
					! Note: netcdf_buff3 is the time position number (i.e., i'th time domain in the array) the netcdf file
					netcdf_buff3 = 0
				end if

				! append data to the opened (or previous) netcdf file
				netcdf_buff3 = netcdf_buff3+1
				call write_netcdf_body
			end if
		end if
	end if
	
	
	!$omp end sections
	!$omp end parallel
end subroutine write_output_files