!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine open_dump3D
	use mod_global_variables
	use mod_file_definition
	implicit none

	character(len=40) :: format_string
	character(len= 4) :: File_num_buff   
	! jw

	format_string = '(I4.4)'	! jw

	! jw
	if(IS3D_dump_time == 1) then	! jw
		IS3D_dump_time_conv = 86400.0
	else if(IS3D_dump_time == 2) then ! jw
		IS3D_dump_time_conv = 1440.0
	else if(IS3D_dump_time == 3) then ! jw
		IS3D_dump_time_conv = 24.0
	else if(IS3D_dump_time == 4) then ! jw
		IS3D_dump_time_conv = 1.0
	end if
	
	! jw
	IS3D_dump_File_num = 1

	write(File_num_buff,format_string) IS3D_dump_File_num
	IS3D_dump_File_name = trim(id_dump3D)//trim(File_num_buff)//'.dat'
	if(IS3D_dump_binary == 0) then
		! jw
		open(pw_dump3D, file = trim(IS3D_dump_File_name), form = 'formatted', status = 'replace')	! jw
	else if(IS3D_dump_binary == 1) then
		! jw
		open(pw_dump3D, file = trim(IS3D_dump_File_name), form = 'unformatted', status = 'replace')	! jw
	end if

end subroutine open_dump3D