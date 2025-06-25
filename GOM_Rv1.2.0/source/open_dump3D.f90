!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine open_dump3D
	use mod_global_variables
	use mod_file_definition
	implicit none

	character(len=40) :: format_string
	character(len= 4) :: File_num_buff   
	! End of local variables ==================================================!

	format_string = '(I4.4)'	! an integer of width 4 with zeros at the left. For example: 0001 ~ 9999

	! define time conversion factor
	if(IS3D_dump_time == 1) then	! second
		IS3D_dump_time_conv = 86400.0
	else if(IS3D_dump_time == 2) then ! minute
		IS3D_dump_time_conv = 1440.0
	else if(IS3D_dump_time == 3) then ! hour
		IS3D_dump_time_conv = 24.0
	else if(IS3D_dump_time == 4) then ! day
		IS3D_dump_time_conv = 1.0
	end if
	
	! set the first file number
	IS3D_dump_File_num = 1

	write(File_num_buff,format_string) IS3D_dump_File_num
	IS3D_dump_File_name = trim(id_dump3D)//trim(File_num_buff)//'.dat'
	if(IS3D_dump_binary == 0) then
		! ASCII format
		open(pw_dump3D, file = trim(IS3D_dump_File_name), form = 'formatted', status = 'replace')	! dump3D_0000.dat
	else if(IS3D_dump_binary == 1) then
		! Bindary format
		open(pw_dump3D, file = trim(IS3D_dump_File_name), form = 'unformatted', status = 'replace')	! dump3D_0000.dat
	end if

end subroutine open_dump3D