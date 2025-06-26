!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! This will skp header lines which contain '!', 'C', or 'c' at the first column of each line
!!
subroutine skip_header_lines(pw_file, id_file)
	use mod_file_definition
	
	implicit none
	integer,intent(in) :: pw_file
	character(len=100),intent(in) :: id_file
	character(len = 10) :: line ! jw
	integer :: EOF
	! jw
	
	do
		read(pw_file,'(A)',iostat=EOF) line  ! jw
		
		if(EOF > 0) then
			write(pw_run_log,*) 	'Error during read: ', trim(id_file)
			write(*,*)				'Error during read: ', trim(id_file)
			stop
		else if(EOF < 0) then
			exit
		end if
		
		! jw
		! jw
		! jw
		! jw
		! jw
		! jw

		! jw
		if(index(line,"!") == 1 .or. index(line,"C") == 1 .or. index(line,"c") == 1) then	! jw
			cycle
		else
			backspace(unit=pw_file)	! jw
			exit
		end if		
	end do
end subroutine skip_header_lines