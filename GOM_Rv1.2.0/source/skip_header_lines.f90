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
	character(len = 10) :: line ! read 10 characters for safety
	integer :: EOF
	! End of local variables ==================================================!
	
	do
		read(pw_file,'(A)',iostat=EOF) line  ! read line by line
		
		if(EOF > 0) then
			write(pw_run_log,*) 	'Error during read: ', trim(id_file)
			write(*,*)				'Error during read: ', trim(id_file)
			stop
		else if(EOF < 0) then
			exit
		end if
		
		! Below this line is for if(EOF == 0), =================================!
		! which means FORTRAN successfully read a line
		! If the line contains '!', which is the statement character, in the first column,
		! then skip the rest part of the loop, and go to the next line.
		! This method is more time consuming than line-by-line reading,
		! but it will have an advantage when you want put more comments.

		! if(index(line,"!") == 1) then	! index function returns the location of the looking character.		
		if(index(line,"!") == 1 .or. index(line,"C") == 1 .or. index(line,"c") == 1) then	! index function returns the location of the looking character.
			cycle
		else
			backspace(unit=pw_file)	! send back to the previos step
			exit
		end if		
	end do
end subroutine skip_header_lines