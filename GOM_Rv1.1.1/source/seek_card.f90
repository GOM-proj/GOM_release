!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will find specific card numbers in input files,
!! and locate the pointer to the first line of the data in the card number.
!! 
subroutine seek_card(card_num)
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	character(len=5),intent(inout) :: card_num	! jw
	character(len=5) :: line
	character(len=5) :: string
	character(len=5) :: lower
	integer :: EOF         ! jw
	! jw
	
	! jw
	string = card_num
	lower = ""
	call string_to_lower(string,lower)
	card_num = lower

	do
		read(pw_main_inp,'(A)',iostat=EOF) line ! jw
		if(EOF > 0) then
			write(*,*) 'Error during read main.inp'
			stop 'seek_card.f90, Error #1'
		else if(EOF < 0) then
			exit
		end if
		
		string = line
		lower = ""
		call string_to_lower(string,lower)
		line = lower		
		
		if(index(trim(line),trim(card_num)) /= 0) then
			exit
		else
			cycle
		end if		
	end do

	! jw
	! jw
	! jw
	do
		read(pw_main_inp,'(A)',iostat=EOF) line  ! jw
		
		if(EOF > 0) then
			write(*,*) 'Error during read main.inp'
			write(*,*) 'Card # = ', card_num
			stop 'seek_card.f90, Error #2'
		else if(EOF < 0) then
			exit
		end if
		
		if(index(line,"!") == 1) then
			cycle
		else if(index(line,"C") == 1) then
			cycle
		else if(index(line,"c") == 1) then
			cycle
		else
			backspace(unit=pw_main_inp)
			exit
		end if		
	end do
	! jw
	! jw
end subroutine seek_card