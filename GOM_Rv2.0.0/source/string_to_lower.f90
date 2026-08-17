!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will change given characters to lower cases.
!! 
subroutine string_to_lower(string,lower)
	implicit none
	
	character(len=5), intent(in) :: string
	character(len=5), intent(out) :: lower
	integer :: i
	! jw
	
	do i=1,len(string)
		if(string(i:i) >= "A" .and. string(i:i) <= "Z") then
			lower(i:i) = achar(iachar(string(i:i)) + 32)
		else
			lower(i:i) = string(i:i)
		end if
	end do
end subroutine string_to_lower