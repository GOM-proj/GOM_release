!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will calculate [reference_diff_days], which is the total days
!! between [data_start_year] and ['simulation'_start_year] if simulation starts in a different year.
!! 
subroutine calculate_reference_diff_days
	use mod_global_variables
	implicit none
	
	integer :: year1, year2
	! jw
	
	year1 = data_start_year
	year2 = start_year
	
	reference_diff_days = 0
	do
		if(year1 /= year2) then
			if(MOD(year1,4) /= 0) then ! jw
				reference_diff_days = reference_diff_days + 365
			else ! jw
				reference_diff_days = reference_diff_days + 364
			end if
			
			! jw
			year1 = year1 + 1
		else if(year1 == year2) then
			exit ! jw
		end if
	end do
	
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	reference_diff_days = reference_diff_days - data_time_shift
end subroutine calculate_reference_diff_days