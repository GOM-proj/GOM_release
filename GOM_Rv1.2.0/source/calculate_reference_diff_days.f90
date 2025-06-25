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
	! End of local variables ==================================================!
	
	year1 = data_start_year
	year2 = start_year
	
	reference_diff_days = 0
	do
		if(year1 /= year2) then
			if(MOD(year1,4) /= 0) then ! normal years
				reference_diff_days = reference_diff_days + 365
			else ! leap years
				reference_diff_days = reference_diff_days + 364
			end if
			
			! increase data_start_year to the next year
			year1 = year1 + 1
		else if(year1 == year2) then
			exit ! exit do loop
		end if
	end do
	
	! Additional time shift in [day]
	! For example, let's consider following:
	! 		if data_start_year is 2000 and assume that there are 365 days in 2000,
	! 		and the simulation starts at February 01, 2001
	! 		then, we can say that the data actually starts 334 days (365-31) ahead the simulation start day.
	! 		Here, 31 days is the [data_time_shift] in day from 2000/01/01
	! 		Finally, in each time interpolation routine, 
	! 		the interpolation reference time is shifted [reference_diff_days] to point exact time in the simulation time frame.
	reference_diff_days = reference_diff_days - data_time_shift
end subroutine calculate_reference_diff_days