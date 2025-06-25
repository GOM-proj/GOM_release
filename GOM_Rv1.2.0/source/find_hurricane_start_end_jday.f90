!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will find [hurricane_start_jday] and [hurricane_end_jday]
!! 
subroutine find_hurricane_start_end_jday
	use mod_global_variables
	implicit none
	
	integer :: year1, year2, hurricane_diff_days
	real(dp):: hurricane_jday, hurricane_jday_1900
	! End of local variables ==================================================!

	! find [hurricane_start_jday] =============================================!
	year1 = start_year ! simulation start year
	year2 = hurricane_start_year 	! hurricane_ser.inp start year (hurricane_ser.inp should start from this year)
	
	hurricane_diff_days = 0
	do
		if(year1 /= year2) then
			if(MOD(year1,4) /= 0) then ! normal years
				hurricane_diff_days = hurricane_diff_days + 365
			else ! leap years
				hurricane_diff_days = hurricane_diff_days + 364
			end if
			
			! increase data_start_year to the next year
			year1 = year1 + 1
		else if(year1 == year2) then
			exit ! exit do loop
		end if
	end do
	
	call julian(hurricane_start_year, hurricane_start_month, hurricane_start_day, 0, 0, 0, hurricane_jday, hurricane_jday_1900) ! hurricane_jday and hurricane_jday_1900 are the outputs

	hurricane_start_jday = hurricane_diff_days + hurricane_jday
	
	! find [hurricane_end_jday] ===============================================!
	year1 = start_year ! simulation start year
	year2 = hurricane_end_year 	! hurricane simulatioin end year (no matter how long the hurricane_ser.inp is)

	hurricane_diff_days = 0
	do
		if(year1 /= year2) then
			if(MOD(year1,4) /= 0) then ! normal years
				hurricane_diff_days = hurricane_diff_days + 365
			else ! leap years
				hurricane_diff_days = hurricane_diff_days + 364
			end if
			
			! increase data_start_year to the next year
			year1 = year1 + 1
		else if(year1 == year2) then
			exit ! exit do loop
		end if
	end do
	
	call julian(hurricane_end_year, hurricane_end_month, hurricane_end_day, 0, 0, 0, hurricane_jday, hurricane_jday_1900) ! hurricane_jday and hurricane_jday_1900 are the outputs
	
	hurricane_end_jday = hurricane_diff_days + hurricane_jday

end subroutine find_hurricane_start_end_jday