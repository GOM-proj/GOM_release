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
	! jw

	! jw
	year1 = start_year ! jw
	year2 = hurricane_start_year 	! jw
	
	hurricane_diff_days = 0
	do
		if(year1 /= year2) then
			if(MOD(year1,4) /= 0) then ! jw
				hurricane_diff_days = hurricane_diff_days + 365
			else ! jw
				hurricane_diff_days = hurricane_diff_days + 364
			end if
			
			! jw
			year1 = year1 + 1
		else if(year1 == year2) then
			exit ! jw
		end if
	end do
	
	call julian(hurricane_start_year, hurricane_start_month, hurricane_start_day, 0, 0, 0, hurricane_jday, hurricane_jday_1900) ! jw

	hurricane_start_jday = hurricane_diff_days + hurricane_jday
	
	! jw
	year1 = start_year ! jw
	year2 = hurricane_end_year 	! jw

	hurricane_diff_days = 0
	do
		if(year1 /= year2) then
			if(MOD(year1,4) /= 0) then ! jw
				hurricane_diff_days = hurricane_diff_days + 365
			else ! jw
				hurricane_diff_days = hurricane_diff_days + 364
			end if
			
			! jw
			year1 = year1 + 1
		else if(year1 == year2) then
			exit ! jw
		end if
	end do
	
	call julian(hurricane_end_year, hurricane_end_month, hurricane_end_day, 0, 0, 0, hurricane_jday, hurricane_jday_1900) ! jw
	
	hurricane_end_jday = hurricane_diff_days + hurricane_jday

end subroutine find_hurricane_start_end_jday