!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This converts julian day to local time, 
!! e.g. if start_year (C3 in main.inp) = 2005 and julian day = 190.0 (in 2005),
!! then, local time will be: 07/09/2005 00:00
!! jw, done
subroutine julian_to_localtime(year, month, day, hour, minute, jday)  
	use mod_global_variables, only : dp
   
   implicit none 
   integer, intent(inout) :: year, month, day, hour, minute
   real(dp), intent(inout) :: jday
	
	integer :: year1
   real(dp) :: day_temp = 0.0 						! jw
   real(dp), dimension(12) :: temp_month = 0.0 	! jw
   ! jw

   year1 = year

   if(mod(year1,4) == 0) then
      if(year1 /= 1900) then
         if(month >= 3) then
         	temp_month( 1) =   0.0_dp
            temp_month( 2) =  32.0_dp
            temp_month( 3) =  60.0_dp
            temp_month( 4) =  91.0_dp
            temp_month( 5) = 121.0_dp
            temp_month( 6) = 152.0_dp
            temp_month( 7) = 182.0_dp
            temp_month( 8) = 213.0_dp
            temp_month( 9) = 244.0_dp
            temp_month(10) = 274.0_dp
            temp_month(11) = 305.0_dp
            temp_month(12) = 335.0_dp
            
            ! jw
            if (jday <= 31.0                     ) month =  1
            if (jday >  31.0 .and. jday <=  60.0 ) month =  2
            if (jday >  60.0 .and. jday <=  91.0 ) month =  3
            if (jday >  91.0 .and. jday <= 121.0 ) month =  4
            if (jday > 121.0 .and. jday <= 152.0 ) month =  5
            if (jday > 152.0 .and. jday <= 182.0 ) month =  6
            if (jday > 182.0 .and. jday <= 213.0 ) month =  7
            if (jday > 213.0 .and. jday <= 244.0 ) month =  8
            if (jday > 244.0 .and. jday <= 274.0 ) month =  9
            if (jday > 274.0 .and. jday <= 305.0 ) month = 10
            if (jday > 305.0 .and. jday <= 335.0 ) month = 11
            if (jday > 335.0                     ) month = 12
         end if
      end if
   else
      temp_month( 1) =   0.0_dp
      temp_month( 2) =  31.0_dp
      temp_month( 3) =  59.0_dp
      temp_month( 4) =  90.0_dp
      temp_month( 5) = 120.0_dp
      temp_month( 6) = 151.0_dp
      temp_month( 7) = 181.0_dp
      temp_month( 8) = 212.0_dp
      temp_month( 9) = 243.0_dp
      temp_month(10) = 273.0_dp
      temp_month(11) = 304.0_dp
      temp_month(12) = 334.0_dp
      
      ! jw
      if (jday <= 31.0                     ) month =  1
      if (jday >  31.0 .and. jday <=  59.0 ) month =  2
      if (jday >  59.0 .and. jday <=  90.0 ) month =  3
      if (jday >  90.0 .and. jday <= 120.0 ) month =  4
      if (jday > 120.0 .and. jday <= 151.0 ) month =  5
      if (jday > 151.0 .and. jday <= 181.0 ) month =  6
      if (jday > 181.0 .and. jday <= 212.0 ) month =  7
      if (jday > 212.0 .and. jday <= 243.0 ) month =  8
      if (jday > 243.0 .and. jday <= 273.0 ) month =  9
      if (jday > 273.0 .and. jday <= 304.0 ) month = 10
      if (jday > 304.0 .and. jday <= 334.0 ) month = 11
      if (jday > 334.0                     ) month = 12
   end if

   day_temp = jday - temp_month(month)   
   
   ! jw
   ! jw
   ! jw
   ! jw
   day = int(day_temp) + 1
   
! jw
! jw
! jw
   
   hour = int((day_temp+1 - day) *24.0_dp)
 	minute = int(((day_temp+1 - day) * 24.0_dp - hour) * 60.0_dp)
 	
end subroutine julian_to_localtime
