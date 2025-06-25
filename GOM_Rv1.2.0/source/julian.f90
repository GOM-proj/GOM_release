!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! Convert simulation time to julian day  
!! jw, done
!! Actually, the two intented inout values, jday & jday_1900, are only output values, 
!! however, they should be set to inout values since this subroutine is called in many other subroutines.
!! 
subroutine julian(year, month, day, hour, minute, second, jday, jday_1900) 
	use mod_global_variables, only : dp
   implicit none     
      
   integer, intent(in)   :: year, month, day, hour, minute, second
   real(dp),intent(inout) :: jday, jday_1900
   integer :: i  
   ! end of local variables ==================================================!
   
   ! As you can see here, this is not actual Julian day of year, but it is elapsed Julian day of year from the first day of each year.
   ! i.e., "jday" starts from 0.0 not from 1.0
   ! I also explained this a bit more later.
   select case(month)
   	case(1)
   		jday =   0.0_dp
   	case(2)
   		jday =  31.0_dp
   	case(3)
   		jday =  59.0_dp
   	case(4)
   		jday =  90.0_dp
   	case(5)
   		jday = 120.0_dp
   	case(6)
   		jday = 151.0_dp
   	case(7)
   		jday = 181.0_dp
   	case(8)
   		jday = 212.0_dp
   	case(9)
   		jday = 243.0_dp
   	case(10)
   		jday = 273.0_dp
   	case(11)
   		jday = 304.0_dp
   	case(12)
   		jday = 334.0_dp 
   	case default
   		! do nothing
   end select
   
   ! include leap year -------------------------------------------------------!
   if(mod(year,4) == 0) then
   	if(year /= 1900) then
   		! we have one more day in February (28 -> 29 days), 
   		! thus we have to include one more day after February.
   		if(month >= 3) then
   			jday = jday + 1
   		end if
      end if
   end if
	
	! Here, jday is the total elapsed time in [days] at the given simulation start time position: start_year/month/day/hour/minute/second,
	! from the beginning of the given starting start_year
	! i.e., if I set the simulation start year/month/day/hour/minute/second as 2020/02/01 - 00/00/00,
	! jday, which is the simulation starting day from the beginning of 2020 year, will be 31.0 [day]
   ! jday = jday + day + hour/24.0_dp + minute/1440.0_dp + second/86400.0_dp
   ! I use (day-1) since elapsed julian day should start from day 0 not 1.
   jday = jday + (day-1) + hour/24.0_dp + minute/1440.0_dp + second/86400.0_dp
   jday_1900 = jday
   
   ! jday_1900 is the total julian day counted from 1900,
   ! thus we have to add the total julian day from 1900 to the (year-1). 
   do i=1900, year-1
      if(mod(i,4) == 0) then 
         if(i == 1900) then
            jday_1900 = jday_1900+365    ! no leap year in 1900
         else
            jday_1900 = jday_1900+366    ! normal leap year
         end if
      else
         jday_1900 = jday_1900+365
      end if
   end do
end subroutine julian
