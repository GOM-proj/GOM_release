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
   ! jw
   
   ! jw
   ! jw
   ! jw
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
   		! jw
   end select
   
   ! jw
   if(mod(year,4) == 0) then
   	if(year /= 1900) then
   		! jw
   		! jw
   		if(month >= 3) then
   			jday = jday + 1
   		end if
      end if
   end if
	
	! jw
	! jw
	! jw
	! jw
   ! jw
   ! jw
   jday = jday + (day-1) + hour/24.0_dp + minute/1440.0_dp + second/86400.0_dp
   jday_1900 = jday
   
   ! jw
   ! jw
   do i=1900, year-1
      if(mod(i,4) == 0) then 
         if(i == 1900) then
            jday_1900 = jday_1900+365    ! jw
         else
            jday_1900 = jday_1900+366    ! jw
         end if
      else
         jday_1900 = jday_1900+365
      end if
   end do
end subroutine julian
