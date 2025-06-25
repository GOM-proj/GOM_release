!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Define spinup function for:
!!		(1) tidal boundary elevation forcing
!!		(2) wind and pressure forcing
!!		(3) tidal potential forcing
!! 
subroutine setup_spinup
   use mod_global_variables
   implicit none
   real(dp) :: spinup_period
   ! end of local variables ==================================================!
   
   ! We can use tanh function to slowly increase the value:
   ! 		tanh(2.0) = 1.0 (0.964)
   ! 		tanh(2.0 * t/T) = 0.0 ~ 1.0
	! 		where T is the spinpu period, and t is the elapsed time)
   ! Spinup_periods are given in [day], thus it should be converted to [sec]
   ! 
   if(tide_spinup == 1) then
   	spinup_period = tide_spinup_period * 86400.0
   	if(elapsed_time < spinup_period) then
      	spinup_function_tide = tanh(2.0*elapsed_time/spinup_period) 
      else
      	spinup_function_tide = 1.0_dp
      end if
   else
      spinup_function_tide = 1.0_dp
   end if

   if(baroclinic_flag == 1 .and. baroclinic_spinup == 1) then
   	spinup_period = baroclinic_spinup_period * 86400.0
   	if(elapsed_time < spinup_period) then
      	spinup_function_baroclinic = tanh(2.0*elapsed_time/spinup_period)
      else
      	spinup_function_baroclinic = 1.0_dp
      end if
   else
      spinup_function_baroclinic = 1.0_dp
   end if

   if(wind_flag == 1 .and. wind_spinup == 1) then
   	spinup_period = wind_spinup_period * 86400.0
   	if(elapsed_time < spinup_period) then
      	spinup_function_wind = tanh(2.0*elapsed_time/spinup_period)
      else
      	spinup_function_wind = 1.0_dp
      end if
   else
      spinup_function_wind = 1.0_dp
   end if
end subroutine setup_spinup
