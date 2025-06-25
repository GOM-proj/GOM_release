!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine heat_term_by_term
	use mod_global_variables
	implicit none
	
	integer :: i,j,k,t
	integer :: t_layer, b_layer
	! real(dp):: cp = 4186.0 ! specific_heat of water [J/KgC]
	
	integer :: iday, hr
	real(dp):: j_day, solar_declination, tau_d, local_hour_angle, EQT, A0
	real(dp):: phi_sn, phi_ac, phi_an
	real(dp),allocatable,dimension(:) :: phi_br, phi_e, phi_c
	real(dp):: Ta, Tdew, w_speed, sol_rad, cloudness
	real(dp):: sol_in, sol_out
	real(dp):: u1,u2,u3,v1,v2,v3
	real(dp):: z0, alpha, bz
	real(dp):: Tw, Ts, es, ea, fWz
	! End of local variables ==================================================!
	! =========================================================================!
	! In this subroutine, I will calculate:
	! 		at global domain:
	! 			phi_sn: short-wave solar radiation
	! 			phi_an: long-wave atmospheric radiation
	! 		at each element:
	! 			phi_br: long-wave back radiation
	! 			phi_e: evaporative heat loss
	! 			phi_c: conductive heat loss
	! =========================================================================!
	! allocation & initialization
	allocate(phi_br(maxele), phi_e(maxele), phi_c(maxele))
	phi_br = 0.0_dp
	phi_e = 0.0_dp
	phi_c = 0.0_dp
	
	
	
	! Before start, find interpolated air data with Linear interpolation ======!
	u2 = julian_day*86400.0 ! julian time [s], current time and/or time for looking
	
	! initialize interpolated values, otherwise it will use garbage values.
	v2 = 0.0_dp ! not used here, but let's just keep this
	Ta = 0.0_dp
	Tdew = 0.0_dp
	w_speed = 0.0_dp
	cloudness = 0.0_dp
	
	do t=2,max_air_data_num
		! Since air_ser_time starts from the data_start_year but u2 is referenced from the simulation start year,
		! u1, and u3 should be shifted toward the simulation start year.
		! That is why I put "- reference_diff_days * 86400.0"				
		u1 = air_ser_time(t-1) - reference_diff_days * 86400.0	! [s], lower bound time
		u3 = air_ser_time(t  ) - reference_diff_days * 86400.0	! [s], upper bound time
		
		if(u2 >= u1 .and. u2 <= u3) then
			! Air temperature, T_air
			v1 = T_air(t-1)		! lower bound
			v3 = T_air(t  )		! upper bound
			Ta = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated value
			
			! Dew point temperature, T_dew
			v1 = T_dew(t-1)		! lower bound
			v3 = T_dew(t  )		! upper bound
			Tdew = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated value
			
			! Wind speed
			v1 = air_wind_speed(t-1) ! lower bound
			v3 = air_wind_speed(t  ) ! upper bound
			w_speed = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated value
			
			! Solar radiation
			v1 = solar_radiation(t-1) ! lower bound
			v3 = solar_radiation(t  ) ! upper bound
			sol_rad = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated value
			
			! Cloudness
			v1 = cloud(t-1) ! lower bound
			v3 = cloud(t  ) ! upper bound
			cloudness = (u2-u1)*(v3-v1)/(u3-u1)+v1	! interpolated value
			
			exit
		end if
	end do
	! write(*,*) 'it, Ta, Tdew, w_speed, sol_rad, cloudness = ', it, Ta, Tdew, w_speed, sol_rad, cloudness
	! end of data preparation for this timestep ===============================!

	
	! calculate short-wave solar radiations: phi_s, phi_sr, phi_sn ============!
	if(sol_swr == 1) then
		! calculate standard meridan:
		! it is calculated in "read_main_inp.f90" since it is only one-time calculation for entire model domain
		! standard_meridian = 15.0_dp * INT(lon/15.0_dp)
		
		! angular fraction of the year in [rad]
		! Here, I used 'julian_day+1' since "julian_day" in GOM starts from 0 not from 1. 
		! i.e., actual julian day in each year is 'julian_day + 1' in GOM
		j_day = julian_day + 1 ! note, I use 'jul_day' since both 'julian_day' and 'jday' are global variables
		tau_d = 2.0_dp*pi*(INT(j_day+1)-1)/365.0_dp
		
		! calculate local hour:
		! note, I can simply use 'hour', which is the global variable, instead calculating it again 
		! hr = (julian_day - INT(julian_day))*24.0_dp ! this is the original equation
		hr = hour ! from global variable
		
		! calculate solar declination angle in [radian]
		solar_declination =    0.006918 &
		&							- 0.399912*cos(tau_d)        + 0.070257*sin(tau_d) &
		&							- 0.006758*cos(2.0_dp*tau_d) + 0.000907*sin(2.0_dp*tau_d) &
		&							- 0.002697*cos(3.0_dp*tau_d) + 0.001480*sin(3.0_dp*tau_d)
		
		! calculate equation of time (EQT) in [hr]:
		EQT = 0.170*sin(4.0*pi*(INT(j_day)-80)/373.0_dp) - 0.129*sin(2.0*pi*(INT(j_day)-8)/355.0_dp)
		
		! calculate local hour angle in [radians]:
		local_hour_angle = 0.261799*(hr-(lon - standard_meridian)*0.0666667 + EQT - 12.0) ! here 0.0261799 = 2pi/24; 0.0666667 = 24/360
		
		! calculate solar altitude in [degree]:
		A0 = ASIN(sin(lat*deg2rad)*sin(solar_declination) &
		&	 + cos(lat*deg2rad)*cos(solar_declination)*cos(local_hour_angle))*rad2deg
		
		! net short-wave solar radiation, i.e., (phi_s - phi_sr):
		! note: we should include phi_sn only if A0 > 0.0
		if(A0 > 0.0) then
			phi_sn = (1.0-0.65*cloudness**2) * 24.0*(2.0444*A0 + 0.1296*A0**2 - 1.941e-3*A0**3 + 7.591e-6*A0**4)*0.1314
		end if
	else if(sol_swr == 2) then
		! use measured value in air_ser.inp
		phi_sn = (1.0-0.65*cloudness**2) * sol_rad
	end if
	
	! calculate long-wave atmospheric radiations: phi_a, phi_ar, phi_an =======!
	! calculate phi_ac, long-wave atmospheric radiation with clear sky, depending on air temperature:
	if(Ta >= 5.0) then
		phi_ac = 5.31e-13*(Ta + 273.15_dp)**6
	else
		! here 5.67e-8 is the Stefan-Boltzmann constant
		phi_ac = 5.67e-8*(Ta + 273.15_dp)**4 * (1.0-0.261*exp(-7.77e-4*Ta**2))
	end if
	
	phi_an = phi_ac*(1.0+0.17*cloudness**2)*0.97_dp
	
	! calculate rest terms at each element: phi_br & phi_e & phi_c ============!
	! Note, from here each term should be calculated at each element
	
	! wind function f(Wz):
	! note, z0 is the wind roughness height
	if(w_speed  <= 2.3) then
		z0 = 0.001
	else
		z0 = 0.005
	end if
	alpha = log(wind_height/z0)/log(2.0_dp/z0) ! conversion factor; note wind_height is given in main.inp
	bz = fWz_b/alpha**fWz_c ! fWz_b & fWz_c are given in main.inp
	fWz = fWz_a + bz*w_speed**fWz_c ! fWz is the evaporative wind speed function
	
	! actual atmospheric vapor pressure, ea:
	! note: this is for entire model domain, i.e., I don't need to calculate this cell by cell.
	ea = exp(2.3026_dp*(7.5*Tdew/(Tdew+237.3_dp) + 0.6609)) 
	
	! the rest terms will be calculated at each cell ==========================!
	!$omp parallel
	!$omp do private (i,k,t_layer,b_layer,Tw,es,ea,sol_in,sol_out)
	do i=1,maxele
		t_layer = top_layer_at_element(i)
		b_layer = bottom_layer_at_element(i)
		
		if(t_layer == 0) then
			do k=1,maxlayer
				phi_sz(k,i) = 0.0_dp
			end do
			! Note: this "cycle" statement will hurt the code if "omp" is used.
			! that is why I used "if ~ else" instead "if ~ cycle"
			! cycle ! skip rest parts of the main do-loop (i.e., go to the next element)			
		else			
			Tw = temp_cell(t_layer,i) ! water surface temperature
			
			! calculate long-wave back radiation, phi_br:
			phi_br(i) = 0.97*5.67e-8*(Tw + 273.15)**4
			
			! calculate evaporative heat loss (latent heat), phi_e:		
			! saturated vapor pressure at the water surface, es:
			if(Tw < 0.0_dp) then
				es = exp(2.3026_dp*(9.5*Tw/(Tw+265.5_dp) + 0.6609))
			else
				es = exp(2.3026_dp*(7.5*Tw/(Tw+237.3_dp) + 0.6609))
			end if
			phi_e(i) = fWz*(es - ea)
			
			! calculate conductive heat loss (sensible heat), phi_c:
			phi_c(i) = 0.47_dp*fWz*(Tw-Ta) ! 0.47 = Bowen's coefficient
			
			! calculate net surface heat flux at the surface -----------------------!
			phi_n(i) = phi_sn + phi_an - phi_br(i) - phi_e(i) - phi_c(i)
			
			! Sediment Heat Exchange, phi_sw: --------------------------------------!
			Tw = temp_cell(b_layer,i) ! water temperature at bottom layer
			
			! think again this part ....	
	! 		if(t_layer == b_layer) then
	! 			sol_sed = phi_sz(i) * sed_temp_coeff 
	! 		end if
			
	! 		sol_in = sol_out
	! 		do k=b_layer,t_layer+1
	! 			sol_out = sol_in*exp(-light_extinction*h_cell(i))
	! 			sol_net = sol_in - sol_out
	! 			sol_sed = sol_out*sed_temp_coeff
	! 			?? = sol_net+sol_sed
	! 		end do
			
			! this sediment solar back radiation will be added at all water columns (i.e., no attenunation)
			phi_sw(i) = sed_water_exchange*(T_sed - Tw) ! add sediment solar back radiation
	
			
			! calculate short-wave radiation decay (i.e., attenuation) =============#
			! The short-wave radiation that reaches to the top skin (not layer) is phi_sn.
			! Then, certain amount of phi_sn will be absorbed at the top skin (and this is included in phi_n), then the rest of the amount will be penerate.
			! So, (1-sol_absorb)*phi_sn will be started to penerate to the depth.
			! So, I need to calculate the short-wave solar radiation reaches to the bottom level of the surface layer,
			! then, this amount should be extracted from the total net solar radiation at the surface layer since phi_n already has the influx at the top skin.
			! then, at the lower layers, I need to include (top layer's bottom level amount - bottom level amount of the current layer)
			! i.e., the following is identical expression:
			! sol_in = (1.0-sol_absorb)*phi_sn
			! sol_out = sol_in*exp(-light_extinction*dz_cell(t_layer,i))
			sol_out = (1.0-sol_absorb)*phi_sn*exp(-light_extinction*dz_cell(t_layer,i))
			! phi_sz(t_layer,i) = phi_sz(t_layer,i) - sol_out
			phi_sz(t_layer,i) = phi_n(i) - sol_out ! check this again...
			do k=(t_layer-1),b_layer,-1
				sol_in = sol_out ! top level value is already calculated at the upper layer
				sol_out = sol_in*exp(-light_extinction*dz_cell(k,i)) ! top level value will be attenuated
				phi_sz(k,i) = sol_in - sol_out
			end do
			
			! at all layer, add solar back radiation & sediment solar back radiation
			! check this again if adding thiese two only in the bottom layer is correct...
			! i.e., either one of the following approaches
			do k=b_layer,t_layer
			 	phi_sz(k,i) = phi_sz(k,i) + sol_out*sed_temp_coeff		
			 	phi_sz(k,i) = phi_sz(k,i) + phi_sw(i)
			end do
			! phi_sz(b_layer,i) = phi_sz(b_layer,i) + sol_out*sed_temp_coeff		
			! phi_sz(b_layer,i) = phi_sz(b_layer,i) + phi_sw(i)
			
		end if
	end do
	!$omp end do
	!$omp end parallel
	! end of term-by-term calculations ========================================!
	
	
	deallocate(phi_br, phi_e, phi_c)
	
end subroutine heat_term_by_term