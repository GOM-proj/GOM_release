!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine heat_term_by_term
	use mod_global_variables
	implicit none
	
	integer :: i,j,k,t
	integer :: t_layer, b_layer
	! jw
	
	integer :: iday, hr
	real(dp):: j_day, solar_declination, tau_d, local_hour_angle, EQT, A0
	real(dp):: phi_sn, phi_ac, phi_an
	real(dp),allocatable,dimension(:) :: phi_br, phi_e, phi_c
	real(dp):: Ta, Tdew, w_speed, sol_rad, cloudness
	real(dp):: sol_in, sol_out
	real(dp):: u1,u2,u3,v1,v2,v3
	real(dp):: z0, alpha, bz
	real(dp):: Tw, Ts, es, ea, fWz
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	allocate(phi_br(maxele), phi_e(maxele), phi_c(maxele))
	phi_br = 0.0_dp
	phi_e = 0.0_dp
	phi_c = 0.0_dp
	
	
	
	! jw
	u2 = julian_day*86400.0 ! jw
	
	! jw
	v2 = 0.0_dp ! jw
	Ta = 0.0_dp
	Tdew = 0.0_dp
	w_speed = 0.0_dp
	cloudness = 0.0_dp
	
	do t=2,max_air_data_num
		! jw
		! jw
		! jw
		u1 = air_ser_time(t-1) - reference_diff_days * 86400.0	! jw
		u3 = air_ser_time(t  ) - reference_diff_days * 86400.0	! jw
		
		if(u2 >= u1 .and. u2 <= u3) then
			! jw
			v1 = T_air(t-1)		! jw
			v3 = T_air(t  )		! jw
			Ta = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
			
			! jw
			v1 = T_dew(t-1)		! jw
			v3 = T_dew(t  )		! jw
			Tdew = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
			
			! jw
			v1 = air_wind_speed(t-1) ! jw
			v3 = air_wind_speed(t  ) ! jw
			w_speed = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
			
			! jw
			v1 = solar_radiation(t-1) ! jw
			v3 = solar_radiation(t  ) ! jw
			sol_rad = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
			
			! jw
			v1 = cloud(t-1) ! jw
			v3 = cloud(t  ) ! jw
			cloudness = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
			
			exit
		end if
	end do
	! jw
	! jw

	
	! jw
	if(sol_swr == 1) then
		! jw
		! jw
		! jw
		
		! jw
		! jw
		! jw
		j_day = julian_day + 1 ! jw
		tau_d = 2.0_dp*pi*(INT(j_day+1)-1)/365.0_dp
		
		! jw
		! jw
		! jw
		hr = hour ! jw
		
		! jw
		solar_declination =    0.006918 &
		&							- 0.399912*cos(tau_d)        + 0.070257*sin(tau_d) &
		&							- 0.006758*cos(2.0_dp*tau_d) + 0.000907*sin(2.0_dp*tau_d) &
		&							- 0.002697*cos(3.0_dp*tau_d) + 0.001480*sin(3.0_dp*tau_d)
		
		! jw
		EQT = 0.170*sin(4.0*pi*(INT(j_day)-80)/373.0_dp) - 0.129*sin(2.0*pi*(INT(j_day)-8)/355.0_dp)
		
		! jw
		local_hour_angle = 0.261799*(hr-(lon - standard_meridian)*0.0666667 + EQT - 12.0) ! jw
		
		! jw
		A0 = ASIN(sin(lat*deg2rad)*sin(solar_declination) &
		&	 + cos(lat*deg2rad)*cos(solar_declination)*cos(local_hour_angle))*rad2deg
		
		! jw
		! jw
		if(A0 > 0.0) then
			phi_sn = (1.0-0.65*cloudness**2) * 24.0*(2.0444*A0 + 0.1296*A0**2 - 1.941e-3*A0**3 + 7.591e-6*A0**4)*0.1314
		end if
	else if(sol_swr == 2) then
		! jw
		phi_sn = (1.0-0.65*cloudness**2) * sol_rad
	end if
	
	! jw
	! jw
	if(Ta >= 5.0) then
		phi_ac = 5.31e-13*(Ta + 273.15_dp)**6
	else
		! jw
		phi_ac = 5.67e-8*(Ta + 273.15_dp)**4 * (1.0-0.261*exp(-7.77e-4*Ta**2))
	end if
	
	phi_an = phi_ac*(1.0+0.17*cloudness**2)*0.97_dp
	
	! jw
	! jw
	
	! jw
	! jw
	if(w_speed  <= 2.3) then
		z0 = 0.001
	else
		z0 = 0.005
	end if
	alpha = log(wind_height/z0)/log(2.0_dp/z0) ! jw
	bz = fWz_b/alpha**fWz_c ! jw
	fWz = fWz_a + bz*w_speed**fWz_c ! jw
	
	! jw
	! jw
	ea = exp(2.3026_dp*(7.5*Tdew/(Tdew+237.3_dp) + 0.6609)) 
	
	! jw
	!$omp parallel
	!$omp do private (i,k,t_layer,b_layer,Tw,es,ea,sol_in,sol_out)
	do i=1,maxele
		t_layer = top_layer_at_element(i)
		b_layer = bottom_layer_at_element(i)
		
		if(t_layer == 0) then
			do k=1,maxlayer
				phi_sz(k,i) = 0.0_dp
			end do
			! jw
			! jw
			! jw
		else			
			Tw = temp_cell(t_layer,i) ! jw
			
			! jw
			phi_br(i) = 0.97*5.67e-8*(Tw + 273.15)**4
			
			! jw
			! jw
			if(Tw < 0.0_dp) then
				es = exp(2.3026_dp*(9.5*Tw/(Tw+265.5_dp) + 0.6609))
			else
				es = exp(2.3026_dp*(7.5*Tw/(Tw+237.3_dp) + 0.6609))
			end if
			phi_e(i) = fWz*(es - ea)
			
			! jw
			phi_c(i) = 0.47_dp*fWz*(Tw-Ta) ! jw
			
			! jw
			phi_n(i) = phi_sn + phi_an - phi_br(i) - phi_e(i) - phi_c(i)
			
			! jw
			Tw = temp_cell(b_layer,i) ! jw
			
			! jw
	! jw
	! jw
	! jw
			
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
			
			! jw
			phi_sw(i) = sed_water_exchange*(T_sed - Tw) ! jw
	
			
			! jw
			! jw
			! jw
			! jw
			! jw
			! jw
			! jw
			! jw
			! jw
			! jw
			sol_out = (1.0-sol_absorb)*phi_sn*exp(-light_extinction*dz_cell(t_layer,i))
			! jw
			phi_sz(t_layer,i) = phi_n(i) - sol_out ! jw
			do k=(t_layer-1),b_layer,-1
				sol_in = sol_out ! jw
				sol_out = sol_in*exp(-light_extinction*dz_cell(k,i)) ! jw
				phi_sz(k,i) = sol_in - sol_out
			end do
			
			! jw
			! jw
			! jw
			do k=b_layer,t_layer
			 	phi_sz(k,i) = phi_sz(k,i) + sol_out*sed_temp_coeff		
			 	phi_sz(k,i) = phi_sz(k,i) + phi_sw(i)
			end do
			! jw
			! jw
			
		end if
	end do
	!$omp end do
	!$omp end parallel
	! jw
	
	
	deallocate(phi_br, phi_e, phi_c)
	
end subroutine heat_term_by_term