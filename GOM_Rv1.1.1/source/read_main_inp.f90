!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine reads main.inp, cell.inp, and node.inp
!! 
subroutine read_main_inp
	use mod_global_variables
	use mod_file_definition
	
	implicit none

	integer :: i, j, k
	integer :: ii
	character(len=200) :: line
	character(len= 10) :: card_num
	integer :: EOF
	
   ! jw
   integer :: index
   ! jw
   integer :: serial_num
   character(len=15) :: ctemp1
   real(dp):: rtemp1, rtemp2
   integer :: ob_eta_type_count0, ob_eta_type_count1, ob_eta_type_count2
   integer :: harmonic_ser_count, eta_ser_count, salt_ser_count, temp_ser_count
	! jw

   ! jw
	write(pw_run_log,*) "	Read main.inp"
	write(pw_run_log,*) "		Now, you are in 'read_input.f90 -> subroutine read_main_inp"
	
	! jw
	open(pw_main_inp, file = id_main_inp, form='formatted', status = 'old')
	
	! jw
	! jw
	open(pw_main_mirr, file = id_main_mirr, form='formatted', status = 'replace')
	
	card_num = "00"
	
	do 
		read(pw_main_inp,'(A)',iostat=EOF) line  ! jw
		
		! jw
		! jw
		! jw
		! jw
		! jw
		if(EOF > 0) then
			write(pw_main_mirr,*) 'Error during read main.inp'
			stop
		else if(EOF < 0) then
			write(pw_main_mirr,*) ! jw
			write(pw_main_mirr,*) 'End of file reached.'
			write(pw_main_mirr,*) 'Succeed to read main.inp'
			exit
		end if
		
		! jw
		! jw
		! jw
		! jw
		! jw
		! jw
		! jw
		if(index(line,"!") == 1 .or. index(line,"C") == 1 .or. index(line,"c") == 1) then	! jw
			write(pw_main_mirr,*) line
			cycle
		else
			backspace(unit=pw_main_inp)	! jw
		end if
		
		! jw
		select case(card_num)
!=============================================================================!
! jw
!=============================================================================!
			case("00")	! jw
				read(pw_main_inp,'(A)') project_name
				write(pw_main_mirr,'(A)') project_name
				card_num = "01"
			case("01")	! jw
				read(pw_main_inp,*) gravity, rho_o, rho_a, max_no_neighbor_node
				write(pw_main_mirr,*) gravity, rho_o, rho_a, max_no_neighbor_node
				card_num = "1"
			case("1")	! jw
				read(pw_main_inp,*) maxnod,maxele,h_node_adjust,ana_depth,coordinate_system,lon,lat,voronoi, node_mirr, cell_mirr
				write(pw_main_mirr,*) maxnod,maxele,h_node_adjust,ana_depth,coordinate_system,lon,lat,voronoi, node_mirr, cell_mirr
				
				! jw
				lon_mid = lon
				lat_mid = lat
				
				! jw
				standard_meridian = 15.0_dp * INT(lon/15.0_dp) ! jw
	
				card_num = "2"
				
			case("2")	! jw
				read(pw_main_inp,*) maxlayer, MSL, z_level(0)
				write(pw_main_mirr,*) maxlayer, MSL, z_level(0)
				card_num = "2_1"
			case("2_1")
			   ! jw
			   ! jw
			   
			   ! jw
			   delta_z_min = 1.0d15

				do k=maxlayer,1,-1
					read(pw_main_inp,*) serial_num, z_level(k)
					write(pw_main_mirr,*) serial_num, z_level(k)

			      ! jw
			      if(dabs(z_level(k) - MSL) < 1.0d-5) then
			      	write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
			         write(pw_run_log,'(A)') 'MSL is too close to level line at z layer number : ', k
			         write(*,*) 			  		'MSL is too close to level line at z layer number : ', k
			         stop
			      end if
				end do
				
				! jw
				do k=1,maxlayer
					delta_z(k) = z_level(k) - z_level(k-1)
					
					! jw
					if(delta_z(k) < delta_z_min) then
			      	delta_z_min = delta_z(k)
			      end if
				end do
				
				! jw
				
				! jw
			   do k = 1, maxlayer
			      if(dabs(z_level(k-1) + delta_z(k) - z_level(k)) > 1.e-6) then
			      	write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
			         write(pw_run_log,*) 'wrong input from vertical grid file at level', k, z_level(k-1), delta_z(k)
			         write(*,*)			  'wrong input from vertical grid file at level', k, z_level(k-1), delta_z(k)
			         stop
			      end if
			   end do
			   if(MSL > z_level(maxlayer)) then
			   	write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
			      write(pw_run_log,'(A)') 'Mean Sea Level(MSL) is above the maximum vertical level, i.e., MSL > z_level(maxlayer)'
			      write(*,*)              'Mean Sea Level(MSL) is above the maximum vertical level, i.e., MSL > z_level(maxlayer)'
			      stop
			   end if
			   if(MSL < z_level(0)) then
			   	write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
			      write(pw_run_log,'(A)') 'Mean Sea Level(MSL) is below the minimum vertical level, i.e., MSL < z_level(0)'
			      write(*,*)              'Mean Sea Level(MSL) is below the minimum vertical level, i.e., MSL < z_level(0)'
			      stop
			   end if
				
				card_num = "3"
			case("3")	! jw
				read(pw_main_inp,*) dt, ndt, data_start_year, data_time_shift, &
				&	start_year, start_month, start_day, start_hour, start_minute, start_second
				write(pw_main_mirr,*) dt, ndt, data_start_year, data_time_shift, &
				&	start_year, start_month, start_day, start_hour, start_minute, start_second
				
				! jw
				call calculate_reference_diff_days ! jw
				
				! jw
				! jw
				call julian(start_year, start_month, start_day, start_hour, start_minute, start_second, jday, jday_1900) ! jw
				
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				
				! jw
				year = start_year
				month = start_month
				day = start_day
				hour = start_hour
				minute = start_minute
				
				card_num = "4"	
			case("4")	! jw
				read(pw_main_inp,*) restart_in, restart_out, restart_freq, ishow, ishow_frequency
				write(pw_main_mirr,*) restart_in, restart_out, restart_freq, ishow, ishow_frequency
				
				card_num = "5_01"
			case("5_01")	! jw
				read(pw_main_inp,*) theta, advection_flag, adv_onoff_depth, &
				&							ELM_backtrace_flag, ELM_sub_iter, ELM_min_iter, ELM_max_iter
				write(pw_main_mirr,*) theta, advection_flag, adv_onoff_depth, &
				&							ELM_backtrace_flag, ELM_sub_iter, ELM_min_iter, ELM_max_iter

				! jw
				! jw
			   if(ELM_backtrace_flag == 0) then
			      do i=1,maxnod
			         do k=1,maxlayer
			            num_sub_elm_iteration(k,i) = ELM_sub_iter
			         end do
			      end do
			   else if(ELM_backtrace_flag == 1) then
			   	! jw
			      do i=1,maxnod
			         do k=1,maxlayer
			            num_sub_elm_iteration(k,i) = 1
			         end do
			      end do
			   else
			   	! jw
			   	do i=1,maxnod
			   		do k=1,maxlayer
			   			num_sub_elm_iteration(k,i) = 0
			   		end do
			   	end do
			   end if
				
				card_num = "5_02"
			case("5_02") ! jw
				read(pw_main_inp,*) bf_flag, bf_varying, rtemp1, rtemp2, von_Karman
				write(pw_main_mirr,*) bf_flag, bf_varying, rtemp1, rtemp2, von_Karman
				
				if(bf_flag == 1) then ! jw
					if(bf_varying == 0) then
						manning = rtemp1 ! jw
					else if(bf_varying == 1) then
						call read_bottom_roughness_inp ! jw
					end if
				else if(bf_flag == 2) then ! jw
					if(bf_varying == 0) then
						bf_height = rtemp2 ! jw
					else if(bf_varying == 1) then
						call read_bottom_roughness_inp ! jw
					end if
				end if
				
				! jw
				if(bf_flag < 0 .or. bf_flag > 2) then
					write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
					write(pw_run_log,'(A)') 'bf_flag is not in range'
					write(*,*) 			  		'bf_flag is not in range'
					stop
				end if
				if(bf_varying < 0 .or. bf_varying > 1) then
					write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
					write(pw_run_log,'(A)') 'bf_varying is not in range'
					write(*,*) 			  		'bf_varying is not in range'
					stop
				end if
								
				card_num = "5_03"
			case("5_03") ! jw
				read(pw_main_inp,*) baroclinic_flag, ana_density, ref_salt, ref_temp, dry_depth, Coriolis_option
				write(pw_main_mirr,*) baroclinic_flag, ana_density, ref_salt, ref_temp, dry_depth, Coriolis_option
				
			   denominator_min_for_matrix = min(dry_depth*(0.5-1.0e-4),delta_z_min*(0.5-1.e-4)) 

				card_num = "5_04"
			case("5_04")	! jw
				read(pw_main_inp,*) smagorinsky_parameter, Ah_0, Kh_0, Av_0, Kv_0
				write(pw_main_mirr,*) smagorinsky_parameter, Ah_0, Kh_0, Av_0, Kv_0
				
				! jw
				do j=1,maxface
					do k=1,maxlayer
						Av(k,j) = Av_0
					end do
				end do
				
				do i=1,maxele
					do k=1,maxlayer
						Kv(k,i) = Kv_0
					end do
				end do
									
			   card_num = "6"
			case("6") ! jw
				read(pw_main_inp,*) max_iteration_pcg, error_tolerance, pcg_result_show
				write(pw_main_mirr,*) max_iteration_pcg, error_tolerance, pcg_result_show
				card_num = "7"
			case("7")	! jw
				read(pw_main_inp,*) transport_flag, transport_solver
				write(pw_main_mirr,*) transport_flag, transport_solver

				if(transport_flag == 0) then
					card_num = "8"
				else if(transport_flag == 1 .and. transport_solver == 1) then
					card_num = "7_1"
				else if(transport_flag == 1 .and. transport_solver == 2) then
					card_num = "7_2"
				end if
			case("7_1")	! jw
				read(pw_main_inp,*) trans_sub_iter, h_flux_limiter, v_flux_limiter
				write(pw_main_mirr,*) trans_sub_iter, h_flux_limiter, v_flux_limiter
				
				card_num = "7_2"
			case("7_2") ! jw
				ii = 0
				do i=1,maxtran ! jw
					read(pw_main_inp,*) is_tran(i), num_tran_ser(i), tran_ser_shape(i)
					write(pw_main_mirr,*) is_tran(i), num_tran_ser(i), tran_ser_shape(i)
					
					if(is_tran(i) == 1) then
						ii = ii + 1
						tran_id(ii) = i ! jw
					end if					
				end do

				! jw
				! jw
				! jw
				if(baroclinic_flag == 1) then
					if(is_tran(1) == 0) then
						salt_cell = ref_salt
						salt_cell_new = ref_salt
						salt_node = ref_salt
						salt_face = ref_salt
					end if
					if(is_tran(2) == 0) then
						temp_cell = ref_temp
						temp_cell_new = ref_temp
						temp_node = ref_temp
						temp_face = ref_temp
					end if					
				end if
				
				
				! jw
				is_salt = is_tran(1)
				is_temp = is_tran(2)
				num_salt_ser = num_tran_ser(1)
				num_temp_ser = num_tran_ser(2)
				salt_ser_shape = tran_ser_shape(1)
				temp_ser_shape = tran_ser_shape(2)
				
				! jw
				if(is_salt == 1 .and. restart_in == 0) then ! jw
					call read_salt_init
				end if
				if(is_temp == 1 .and. restart_in == 0) then ! jw
					call read_temp_init
				end if
					
				if(is_salt == 1 .and. num_salt_ser > 0) then
					call read_salt_ser
				end if
				if(is_temp == 1 .and. num_temp_ser > 0) then
					call read_temp_ser
				end if
				if(is_temp == 1) then
					call read_air_ser
				end if
				
				if(is_tran(2) == 1) then
					card_num = "7_3"
				else
					card_num = "8"	
				end if
				
			case("7_3")
				read(pw_main_inp,*) heat_option, sol_swr, fWz_a, fWz_b, fWz_c, wind_height
				write(pw_main_mirr,*) heat_option, sol_swr, fWz_a, fWz_b, fWz_c, wind_height
				
				card_num = "7_4"
			
			case("7_4")
				read(pw_main_inp,*) light_extinction, sol_absorb, sed_water_exchange, T_sed, sed_temp_coeff
				write(pw_main_mirr,*) light_extinction, sol_absorb, sed_water_exchange, T_sed, sed_temp_coeff
				
				card_num = "8"
			
			case("8") ! jw
				read(pw_main_inp,*) &
				&	terminate_check_freq, eta_min_terminate, eta_max_terminate, uv_terminate, salt_terminate, temp_terminate
				write(pw_main_mirr,*) &
				& 	terminate_check_freq, eta_min_terminate, eta_max_terminate, uv_terminate, salt_terminate, temp_terminate
				
				card_num = "20"
				
!=============================================================================!
! jw
!=============================================================================!
			case("20")
				read(pw_main_inp,*) num_ob_cell, ob_info_option
				read(pw_main_inp,*) num_tide_bc, num_harmonic_ser, check_tide
				read(pw_main_inp,*) num_eta_bc, num_eta_ser, check_eta, eta_ser_shape
				write(pw_main_mirr,*) num_ob_cell, ob_info_option
				write(pw_main_mirr,*) num_tide_bc, num_harmonic_ser, check_tide
				write(pw_main_mirr,*) num_eta_bc, num_eta_ser, check_eta, eta_ser_shape
				
				if(num_ob_cell /= (num_tide_bc + num_eta_bc)) then
					write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
					write(pw_run_log,'(A)') 'num_tide_bc is not equal to (num_tide_bc + num_eta_bc)'
					stop 			  				'num_tide_bc is not equal to (num_tide_bc + num_eta_bc)'
				end if
				if(num_tide_bc == 0) then
					if(num_harmonic_ser /= 0) then
						write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
						write(pw_run_log,'(A)') 'If [num_tide_bc] == 0, [num_harmonic_ser] should be 0.'
						write(*,*)					'If [num_tide_bc] == 0, [num_harmonic_ser] should be 0.'
						stop
					end if
				end if
				if(num_eta_bc == 0) then
					if(num_eta_ser /= 0) then
						write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
						write(pw_run_log,'(A)') 'If [num_eta_bc] == 0, [num_eta_ser] should be 0.'
						write(*,*)					'If [num_eta_bc] == 0, [num_eta_ser] should be 0.'
						stop
					end if					
				end if
				
				
				if(num_ob_cell == 0) then
					! jw
					! jw
					! jw
					! jw
					call set_geometry_2					
					card_num = "21"
				else if(num_ob_cell > 0) then
					! jw
					! jw
					if(num_tide_bc > 0) then
						if(num_harmonic_ser > 0) then
							call read_harmonic_ser
						else
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') 'If [num_tide_bc] > 0, [num_harmonic_ser] should be greater than 0.'
							write(*,*)					'If [num_tide_bc] > 0, [num_harmonic_ser] should be greater than 0.'
							stop
						end if
					end if	
					
					! jw
					! jw
					if(num_eta_bc > 0) then
						if(num_eta_ser > 0) then
							if(eta_ser_shape == 1 .or. eta_ser_shape == 2) then
								call read_eta_ser
							else
								write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
								write(pw_run_log,'(A)') '[eta_ser_shape] should be either 1 or 2'
								write(*,*) 					'[eta_ser_shape] should be either 1 or 2'
								stop
							end if
						else
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') 'If [num_eta_bc] > 0, [num_eta_ser] should be greater than 0.'
							write(*,*) 					'If [num_eta_bc] > 0, [num_eta_ser] should be greater than 0.'
							stop
						end if
						
						! jw
						! jw
						! jw
						! jw
					end if
					
					if(ob_info_option == 0) then
						card_num = "20_1"
					else if(ob_info_option == 1) then
						! jw
						! jw
						open(pw_Obc_info,file=id_Obc_info,form='formatted',status='old')

						! jw
						call skip_header_lines(pw_Obc_info,id_Obc_info)
						
						! jw
						ob_eta_type_count0 = 0 ! jw
						ob_eta_type_count1 = 0 ! jw
						ob_eta_type_count2 = 0 ! jw
						harmonic_ser_count = 0
						eta_ser_count = 0
						salt_ser_count = 0
						temp_ser_count = 0
						
						do i=1,num_ob_cell
							read(pw_Obc_info,*) serial_num, ob_nodes(i,1), ob_nodes(i,2), ob_eta_type(i), harmonic_ser_id(i), &
							&	eta_ser_id(i), salt_ser_id(i), temp_ser_id(i) ! jw
							write(pw_main_mirr,*) serial_num, ob_nodes(i,1), ob_nodes(i,2), ob_eta_type(i), harmonic_ser_id(i), &
							&	eta_ser_id(i), salt_ser_id(i), temp_ser_id(i) ! jw
							
							! jw
							if(ob_eta_type(i) == -1) then
								ob_eta_type_count0 = ob_eta_type_count0 + 1 ! jw
							else if(ob_eta_type(i) == 1) then
								ob_eta_type_count1 = ob_eta_type_count1 + 1 ! jw
							else if(ob_eta_type(i) == 2) then
								ob_eta_type_count2 = ob_eta_type_count2 + 1 ! jw
							else
								write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
								write(pw_run_log,'(A)') '[ob_eta_type] should be one of: -1, 1, 2'
								write(*,*) 					'[ob_eta_type] should be one of: -1, 1, 2'
								stop
							end if
							
							! jw
							if(harmonic_ser_id(i) > 0) then
								harmonic_ser_count = harmonic_ser_count + 1
							end if
							if(eta_ser_id(i) > 0) then
								eta_ser_count = eta_ser_count + 1
							end if
							if(salt_ser_id(i) > 0) then
								salt_ser_count = salt_ser_count + 1
							end if
							if(temp_ser_id(i) > 0) then
								temp_ser_count = temp_ser_count + 1
							end if
							
							
							! jw
							if(ob_eta_type(i) == -1) then
								if(harmonic_ser_id(i) > 0) then
									write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
									write(pw_run_log,'(A)') '[ob_eta_type == -1], but [harmonic_ser_id > 0]'
									write(*,*) 					'[ob_eta_type == -1], but [harmonic_ser_id > 0]'
									stop
								end if
								if(eta_ser_id(i) > 0) then
									write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
									write(pw_run_log,'(A)') '[ob_eta_type == -1], but [eta_ser_id > 0]'
									write(*,*) 					'[ob_eta_type == -1], but [eta_ser_id > 0]'
									stop
								end if
								if(salt_ser_id(i) > 0) then
									write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
									write(pw_run_log,'(A)') '[ob_eta_type == -1], but [salt_ser_id > 0]'
									write(*,*) 					'[ob_eta_type == -1], but [salt_ser_id > 0]'
									stop
								end if
								if(temp_ser_id(i) > 0) then
									write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
									write(pw_run_log,'(A)') '[ob_eta_type == -1], but [temp_ser_id > 0]'
									write(*,*) 					'[ob_eta_type == -1], but [temp_ser_id > 0]'
									stop
								end if											
							end if
							
							! jw
							if(ob_eta_type(i) == 1) then
								if(harmonic_ser_id(i) < 1 .or. harmonic_ser_id(i) > num_harmonic_ser) then
									write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
									write(pw_run_log,'(A)') &
									& '[ob_eta_type == 1], but [harmonic_ser_id < 1] or [harmonic_ser_id > num_harmonic_ser]'
									write(*,*) &
									& '[ob_eta_type == 1], but [harmonic_ser_id < 1] or [harmonic_ser_id > num_harmonic_ser]'
									stop
								end if
							end if
							
							! jw
							if(ob_eta_type(i) == 2) then
								if(eta_ser_id(i) < 1 .or. eta_ser_id(i) > num_eta_ser) then
									write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
									write(pw_run_log,'(A)') '[ob_eta_type == 2], but [eta_ser_id < 1] or [eta_ser_id > num_eta_ser]'
									write(*,*) 					'[ob_eta_type == 2], but [eta_ser_id < 1] or [eta_ser_id > num_eta_ser]'
									stop
								end if
							end if
						end do ! jw
						
						
						
						! jw
						if(num_tide_bc == 0) then
							if(harmonic_ser_count > 0) then
								write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
								write(pw_run_log,'(A)') '[num_tide_bc == 0], but sum(harmonic_ser_id) > 0'
								write(*,*) 					'[num_tide_bc == 0], but sum(harmonic_ser_id) > 0'
								stop
							end if
						end if
						if(num_eta_bc == 0) then
							if(eta_ser_count > 0) then
								write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
								write(pw_run_log,'(A)') '[num_eta_bc == 0], but sum(eta_ser_id) > 0'
								write(*,*) 					'[num_eta_bc == 0], but sum(eta_ser_id) > 0'
								stop
							end if
						end if
						
						if(transport_flag == 0) then
							if(salt_ser_count > 0 .or. temp_ser_count > 0) then
								write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
								write(pw_run_log,'(A)') '[transport_flag == 0], but sum(salt_ser_id) > 0 or sum(temp_ser_id) > 0'
								write(*,*) 					'[transport_flag == 0], but sum(salt_ser_id) > 0 or sum(temp_ser_id) > 0'
								stop
							end if
						end if
						if(transport_flag > 0) then
							if(is_tran(1) > 0) then ! jw
								if(maxval(salt_ser_id) > num_tran_ser(1)) then
									write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
									write(pw_run_log,'(A)') &
									& '[transport_flag == 1 & is_tran(1) > 0], then maxval(salt_ser_id) > num_tran_ser(1)'
									write(*,*) &
									& '[transport_flag == 1 & is_tran(1) > 0], then maxval(salt_ser_id) > num_tran_ser(1)'
									stop
								end if
							end if
							if(is_tran(2) > 0) then ! jw
								if(maxval(temp_ser_id) > num_tran_ser(2)) then
									write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
									write(pw_run_log,'(A)') &
									& '[transport_flag == 1 & is_tran(2) > 0], then maxval(temp_ser_id) > num_tran_ser(2)'
									write(*,*) &
									& '[transport_flag == 1 & is_tran(2) > 0], then maxval(temp_ser_id) > num_tran_ser(2)'
									stop
								end if
							end if
						end if
						
						! jw
						close(pw_Obc_info)
						call set_geometry_2
						
						! jw
						card_num = "21"
					end if
				end if
			case("20_1")
				ob_eta_type_count0 = 0 ! jw
				ob_eta_type_count1 = 0 ! jw
				ob_eta_type_count2 = 0 ! jw
				harmonic_ser_count = 0
				eta_ser_count = 0
				salt_ser_count = 0
				temp_ser_count = 0
				
				do i=1,num_ob_cell
					read(pw_main_inp,*) serial_num, ob_nodes(i,1), ob_nodes(i,2), ob_eta_type(i), harmonic_ser_id(i), &
					&	eta_ser_id(i), salt_ser_id(i), temp_ser_id(i) ! jw
					write(pw_main_mirr,*) serial_num, ob_nodes(i,1), ob_nodes(i,2), ob_eta_type(i), harmonic_ser_id(i), &
					&	eta_ser_id(i), salt_ser_id(i), temp_ser_id(i) ! jw
					
					! jw
					if(ob_eta_type(i) == -1) then
						ob_eta_type_count0 = ob_eta_type_count0 + 1 ! jw
					else if(ob_eta_type(i) == 1) then
						ob_eta_type_count1 = ob_eta_type_count1 + 1 ! jw
					else if(ob_eta_type(i) == 2) then
						ob_eta_type_count2 = ob_eta_type_count2 + 1 ! jw
					else
						write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
						write(pw_run_log,'(A)') '[ob_eta_type] should be one of: -1, 1, 2'
						write(*,*) 					'[ob_eta_type] should be one of: -1, 1, 2'
						stop
					end if
					
					! jw
					if(harmonic_ser_id(i) > 0) then
						harmonic_ser_count = harmonic_ser_count + 1
					end if
					if(eta_ser_id(i) > 0) then
						eta_ser_count = eta_ser_count + 1
					end if
					if(salt_ser_id(i) > 0) then
						salt_ser_count = salt_ser_count + 1
					end if
					if(temp_ser_id(i) > 0) then
						temp_ser_count = temp_ser_count + 1
					end if
					
					
					! jw
					if(ob_eta_type(i) == -1) then
						if(harmonic_ser_id(i) > 0) then
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') '[ob_eta_type == -1], but [harmonic_ser_id > 0]'
							write(*,*) 					'[ob_eta_type == -1], but [harmonic_ser_id > 0]'
							stop
						end if
						if(eta_ser_id(i) > 0) then
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') '[ob_eta_type == -1], but [eta_ser_id > 0]'
							write(*,*) 					'[ob_eta_type == -1], but [eta_ser_id > 0]'
							stop
						end if
						if(salt_ser_id(i) > 0) then
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') '[ob_eta_type == -1], but [salt_ser_id > 0]'
							write(*,*) 					'[ob_eta_type == -1], but [salt_ser_id > 0]'
							stop
						end if
						if(temp_ser_id(i) > 0) then
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') '[ob_eta_type == -1], but [temp_ser_id > 0]'
							write(*,*) 					'[ob_eta_type == -1], but [temp_ser_id > 0]'
							stop
						end if											
					end if
					
					! jw
					if(ob_eta_type(i) == 1) then
						if(harmonic_ser_id(i) < 1 .or. harmonic_ser_id(i) > num_harmonic_ser) then
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') &
							& '[ob_eta_type == 1], but [harmonic_ser_id < 1] or [harmonic_ser_id > num_harmonic_ser]'
							write(*,*) &
							& '[ob_eta_type == 1], but [harmonic_ser_id < 1] or [harmonic_ser_id > num_harmonic_ser]'
							stop
						end if
					end if
					
					! jw
					if(ob_eta_type(i) == 2) then
						if(eta_ser_id(i) < 1 .or. eta_ser_id(i) > num_eta_ser) then
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') '[ob_eta_type == 2], but [eta_ser_id < 1] or [eta_ser_id > num_eta_ser]'
							write(*,*) 					'[ob_eta_type == 2], but [eta_ser_id < 1] or [eta_ser_id > num_eta_ser]'
							stop
						end if
					end if
				end do ! jw
				
				
				
				! jw
				if(num_tide_bc == 0) then
					if(harmonic_ser_count > 0) then
						write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
						write(pw_run_log,'(A)') '[num_tide_bc == 0], but sum(harmonic_ser_id) > 0'
						write(*,*) 					'[num_tide_bc == 0], but sum(harmonic_ser_id) > 0'
						stop
					end if
				end if
				if(num_eta_bc == 0) then
					if(eta_ser_count > 0) then
						write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
						write(pw_run_log,'(A)') '[num_eta_bc == 0], but sum(eta_ser_id) > 0'
						write(*,*) 					'[num_eta_bc == 0], but sum(eta_ser_id) > 0'
						stop
					end if
				end if
				
				if(transport_flag == 0) then
					if(salt_ser_count > 0 .or. temp_ser_count > 0) then
						write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
						write(pw_run_log,'(A)') '[transport_flag == 0], but sum(salt_ser_id) > 0 or sum(temp_ser_id) > 0'
						write(*,*) 					'[transport_flag == 0], but sum(salt_ser_id) > 0 or sum(temp_ser_id) > 0'
						stop
					end if
				end if
				if(transport_flag > 0) then
					if(is_tran(1) > 0) then ! jw
						if(maxval(salt_ser_id) > num_tran_ser(1)) then
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') &
							& '[transport_flag == 1 & is_tran(1) > 0], then maxval(salt_ser_id) > num_tran_ser(1)'
							write(*,*) &
							& '[transport_flag == 1 & is_tran(1) > 0], then maxval(salt_ser_id) > num_tran_ser(1)'
							stop
						end if
					end if
					if(is_tran(2) > 0) then ! jw
						if(maxval(temp_ser_id) > num_tran_ser(2)) then
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') &
							& '[transport_flag == 1 & is_tran(2) > 0], then maxval(temp_ser_id) > num_tran_ser(2)'
							write(*,*) &
							& '[transport_flag == 1 & is_tran(2) > 0], then maxval(temp_ser_id) > num_tran_ser(2)'
							stop
						end if
					end if
				end if
				
				call set_geometry_2

				card_num = "21"
			case("21")	! jw
				read(pw_main_inp,*) tide_spinup, tide_spinup_period
				read(pw_main_inp,*) baroclinic_spinup, baroclinic_spinup_period
				write(pw_main_mirr,*)  tide_spinup, tide_spinup_period
				write(pw_main_mirr,*)  baroclinic_spinup, baroclinic_spinup_period
				
				card_num = "22"
			case("22")
				read(pw_main_inp,*) num_Qb_cell, Qbc_info_option, num_Q_ser, check_Q, Q_ser_shape
				write(pw_main_mirr,*) num_Qb_cell, Qbc_info_option, num_Q_ser, check_Q, Q_ser_shape
				
				if(num_Qb_cell == 0) then
					card_num = "23"
				else if(num_Qb_cell > 0) then
					if(Qbc_info_option == 0) then
						card_num = "22_1"
					else if(Qbc_info_option == 1) then
						open(pw_Qbc_info,file=id_Qbc_info,form='formatted',status='old')

						! jw
						call skip_header_lines(pw_Qbc_info,id_Qbc_info)

						do i=1,num_Qb_cell
							! jw

							read(pw_Qbc_info,*) serial_num, Q_nodes(i,1), Q_nodes(i,2), Q_ser_id(i), Q_portion(i), &
							&	Q_salt_ser_id(i), Q_temp_ser_id(i)
							! jw
							! jw

							! jw
							call find_element_and_face(Q_nodes(i,1),Q_nodes(i,2))
							! jw
												
							if(element_id == 0 .or. face_id == 0) then
								write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
								write(pw_run_log,'(A,I5)') &
								&	"Fail to find the corresponding element at river number: ", i
								write(pw_run_log,'(A,I5)') &
								&	"Check whether given node numvers are in counterclockwise direction..."
								write(*,*) "Fail to find the corresponding element at river number: ", i
								write(*,*) "Check whether given node numvers are in counterclockwise direction or not ..."
								stop
							end if
							
							! jw
							Q_boundary(i,1) = element_id 			! jw
							Q_boundary(i,2) = face_id 				! jw
							isflowside3(Q_boundary(i,2)) = i 	! jw
							Qb_element_flag(element_id) = i 		! jw
						end do
						
						! jw
						! jw
						! jw
						! jw
						if(maxval(Q_ser_id) > num_Q_ser) then
							write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
							write(pw_run_log,'(A)') &
							&	"One of Q_ser_id(i) is geater than num_Q_ser &
							&	 Q_ser_id(i) can be smaller than num_Q_ser, but it should not be bigger than num_Q_ser, i.e., &
							&	 Q_ser_id(i) <= num_Q_ser	-> correct &
							&	 Q_ser_id(i) >  num_Q_ser	-> incorrect"
							write(*,*) &
							&	"One of Q_ser_id(i) is geater than num_Q_ser &
							&   Q_ser_id(i) can be smaller than num_Q_ser, but it should not be bigger than num_Q_ser, i.e., &
							&   Q_ser_id(i) <= num_Q_ser	-> correct &
							&   Q_ser_id(i) >  num_Q_ser	-> incorrect"
							stop
						end if
						
						close(pw_Qbc_info)
						
						! jw
						call read_q_ser ! jw
						! jw

						! jw
						card_num = "23"
					end if
				end if
			case("22_1")
				do i=1,num_Qb_cell
					! jw

					read(pw_main_inp,*) serial_num, Q_nodes(i,1), Q_nodes(i,2), Q_ser_id(i), Q_portion(i), &
					&	Q_salt_ser_id(i), Q_temp_ser_id(i)
					write(pw_main_mirr,*) serial_num, Q_nodes(i,1), Q_nodes(i,2), Q_ser_id(i), Q_portion(i), &
					&	Q_salt_ser_id(i), Q_temp_ser_id(i)

					! jw
					call find_element_and_face(Q_nodes(i,1),Q_nodes(i,2))
					! jw
										
					if(element_id == 0 .or. face_id == 0) then
						write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
						write(pw_run_log,'(A,I5)') &
						&	"Fail to find the corresponding element. Stop.; read_main_inp.f90; C22_1; at river number: ", i
						write(pw_run_log,'(A,I5)') &
						&	"Check whether given node numvers are in counterclockwise direction..."
						write(*,*) "Fail to find the corresponding element. Stop.; read_main_inp.f90; C22_1; at river number: ", i
						write(*,*) "Check whether given node numbers are in counterclockwise direction or not ..."
						stop
					end if
					
					! jw
					Q_boundary(i,1) = element_id 			! jw
					Q_boundary(i,2) = face_id 				! jw
					isflowside3(Q_boundary(i,2)) = i 	! jw
					Qb_element_flag(element_id) = i 		! jw
					
					! jw
					! jw
				end do
				
				! jw
				! jw
				! jw
				! jw
				if(maxval(Q_ser_id) > num_Q_ser) then
					write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
					write(pw_run_log,'(A)') &
					&	"One of Q_ser_id(i) is geater than num_Q_ser &
					&	 Q_ser_id(i) can be smaller than num_Q_ser, but it should not be bigger than num_Q_ser, i.e., &
					&	 Q_ser_id(i) <= num_Q_ser	-> correct &
					&	 Q_ser_id(i) >  num_Q_ser	-> incorrect"
					write(*,*) &
					&	"One of Q_ser_id(i) is geater than num_Q_ser &
					&   Q_ser_id(i) can be smaller than num_Q_ser, but it should not be bigger than num_Q_ser, i.e., &
					&   Q_ser_id(i) <= num_Q_ser	-> correct &
					&   Q_ser_id(i) >  num_Q_ser	-> incorrect"
					stop
				end if
				
				! jw
				call read_q_ser ! jw
				! jw
				
				card_num = "23"
			case("23")
				read(pw_main_inp,*) num_WR_cell
				write(pw_main_mirr,*) num_WR_cell
				
				if(num_WR_cell > 0) then
					card_num = "23_1"
				else if(num_WR_cell == 0) then
					card_num = "24"
				end if
			case("23_1")
				do i=1,num_WR_cell
					read(pw_main_inp,*) serial_num, WR_nodes(i,1), WR_nodes(i,2), WR_layer(i), WR_Q_ser_id(i), WR_portion(i), &
					&	WR_salt_ser_id(i), WR_temp_ser_id(i)	! jw
					write(pw_main_mirr,*) serial_num, WR_nodes(i,1), WR_nodes(i,2), WR_layer(i), WR_Q_ser_id(i), WR_portion(i), &
					&	WR_salt_ser_id(i), WR_temp_ser_id(i)	! jw

					! jw
					call find_element_and_face(WR_nodes(i,1),WR_nodes(i,2))
										
					if(element_id == 0 .or. face_id == 0) then
						write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
						write(pw_run_log,'(A,I5)') &
						&	"Fail to find the corresponding element. Stop.; read_main_inp.f90; C23_1; at WR number: ", i
						write(pw_run_log,'(A,I5)') &
						&	"Check whether given node numvers are in counterclockwise direction..."
						write(*,*) "Fail to find the corresponding element. Stop.; read_main_inp.f90; C23_1; at WR number: ", i
						write(*,*) "Check whether given node numbers are in counterclockwise direction or not ..."
						stop
					end if

					! jw
					WR_boundary(i,1) = element_id 			! jw
					WR_boundary(i,2) = face_id 				! jw
					isflowside4(WR_boundary(i,2)) = i 	! jw
					WR_element_flag(element_id) = i 		! jw
				end do
				
				card_num = "24"
			case("24")
				read(pw_main_inp,*) num_SS_cell
				write(pw_main_mirr,*) num_SS_cell
				
				if(num_SS_cell > 0) then
					card_num = "24_1"
				else if(num_SS_cell == 0) then
					card_num = "25"
				end if
			case("24_1")
				do i=1,num_SS_cell
					read(pw_main_inp,*) serial_num, SS_cell(i), SS_layer(i), SS_Q_ser_id(i), SS_portion(i), &
					&	SS_salt_ser_id(i), SS_temp_ser_id(i)	! jw
					write(pw_main_mirr,*) serial_num, SS_cell(i), SS_layer(i), SS_Q_ser_id(i), SS_portion(i), &
					&	SS_salt_ser_id(i), SS_temp_ser_id(i)	! jw
					
					! jw
					SS_element_flag(SS_cell(i)) = i
				end do
				card_num = "25"
			case("25")	! jw
				read(pw_main_inp,*) no_Etide_species, Etide_cutoff_depth
				write(pw_main_mirr,*) no_Etide_species, Etide_cutoff_depth
				
				if(no_Etide_species > 0) then
					card_num = "25_1"
				else if(no_Etide_species == 0) then
					card_num = "30"
				end if
			case("25_1") ! jw
		      do i = 1, no_Etide_species
		         read(pw_main_inp,*) &
		         &	Etide_species(i), Etide_amplitude(i), Etide_frequency(i), Etide_nodal_factor(i), Etide_astro_arg_degree(i)
		         write(pw_main_mirr,*) &
		         &	Etide_species(i), Etide_amplitude(i), Etide_frequency(i), Etide_nodal_factor(i), Etide_astro_arg_degree(i)
		           
		         Etide_astro_arg_degree(i) = Etide_astro_arg_degree(i)*deg2rad
		      end do 			      
				
				! jw
			  	open(32,file='./input/long_lati_mesh.dat',status='old')
		      read(32,*)
		      read(32,*)
		      
		      ! jw
		      do i=1,maxnod
		         read(32,*) j, lon_node(i), lat_node(i)
		         lon_node(i) = lon_node(i) * deg2rad
		         lat_node(i) = lat_node(i)  * deg2rad
		         Etide_species_coef_at_node(0,i) = 3.0*dsin(    lat_node(i) )**2-1.0
		         Etide_species_coef_at_node(1,i) =     dsin(2.0*lat_node(i) )
		         Etide_species_coef_at_node(2,i) =     dcos(    lat_node(i) )**2
		      end do
		      close(32)
		      
			  	do i = 1, maxele
			   	lon_cell(i) = 0.0
			     	lat_cell(i)  = 0.0
			     	do j = 1, tri_or_quad(i)
			        	lon_cell(i) = lon_cell(i) + lon_node(nodenum_at_cell(j,i)) / tri_or_quad(i)
			        	lat_cell(i)  = lat_cell(i)  + lat_node(nodenum_at_cell(j,i)) / tri_or_quad(i)
			     	end do
			     	Etide_species_coef_at_element(0,i) = 3.0*dsin(    lat_cell(i) )**2-1.0
			     	Etide_species_coef_at_element(1,i) =     dsin(2.0*lat_cell(i) )
			     	Etide_species_coef_at_element(2,i) =     dcos(    lat_cell(i) )**2
			  	end do
		      ! jw
				
				card_num = "30"
				

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

!=============================================================================!
! jw
!=============================================================================!			
			case("30")	! jw
				read(pw_main_inp,*) wind_flag, airp_flag, num_windp_ser, wind_formular, wind_spinup, wind_spinup_period
				write(pw_main_mirr,*) wind_flag, airp_flag, num_windp_ser, wind_formular, wind_spinup, wind_spinup_period
				
				if(wind_flag == 1 .or. airp_flag == 1) then
					if(num_windp_ser == 0) then
						write(pw_run_log,'(A,A)') 'Error in main.inp, Card#: ', card_num
						write(pw_run_log,'(A)') '[num_windp_ser] should be greater than 0'
						write(*,*) 					'[num_windp_ser] should be greater than 0'
						stop
					else if(num_windp_ser > 0) then
						call read_windp_ser
					end if
				end if
				
				card_num = "31"
			case("31")	! jw
				read(pw_main_inp,*) hurricane_flag, hurricane_data_type, hurricane_interp_method, hurricane_dt
				write(pw_main_mirr,*) hurricane_flag, hurricane_data_type, hurricane_interp_method, hurricane_dt
				
				if(hurricane_flag == 1) then
					! jw
					call read_hurricane_ser_0
					
					if(hurricane_interp_method == 2) then
					   allocate(wind_u1_new(maxnod),wind_v1_new(maxnod))
					   allocate(wind_u2_new(maxnod),wind_v2_new(maxnod))
					   allocate(air_p1_new(maxnod),air_p2_new(maxnod))
					   allocate(shiftx(maxnod),shifty(maxnod))
					   wind_u1_new = 0.0_dp
					   wind_v1_new = 0.0_dp
					   wind_u2_new = 0.0_dp
					   wind_v2_new = 0.0_dp
					   air_p1_new 	= 0.0_dp 
					   air_p2_new 	= 0.0_dp
					   shiftx 		= 0.0_dp
					   shifty 		= 0.0_dp
					end if
				end if
				
				if(hurricane_flag == 0) then				
					card_num = "32"
				else if(hurricane_flag == 1) then
					card_num = "31_1"
				end if
			case("31_1") ! jw
				read(pw_main_inp,*) hurricane_start_year, hurricane_start_month, hurricane_start_day
				read(pw_main_inp,*) hurricane_end_year, hurricane_end_month, hurricane_end_day
				write(pw_main_mirr,*) hurricane_start_year, hurricane_start_month, hurricane_start_day
				write(pw_main_mirr,*) hurricane_end_year, hurricane_end_month, hurricane_end_day
				
				! jw
				call find_hurricane_start_end_jday
				
				card_num = "32"
			case("32")	! jw
				read(pw_main_inp,*) holland_flag, holland_start_jday, holland_end_jday
				write(pw_main_mirr,*) holland_flag, holland_start_jday, holland_end_jday
				
				if(holland_flag == 1) then
			      open(650,file='./input/fort.650')   ! jw
			      read(650,*)
			      ! jw
				end if
				
				card_num = "40"
!=============================================================================!
! jw
!=============================================================================!
      	case("40")	! jw
      		read(pw_main_inp,*) &
      		&	check_grid_2D, check_grid_2DO, check_grid_3D, check_grid_format, check_grid_info, check_grid_unit_conv
      		write(pw_main_mirr,*) &
      		&	check_grid_2D, check_grid_2DO, check_grid_3D, check_grid_format, check_grid_info, check_grid_unit_conv
      		card_num = "41"
         case("41")
				read(pw_main_inp,*) tser_station_num, tser_info_option, tser_hloc, tser_time, tser_frequency, tser_format
				write(pw_main_mirr,*) tser_station_num, tser_info_option, tser_hloc, tser_time, tser_frequency, tser_format
				
				if(tser_station_num == 0) then
					card_num = "42"
				else if(tser_station_num > 0) then
					if(tser_info_option == 0) then
						card_num = "41_1"
					else if(tser_info_option == 1) then
						open(pw_tser_station_info,file=id_tser_station_info,form='formatted',status='old')
						! jw
						call skip_header_lines(pw_tser_station_info,id_tser_station_info)
						do i=1, tser_station_num
							read(pw_tser_station_info,*) serial_num, tser_station_cell(i), tser_station_node(i), tser_station_name(i)
						end do						
						close(pw_tser_station_info)
						
						! jw
						card_num = "41_2"
					end if
				end if
			case("41_1")
				do i=1, tser_station_num
					read(pw_main_inp,*) serial_num, tser_station_cell(i), tser_station_node(i), tser_station_name(i)
					write(pw_main_mirr,*) serial_num, tser_station_cell(i), tser_station_node(i), tser_station_name(i)
				end do
				card_num = "41_2"
			case("41_2")
				read(pw_main_inp,*) tser_eta, tser_H, tser_u, tser_v, tser_salt, tser_temp, tser_airp
				write(pw_main_mirr,*) tser_eta, tser_H, tser_u, tser_v, tser_salt, tser_temp, tser_airp
				card_num = "42"
			case("42")
				read(pw_main_inp,*) IS2D_switch, IS2D_flood_map
				read(pw_main_inp,*) IS2D_format, IS2D_binary, IS2D_File_freq
				read(pw_main_inp,*) IS2D_time, IS2D_frequency, IS2D_start, IS2D_end, IS2D_unit_conv
				
				write(pw_main_mirr,*) IS2D_switch, IS2D_flood_map
				write(pw_main_mirr,*) IS2D_format, IS2D_binary, IS2D_File_freq
				write(pw_main_mirr,*) IS2D_time, IS2D_frequency, IS2D_start, IS2D_end, IS2D_unit_conv
				
				! jw
				card_num = "42_1"
			case("42_1")
				do i=1,6
					read(pw_main_inp,*) ctemp1, IS2D_variable(i)
					write(pw_main_mirr,*) ctemp1, IS2D_variable(i)	
				end do
				card_num = "43"
			case("43")	
				read(pw_main_inp,*) IS3D_full_switch, IS3D_grid_format, IS3D_surf_switch
				read(pw_main_inp,*) IS3D_format, IS3D_binary, IS3D_File_freq
				read(pw_main_inp,*) IS3D_time, IS3D_frequency, IS3D_start, IS3D_end, IS3D_unit_conv
				
				write(pw_main_mirr,*) IS3D_full_switch, IS3D_grid_format, IS3D_surf_switch
				write(pw_main_mirr,*) IS3D_format, IS3D_binary, IS3D_File_freq
				write(pw_main_mirr,*) IS3D_time, IS3D_frequency, IS3D_start, IS3D_end, IS3D_unit_conv
				
				! jw
				card_num = "43_1"
			case("43_1")
				do i=1,6
					read(pw_main_inp,*) ctemp1, IS3D_variable(i)
					write(pw_main_mirr,*) ctemp1, IS3D_variable(i)	
				end do
				card_num = "44"
			case("44")
				read(pw_main_inp,*) IS2D_dump_switch, IS2D_dump_binary, IS2D_dump_time, IS2D_dump_File_freq
				read(pw_main_inp,*) IS2D_dump_frequency, IS2D_dump_start, IS2D_dump_end
				
				write(pw_main_mirr,*) IS2D_dump_switch, IS2D_dump_binary, IS2D_dump_time, IS2D_dump_File_freq
				write(pw_main_mirr,*) IS2D_dump_frequency, IS2D_dump_start, IS2D_dump_end
				
				card_num = "45"
			case("45")
				read(pw_main_inp,*) IS3D_dump_switch, IS3D_dump_binary, IS3D_dump_time, IS3D_dump_File_freq
				read(pw_main_inp,*) IS3D_dump_frequency, IS3D_dump_start, IS3D_dump_end
				
				write(pw_main_mirr,*) IS3D_dump_switch, IS3D_dump_binary, IS3D_dump_time, IS3D_dump_File_freq
				write(pw_main_mirr,*) IS3D_dump_frequency, IS3D_dump_start, IS3D_dump_end
				
				card_num = "46"
			case("46")
				! jw
				read(pw_main_inp,*) dia_advection, dia_momentum, dia_freesurface, dia_eta_at_ob
				write(pw_main_mirr,*) dia_advection, dia_momentum, dia_freesurface, dia_eta_at_ob
				
				! jw
				read(pw_main_inp,*) dia_bottom_friction, dia_face_velocity, dia_node_velocity
				write(pw_main_mirr,*) dia_bottom_friction, dia_face_velocity, dia_node_velocity
				
				card_num = "999"
			case("999")
				! jw
				exit		
			case default
				! jw
		end select
	end do
	
  	close(pw_main_inp)		! jw
  	close(pw_main_mirr)		! jw
end subroutine read_main_inp
