!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! This subroutine will read some important variables from given input files,
!! and the varialbes will be used to allocate arrays. 
subroutine scan_gom
	use omp_lib
	use mod_global_variables
	use mod_file_definition	
	implicit none

	integer :: i, j
	integer,allocatable :: nd(:,:)
	integer :: serial_num, node_count
	real	  :: r1, r2, r3
	character(len=5) :: card_num
	! End of local variables ==================================================!

	! Scan main.inp ===========================================================!
	write(pw_run_log,*) "	Scanning main.inp"
		
	! Open 'main.inp, which should exist and should be a read-only file. 
	open(pw_main_inp, file = id_main_inp, form='formatted', status = 'old')

!=============================================================================!
! C0 ~ C19: General Model Configuration Information									!
!=============================================================================!
	! find: max_no_neighbor_node ----------------------------------------------!
	card_num = "C01"
	call seek_card(card_num)
	read(pw_main_inp,*) r1, r2, r3, max_no_neighbor_node
	write(*,'(A30,I10)') 'max_no_nieghbor_node = ', max_no_neighbor_node
	
	! find: maxnod, maxele ----------------------------------------------------!
	card_num = "C1"
	call seek_card(card_num)
	read(pw_main_inp,*) maxnod, maxele, h_node_adjust, ana_depth, coordinate_system, lon, lat, voronoi, node_mirr, cell_mirr
	write(*,'(A30,I10)') 'maxnod = ', maxnod
	write(*,'(A30,I10)') 'maxele = ', maxele
	
	! find: maxlayer ----------------------------------------------------------!
	card_num = "C2"
	call seek_card(card_num)
	read(pw_main_inp,*) maxlayer, MSL				
	write(*,'(A30,I10)') 'maxlayer = ', maxlayer
				
	! find: maxface -----------------------------------------------------------!
	! To find maxface, I need to read node.inp and cell.inp
	! Note: should read node.inp first
	call read_node_inp
	call read_cell_inp
	
	! find: max_salt_data_num, max_salt_ser, max_temp_data_num, max_temp_ser
	card_num = "C7"
	call seek_card(card_num)
	read(pw_main_inp,*) transport_flag ! don't need to read the rest part
	
	if(transport_flag == 1) then
		card_num = "C7_2"
		call seek_card(card_num)
		
		maxtran2 = 0
		do i=1,maxtran ! currently, maxtran = 2
			read(pw_main_inp,*) is_tran(i), num_tran_ser(i), tran_ser_shape(i)
			if(is_tran(i) == 1) then
				maxtran2 = maxtran2 + 1
			end if
		end do

		is_salt = is_tran(1)
		is_temp = is_tran(2)
		num_salt_ser = num_tran_ser(1)
		num_temp_ser = num_tran_ser(2)
		salt_ser_shape = tran_ser_shape(1)
		temp_ser_shape = tran_ser_shape(2)
		
		! allocate tran_id and initialize it
		! tran_id will be find in read_main_inp.f90
		allocate(tran_id(maxtran2))
		tran_id = 0
		
		if(is_tran(1) == 1 .and. num_salt_ser > 0) then
			call scan_salt_ser ! update: max_salt_data_num
		end if
			
		if(is_tran(2) == 1 .and. num_temp_ser > 0) then
			call scan_temp_ser ! update: max_temp_data_num
		end if
		
		if(is_tran(2) == 1) then
			call scan_air_ser ! update: max_air_data_num
		end if		
	end if
	
	max_salt_ser 		= max(1,max_salt_ser)
	max_temp_ser 		= max(1,max_temp_ser)
	max_salt_data_num = max(1,max_salt_data_num)
	max_temp_data_num = max(1,max_temp_data_num)
	max_air_data_num  = max(1,max_air_data_num)
	
	write(*,'(A30,I10)') 'max_salt_ser = ', max_salt_ser
	write(*,'(A30,I10)') 'max_temp_ser = ', max_temp_ser
	write(*,'(A30,I10)') 'max_salt_data_num = ', max_salt_data_num
	write(*,'(A30,I10)') 'max_temp_data_num = ', max_temp_data_num
	write(*,'(A30,I10)') 'max_air_data_num = ', max_air_data_num
	
!=============================================================================!
! C20 ~ C29, Tidal & River Boundary Conditions											!
!=============================================================================!
	! find: max_ob_node, max_ob_element, max_ob_face --------------------------!
	card_num = "C20"
	call seek_card(card_num)
	read(pw_main_inp,*) num_ob_cell, ob_info_option
	read(pw_main_inp,*) num_tide_bc, num_harmonic_ser, check_tide
	read(pw_main_inp,*) num_eta_bc, num_eta_ser, check_eta, eta_ser_shape

	if(num_ob_cell /= (num_tide_bc + num_eta_bc)) then
		write(*,*) 'num_tide_bc is not equal to (num_tide_bc + num_eta_bc)'
		write(*,*) 'scan_gom.f90'
		stop
	end if
	
	if(num_ob_cell > 0 .and. num_tide_bc > 0 .and. num_harmonic_ser > 0) then
		! find maximum numuber of tidal constituent which will be using:
		! 		"max_no_tidal_constituent"
		call scan_harmonic_ser 
	end if
	
	if(num_ob_cell > 0 .and. num_eta_bc > 0 .and. num_eta_ser > 0) then
		! find "max_eta_data_num" and "max_eta_ser"
		call scan_eta_ser
	end if
	
	! find max_ob_node
	if(num_ob_cell > 0) then
		allocate(nd(num_ob_cell,2))
		nd = 0
		node_count = 0
		
		if(ob_info_option == 0) then
			card_num = "C20_1"
			call seek_card(card_num)
			do i=1,num_ob_cell
				read(pw_main_inp,*) serial_num, nd(i,1), nd(i,2)
				node_count = node_count + 2
				
				if(i>1) then
					do j=1,(i-1)
						if(nd(i,1) == nd(j,1)) then
							node_count = node_count -1
						else if(nd(i,1) == nd(j,2)) then
							node_count = node_count -1
						else if(nd(i,2) == nd(j,1)) then
							node_count = node_count -1
						else if(nd(i,2) == nd(j,2)) then
							node_count = node_count -1
						end if
					end do
				end if
			end do
		else if(ob_info_option == 1) then
			open(pw_Obc_info,file=id_Obc_info,form='formatted',status='old')

			! skip header lines 
			call skip_header_lines(pw_Obc_info,id_Obc_info)
						
			do i=1,num_ob_cell
				read(pw_Obc_info,*) serial_num, nd(i,1), nd(i,2)
				node_count = node_count + 2
				
				if(i>1) then
					do j=1,(i-1)
						if(nd(i,1) == nd(j,1)) then
							node_count = node_count -1
						else if(nd(i,1) == nd(j,2)) then
							node_count = node_count -1
						else if(nd(i,2) == nd(j,1)) then
							node_count = node_count -1
						else if(nd(i,2) == nd(j,2)) then
							node_count = node_count -1
						end if
					end do
				end if
			end do
			close(pw_Obc_info)
		end if
		deallocate(nd)
	end if
		
	max_ob_node 					= MAX(1,node_count)
	max_ob_element 				= MAX(1,num_ob_cell)
	max_ob_face						= MAX(1,num_ob_cell)
	max_no_tidal_constituent	= MAX(1,max_no_tidal_constituent)
	max_eta_ser 					= MAX(1,max_eta_ser)
	max_eta_data_num 				= MAX(1,max_eta_data_num)
	
	write(*,'(A30,I10)') 'max_ob_node = ', max_ob_node
	write(*,'(A30,I10)') 'max_ob_element = ', max_ob_element
	write(*,'(A30,I10)') 'max_ob_face = ', max_ob_face
	write(*,'(A30,I10)') 'max_no_tidal_constituent = ', max_no_tidal_constituent ! from "scan_harmonic_ser"
	write(*,'(A30,I10)') 'max_eta_ser = ', max_eta_ser
	write(*,'(A30,I10)') 'max_eta_data_num = ', max_eta_data_num
	
	! find: max_Q_bc, max_Q_ser, max_q_data_num -------------------------------!
	card_num = "C22"
	max_q_data_num = 1
	call seek_card(card_num)
	read(pw_main_inp,*) num_Qb_cell, Qbc_info_option, num_Q_ser, check_Q, Q_ser_shape
	max_Q_bc 	= MAX(1,num_Qb_cell)
	max_Q_ser 	= MAX(1,num_Q_ser)
	
	if(num_Qb_cell > 0) then
		if(num_Q_ser > 0) then
			call scan_q_ser ! update max_q_data_num
		end if
	end if
	max_q_data_num = MAX(1,max_q_data_num)
	
	write(*,'(A30,I10)') 'max_Q_bc = ', max_Q_bc
	write(*,'(A30,I10)') 'max_Q_ser = ', max_Q_ser
	write(*,'(A30,I10)') 'max_q_data_num = ', max_q_data_num

	! find: max_WR_bc ---------------------------------------------------------!
	card_num = "C23"
	call seek_card(card_num)
	read(pw_main_inp,*) num_WR_cell
	max_WR_bc = MAX(1,num_WR_cell)
	write(*,'(A30,I10)') 'max_WR_bc = ', max_WR_bc
	
	! find: max_SS_bc ---------------------------------------------------------!
	card_num = "C24"
	call seek_card(card_num)
	read(pw_main_inp,*) num_SS_cell
	max_SS_bc = MAX(1,num_SS_cell)
	write(*,'(A30,I10)') 'max_SS_bc = ', max_SS_bc
	
	! find: max_Etide_species -------------------------------------------------!
	card_num = "C25"
	call seek_card(card_num)
	read(pw_main_inp,*) no_Etide_species, Etide_cutoff_depth
	max_Etide_species = MAX(1,no_Etide_species)
	write(*,'(A30,I10)') 'max_Etide_species = ', max_Etide_species	

!=============================================================================!
! C30 ~ C39, Air/Sea Boundary Conditions													!
!=============================================================================!				
	! find: max_windp_station -------------------------------------------------!
	card_num = "C30"
	call seek_card(card_num)
	read(pw_main_inp,*) wind_flag, airp_flag, num_windp_ser

	max_windp_station = MAX(1,num_windp_ser)
	max_windp_data_num = 1 ! to avoid allocation of 0 length array

	! find: max_windp_data_num
	if(wind_flag == 1 .or. airp_flag == 1) then
		if(num_windp_ser > 0) then
			call scan_windp_ser
		end if
	end if
	write(*,'(A30,I10)') 'max_windp_station = ', max_windp_station
	write(*,'(A30,I10)') 'max_windp_data_num = ', max_windp_data_num
	
!=============================================================================!
! C40 ~ C49: Output Control Options															!
!=============================================================================!	
	! find: max_tser_station --------------------------------------------------!
	card_num = "C41"
	call seek_card(card_num)
	read(pw_main_inp,*) tser_station_num
	max_tser_station = MAX(1,tser_station_num)
	write(*,'(A30,I10)') 'max_tser_station = ', max_tser_station
! End of reading main.inp ====================================================!
! ============================================================================!	
! 	write(*,*) '.'
! 	write(*,*) '.'
! 	write(*,*) '.'
  	close(pw_main_inp)		! Close 'main.inp'
end subroutine scan_gom