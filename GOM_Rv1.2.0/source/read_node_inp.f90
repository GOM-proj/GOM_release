!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! This subroutine reads node.inp
!! 
subroutine read_node_inp
	use mod_global_variables
	use mod_file_definition	
	
	implicit none
	integer :: i																			    	 ! line number in a file
	integer :: node_num	! number of nodes, this number should be identical with maxnod
	integer :: ibuff
	real(dp):: rtemp
	real(dp):: lon_min, lon_max, lat_min, lat_max
	real(dp):: eastern_boundary, western_boundary
	integer :: utm_zone
   ! End of local variables ==================================================!
	   
   ! allocate & initialize some variables ====================================!
	allocate(x_node(maxnod), 					& ! [m]
	&			y_node(maxnod),					& ! [m]
	&			h_node(maxnod),					& ! [m]
	&			initial_wetdry_node(maxnod),	& ! 0/1
	&			lon_node(maxnod),					& 
	&			lat_node(maxnod))
	x_node = 0.0_dp
	y_node = 0.0_dp
	h_node = 0.0_dp
	lon_node = 0.0_dp
	lat_node = 0.0_dp
	
	! If the node is wet, initial_wetdry_node(i) = 0  (i.e., false)
	! If the node is dry, initial_wetdry_node(i) = 1  (i.e., true)
   initial_wetdry_node = 0		! initially set to wet.
	! end of allocation =======================================================!
	
	
	write(pw_run_log,*) "	Read node.inp"
	write(pw_run_log,*) "		Now, you are in 'read_input.f90 -> subroutine read_node_inp"
	
	open(pw_node_inp, file = id_node_inp, form='formatted', status = 'old')
	
	if(node_mirr == 1) then
		open(pw_node_mirr, file = id_node_mirr, form = 'formatted', status = 'replace')
	end if
	
	! skip header lines 
	call skip_header_lines(pw_node_inp,id_node_inp)
	
	! read main body ==========================================================!
	! read total_node_number
	read(pw_node_inp,*) node_num
		
	if(node_mirr == 1) then
		write(pw_node_mirr,*) node_num
	end if
	
	if(node_num /= MAXNOD) then	
		write(pw_run_log,*) "Total node number does not mathch: STOP"
		stop "Total node number does not mathch: STOP"
	end if

	! read main body of node.inp ==============================================!
	if(coordinate_system == 1) then			! Cartesian
		do i=1,maxnod
			read(pw_node_inp,*) ibuff, x_node(i), y_node(i), rtemp ! h_node(i)
			h_node(i) = rtemp + h_node_adjust ! h_node_adjust is defined in C1 in main.inp
			
			! set dry node ------------------------------------------------------!
	      if(h_node(i) <= 0.0) then
	         initial_wetdry_node(i) = 1
	      end if
		end do
	else if(coordinate_system == 2) then	! Lon/Lat
		! in this case, I need two steps: (1) read data, (2) find mid lon/lat and do other process
		! first, just read the data
		do i=1,maxnod
         read(pw_node_inp,*) ibuff, lon_node(i), lat_node(i), rtemp ! h_node(i)
         h_node(i) = rtemp + h_node_adjust ! h_node_adjust is defined in C1 in main.inp
         
			! set dry node ------------------------------------------------------!
	      if(h_node(i) <= 0.0) then
	         initial_wetdry_node(i) = 1
	      end if
		end do
		
		! second, find the mid lon/lat value -----------------------------------!
		lon_min = minval(lon_node)
		lon_max = maxval(lon_node)
		lat_min = minval(lat_node)
		lat_max = maxval(lat_node)
		
		! mid-point of the model domain:
		lon_mid = (lon_min + lon_max) * 0.5_dp
		lat_mid = (lat_min + lat_max) * 0.5_dp
		
		! find UTM zone: -------------------------------------------------------!
		do utm_zone=1,60 ! UTM zone has 1 ~ 60 (each 6 degree)
			eastern_boundary = (utm_zone*6.0_dp) - 180.0_dp 	! in [degree]
			western_boundary = eastern_boundary - 6.0_dp 		! in [degree]
			
			if(lon_mid >= western_boundary .and. lon_mid <= eastern_boundary) then
				! if mid-point of the model domain is in between eastern and western boundary of the UTM zone,
				! it is the zone for our model domain.
				utm_projection_zone = utm_zone
				exit
			end if
		end do
		
		! now, convert given lon/lat to UTM coordinate (i.e., in Cartesian coordinate)
		! conversion option in coordinate_conversion.f90:
		! 		1: lon/lat -> UTM
		! 		2: UTM -> lon/lat
		!$omp parallel do private(i)
		do i=1,maxnod
			call coordinate_conversion(lon_node(i),lat_node(i),utm_projection_zone,1, x_node(i),y_node(i))
		end do
		!$omp end parallel do
	end if	
	
	xn_min = minval(x_node)	! use this value to shift origin to zero
	yn_min = minval(y_node)	! use this value to shift origin to zero
	
   close(pw_node_inp)		! Close 'node.inp'
   ! end of reading main body of node.inp ====================================!
   
	! overwrite water depth with analytical information if ana_depth is on ====!
	! this will overwrite:
	! 		h_node & initial_wetdry_node
	if(ana_depth == 1) then
		call ana_depth_adjustment
	end if
   
	! write mirror image file, node_mirr.dat ==================================!
	if(node_mirr == 1 .and. ana_depth == 0) then
		! if ana_depth == 1, this information will be written at ana_depth.f90
		if(coordinate_system == 1) then
			do i=1,maxnod
				write(pw_node_mirr,*) i, x_node(i), y_node(i), h_node(i)
			end do
		else if(coordinate_system == 2) then
			do i=1,maxnod
				write(pw_node_mirr,*) i, lon_node(i), lat_node(i), h_node(i)
			end do
		end if
			
		close(pw_node_mirr)
	end if	
end subroutine read_node_inp
