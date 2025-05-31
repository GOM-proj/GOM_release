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
	integer :: i																			    	 ! jw
	integer :: node_num	! jw
	integer :: ibuff
	real(dp):: rtemp
	real(dp):: lon_min, lon_max, lat_min, lat_max
	real(dp):: eastern_boundary, western_boundary
	integer :: utm_zone
   ! jw
	   
   ! jw
	allocate(x_node(maxnod), 					& ! jw
	&			y_node(maxnod),					& ! jw
	&			h_node(maxnod),					& ! jw
	&			initial_wetdry_node(maxnod),	& ! jw
	&			lon_node(maxnod),					& 
	&			lat_node(maxnod))
	x_node = 0.0_dp
	y_node = 0.0_dp
	h_node = 0.0_dp
	lon_node = 0.0_dp
	lat_node = 0.0_dp
	
	! jw
	! jw
   initial_wetdry_node = 0		! jw
	! jw
	
	
	write(pw_run_log,*) "	Read node.inp"
	write(pw_run_log,*) "		Now, you are in 'read_input.f90 -> subroutine read_node_inp"
	
	open(pw_node_inp, file = id_node_inp, form='formatted', status = 'old')
	
	if(node_mirr == 1) then
		open(pw_node_mirr, file = id_node_mirr, form = 'formatted', status = 'replace')
	end if
	
	! jw
	call skip_header_lines(pw_node_inp,id_node_inp)
	
	! jw
	! jw
	read(pw_node_inp,*) node_num
		
	if(node_mirr == 1) then
		write(pw_node_mirr,*) node_num
	end if
	
	if(node_num /= MAXNOD) then	
		write(pw_run_log,*) "Total node number does not mathch: STOP"
		stop "Total node number does not mathch: STOP"
	end if

	! jw
	if(coordinate_system == 1) then			! jw
		do i=1,maxnod
			read(pw_node_inp,*) ibuff, x_node(i), y_node(i), rtemp ! jw
			h_node(i) = rtemp + h_node_adjust ! jw
			
			! jw
	      if(h_node(i) <= 0.0) then
	         initial_wetdry_node(i) = 1
	      end if
		end do
	else if(coordinate_system == 2) then	! jw
		! jw
		! jw
		do i=1,maxnod
         read(pw_node_inp,*) ibuff, lon_node(i), lat_node(i), rtemp ! jw
         h_node(i) = rtemp + h_node_adjust ! jw
         
			! jw
	      if(h_node(i) <= 0.0) then
	         initial_wetdry_node(i) = 1
	      end if
		end do
		
		! jw
		lon_min = minval(lon_node)
		lon_max = maxval(lon_node)
		lat_min = minval(lat_node)
		lat_max = maxval(lat_node)
		
		! jw
		lon_mid = (lon_min + lon_max) * 0.5_dp
		lat_mid = (lat_min + lat_max) * 0.5_dp
		
		! jw
		do utm_zone=1,60 ! jw
			eastern_boundary = (utm_zone*6.0_dp) - 180.0_dp 	! jw
			western_boundary = eastern_boundary - 6.0_dp 		! jw
			
			if(lon_mid >= western_boundary .and. lon_mid <= eastern_boundary) then
				! jw
				! jw
				utm_projection_zone = utm_zone
				exit
			end if
		end do
		
		! jw
		! jw
		! jw
		! jw
		!$omp parallel do private(i)
		do i=1,maxnod
			call coordinate_conversion(lon_node(i),lat_node(i),utm_projection_zone,1, x_node(i),y_node(i))
		end do
		!$omp end parallel do
	end if	
	
	xn_min = minval(x_node)	! jw
	yn_min = minval(y_node)	! jw
	
   close(pw_node_inp)		! jw
   ! jw
   
	! jw
	! jw
	! jw
	if(ana_depth == 1) then
		call ana_depth_adjustment
	end if
   
	! jw
	if(node_mirr == 1 .and. ana_depth == 0) then
		! jw
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
