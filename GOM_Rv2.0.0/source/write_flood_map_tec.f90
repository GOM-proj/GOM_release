!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_flood_map_tec
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i,k
	character(len=15) :: zonetype = 'FEQUADRILATERAL'
	! jw
	
	! jw
	write(pw_flood_map_tec,'(A)') 'Title = "2D flood map contour plot"'
	if(IS2D_unit_conv == 1.0_dp) then
		write(pw_flood_map_tec,'(A)',advance = 'no') 'Variables = "X [m]", "Y [m]"'
	else if(IS2D_unit_conv == 0.001_dp) then
		write(pw_flood_map_tec,'(A)',advance = 'no') 'Variables = "X [km]", "Y [km]"'
	end if
	write(pw_flood_map_tec,'(A)') ', "max_eta", "max_flood_time", "flood_id"'
	
	
	! jw
	write(pw_flood_map_tec,'(A,F15.5,A,I10,A,I10,A,A,A,F15.5,A)') &
	&	'ZONE T = "', elapsed_time * IS2D_time_conv, '", N =', MAXNOD, ', E = ', MAXELE, ', DATAPACKING=POINT, ZONETYPE=',zonetype,	&
	&	', SOLUTIONTIME=', elapsed_time * IS2D_time_conv, ', STRANDID=1' 
	
	! jw
	do i=1,maxnod
		write(pw_flood_map_tec,'(4F15.5, I3)') &
		&	x_node(i) * IS2D_unit_conv, y_node(i) * IS2D_unit_conv, max_eta_node(i), max_flood_time(i), flood_id(i)
	end do
	
	! jw
	do i=1,MAXELE
		write(pw_flood_map_tec,'(4I10)') (nodenum_at_cell_tec(k,i),k=1,4)
	end do
	
	! jw
	if(IS2D_time == 1) then	! jw
		write(pw_flood_map_tec,'(A,F15.5,A,I5)') 'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, &
		&	T="TIME = ',elapsed_time * IS2D_time_conv,' sec", ZN =', 1
	else if(IS2D_time == 2) then ! jw
		write(pw_flood_map_tec,'(A,F15.5,A,I5)') 'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, &
		&	T="TIME = ',elapsed_time * IS2D_time_conv,' min", ZN =', 1
	else if(IS2D_time == 3) then ! jw
		write(pw_flood_map_tec,'(A,F15.5,A,I5)') 'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, &
		&	T="TIME = ',elapsed_time * IS2D_time_conv,' hr", ZN =', 1
	else if(IS2D_time == 4) then ! jw
		write(pw_flood_map_tec,'(A,F15.5,A,I5)') 'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, &
		&	T="TIME = ',elapsed_time * IS2D_time_conv,' day", ZN =', 1
	end if
	
	close(pw_flood_map_tec)
end subroutine write_flood_map_tec