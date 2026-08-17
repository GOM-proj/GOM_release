!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_tec2D_head
	use mod_global_variables
	use mod_file_definition, only : pw_tec2D
	
	implicit none
	integer :: i, j, k
	integer :: quotient, remainder
   character(len=15) :: zonetype = 'FEQUADRILATERAL'
	character(len= 4) :: time_char
	character(len=40) :: format1, format2
	! jw
	
	! jw
	write(pw_tec2D,'(A)') 'Title = "2D contour plot"'
	
	! jw
	! jw
	! jw
	if(IS2D_unit_conv > 0.5_dp) then	! jw
		write(pw_tec2D,'(A)',advance = 'no') 'Variables = "X [m]", "Y [m]"'
	else										! jw
		write(pw_tec2D,'(A)',advance = 'no') 'Variables = "X [km]", "Y [km]"'
	end if
	
	do i=1,6
		if(IS2D_variable(i) == 1 .and. i == 1) then
			write(pw_tec2D,'(A)',advance = 'no') ', "eta [m]"'
		else if(IS2D_variable(i) == 1 .and. i == 2) then
			write(pw_tec2D,'(A)',advance = 'no') ', "u [m/s]"'
		else if(IS2D_variable(i) == 1 .and. i == 3) then
			write(pw_tec2D,'(A)',advance = 'no') ', "v [m/s]"'
		else if(IS2D_variable(i) == 1 .and. i == 4) then
			write(pw_tec2D,'(A)',advance = 'no') ', "salt [ppt]"'
		else if(IS2D_variable(i) == 1 .and. i == 5) then
			write(pw_tec2D,'(A)',advance = 'no') ', "temp [C]"'
		else if(IS2D_variable(i) == 1 .and. i == 6) then
			write(pw_tec2D,'(A)',advance = 'no') ', "rho [kg/m3]"'
		end if
	end do
	write(pw_tec2D,*) ! jw
	write(pw_tec2D,'(A,F15.5, A,I10, A,I10, A,A, A,F15.5, A)') &
	& 	'ZONE T = "', julian_day * IS2D_time_conv, &
	&	'", N =', MAXNOD, &
	&	', E = ', 	MAXELE, &
	&	', DATAPACKING=BLOCK, ZONETYPE=', zonetype, &
 	&	', SOLUTIONTIME=', julian_day * IS2D_time_conv, &
 	&	', STRANDID=1' 
	
	! jw
	quotient = int(maxnod/100)
	remainder = mod(maxnod,100)
	
	! jw
	! jw
	! jw
	! jw
	
	! jw
	write(format1,'(A,I0,A)') '(',100,'(F0.5,1x))'
	write(format2,'(A,I0,A)') '(',remainder,'(F0.5,1x))'

	! jw
	! jw
	
	! jw
	do i=1,quotient
		write(pw_tec2D,format1) ((x_node((i-1)*100+j)-xn_min) * IS2D_unit_conv, j=1,100) ! jw
	end do
	if(remainder > 0) then
		write(pw_tec2D,format2) ((x_node(quotient*100+j)-xn_min) * IS2D_unit_conv, j=1,remainder) ! jw
	end if
	write(pw_tec2D,*) ! jw

	! jw
	do i=1,quotient
		write(pw_tec2D,format1) ((y_node((i-1)*100+j)-yn_min) * IS2D_unit_conv, j=1,100) ! jw
	end do
	if(remainder > 0) then
		write(pw_tec2D,format2) ((y_node(quotient*100+j)-yn_min) * IS2D_unit_conv, j=1,remainder) ! jw
	end if
	write(pw_tec2D,*) ! jw
	
	! jw
	! jw
	if(IS2D_variable(1) == 1) then
		do i=1,quotient
			! jw
			write(pw_tec2D,format1) (eta_node((i-1)*100+j), j=1,100)
		end do
		if(remainder > 0) then
			! jw
			write(pw_tec2D,format2) (eta_node(quotient*100+j), j=1,remainder)
		end if
		write(pw_tec2D,*) ! jw
	end if
	
	! jw
	if(IS2D_variable(2) == 1) then
		do i=1,quotient
			write(pw_tec2D,format1) (ubar_node((i-1)*100+j), j=1,100)
		end do
		if(remainder > 0) then
			write(pw_tec2D,format2) (ubar_node(quotient*100+j), j=1,remainder)
		end if
		write(pw_tec2D,*) ! jw
	end if
	
	! jw
	if(IS2D_variable(3) == 1) then
		do i=1,quotient
			write(pw_tec2D,format1) (vbar_node((i-1)*100+j), j=1,100)
		end do
		if(remainder > 0) then
			write(pw_tec2D,format2) (vbar_node(quotient*100+j), j=1,remainder)
		end if
		write(pw_tec2D,*) ! jw
	end if
	
	! jw
	if(IS2D_variable(4) == 1) then
		do i=1,quotient
			write(pw_tec2D,format1) (sbar_node((i-1)*100+j), j=1,100)
		end do
		if(remainder > 0) then
			write(pw_tec2D,format2) (sbar_node(quotient*100+j), j=1,remainder)
		end if
		write(pw_tec2D,*) ! jw
	end if
	
	! jw
	if(IS2D_variable(5) == 1) then
		do i=1,quotient
			write(pw_tec2D,format1) (tbar_node((i-1)*100+j), j=1,100)
		end do
		if(remainder > 0) then
			write(pw_tec2D,format2) (tbar_node(quotient*100+j), j=1,remainder)
		end if
		write(pw_tec2D,*) ! jw
	end if

	! jw
	if(IS2D_variable(6) == 1) then
		do i=1,quotient
			write(pw_tec2D,format1) (rbar_node((i-1)*100+j), j=1,100)
		end do
		if(remainder > 0) then
			write(pw_tec2D,format2) (rbar_node(quotient*100+j), j=1,remainder)
		end if
		write(pw_tec2D,*) ! jw
	end if

	! jw
		
	! jw
	do i=1,MAXELE
		write(pw_tec2D,'(4I10)') (nodenum_at_cell_tec(k,i),k=1,4)
	end do
	
	! jw
	if(IS2D_time == 1) then	! jw
		time_char = ' sec'
	else if(IS2D_time == 2) then ! jw
		time_char = ' min'
	else if(IS2D_time == 3) then ! jw
		time_char = '  hr'
	else if(IS2D_time == 4) then ! jw
		time_char = ' day'
	end if
	write(pw_tec2D,'(A,F15.5,A,A,I5)') &
	&	'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, T="TIME = ', &
	&	julian_day * IS2D_time_conv, time_char, '", ZN =', zone_num_2D

	write(pw_tec2D,*)	! jw
end subroutine write_tec2D_head