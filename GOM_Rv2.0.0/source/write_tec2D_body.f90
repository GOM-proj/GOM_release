!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_tec2D_body
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, j
	integer :: quotient, remainder
	character(len=15) :: zonetype = 'FEQUADRILATERAL'
	character(len= 4) :: time_char
	character(len=40) :: format1, format2
	! jw
	
	! jw
	write(pw_tec2D,'(A,F15.5, A,I10, A,I10, A,A, A, F15.5, A)') &
	&	'ZONE T = "', julian_day * IS2D_time_conv, &
	&	'", N =', MAXNOD, &
	&	', E = ', MAXELE, &
	&	', DATAPACKING=BLOCK, ZONETYPE=',zonetype, &
	&	', VARSHARELIST=([1-2]=1), CONNECTIVITYSHAREZONE=1, SOLUTIONTIME=', &
	&	julian_day * IS2D_time_conv, ', &
	&	STRANDID=1' 
			
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
	if(IS2D_variable(1) == 1) then
		do i=1,quotient
			write(pw_tec2D,format1) (eta_node((i-1)*100+j), j=1,100)
		end do
		if(remainder > 0) then
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
end subroutine write_tec2D_body