!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_tec3D_surf_body
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
	write(pw_tec3D_surf,'(A,F15.5,A,I10,A,I10,A,A,A,F15.5,A)') &
	& 'ZONE T = "',elapsed_time * IS3D_time_conv,'", N =', MAXNOD, ', E = ', MAXELE, &
	& ', DATAPACKING=BLOCK, ZONETYPE=',zonetype,', VARSHARELIST=([1-2]=1), &
	& CONNECTIVITYSHAREZONE=1, SOLUTIONTIME=', julian_day * IS3D_time_conv,', STRANDID=1'
	
	quotient = int(maxnod/100)
	remainder = mod(maxnod,100)

	! jw
	write(format1,'(A,I0,A)') '(',100,'E15.5E4)'
	write(format2,'(A,I0,A)') '(',remainder,'E15.5E4)'
	
	! jw
	do i=1,quotient
		write(pw_tec3D_surf,format1) (eta_node((i-1)*100+j), j=1,100)
	end do
	if(remainder > 0) then
		write(pw_tec3D_surf,format2) (eta_node(quotient*100+j), j=1,remainder)
	end if
	
	! jw
	if(IS3D_time == 1) then	! jw
		time_char = ' sec'
	else if(IS3D_time == 2) then ! jw
		time_char = ' min'
	else if(IS3D_time == 3) then ! jw
		time_char = '  hr'
	else if(IS3D_time == 4) then ! jw
		time_char = ' day'
	end if
	write(pw_tec3D_surf,'(A,F15.5,A,A,I5)') &
	&	'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, T="TIME = ', &
	&	julian_day * IS3D_time_conv, time_char, '", ZN =', zone_num_3D_surf

	write(pw_tec3D_surf,*)	! jw
	
	if(IS3D_binary == 1) then
		! jw
	end if	
end subroutine write_tec3D_surf_body