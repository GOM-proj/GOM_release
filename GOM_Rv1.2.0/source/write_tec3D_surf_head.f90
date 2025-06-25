!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_tec3D_surf_head
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i,j,k
	integer :: quotient, remainder
	character(len=15) :: zonetype = 'FEQUADRILATERAL'
	character(len= 4) :: time_char
	character(len=40) :: format1, format2
	! End of local variables ==================================================!
	
	! Write header lines ======================================================!
	write(pw_tec3D_surf,'(A)') 'Title = "3D contour plot"'

	! note: do not use:
	!   	if(check_grid_unit_conv == 0.1)
	! 		since, real number annnot be equal to other real number
	if(IS3D_unit_conv > 0.5_dp) then	! 1.0
		write(pw_tec3D_surf,'(A)') 'Variables = "X [m]", "Y [m]", "Z [m]"'
	else										! 0.001
		write(pw_tec3D_surf,'(A)') 'Variables = "X [km]", "Y [km]", "Z [m]"'
	end if
	
	write(pw_tec3D_surf,'(A,F15.5,A,I10,A,I10,A,A,A,F15.5,A)') 'ZONE T = "', julian_day * IS3D_time_conv,	&
	&	'", N =', MAXNOD, ', E = ', MAXELE, ', DATAPACKING=BLOCK, ZONETYPE=',zonetype,	&
	&	', SOLUTIONTIME=', julian_day * IS3D_time_conv, ', STRANDID=1' 
		
	! Write x & y coordinates & variables =====================================!
	quotient = int(MAXNOD/100)
	remainder = mod(MAXNOD,100)

	! prepare the format
	write(format1,'(A,I0,A)') '(',100,'E15.5E4)' ! ex, (100F15.5)
	write(format2,'(A,I0,A)') '(',remainder,'E15.5E4)' ! ex, (100F15.5)
	
	! write x coordiante
	do i=1,quotient
		write(pw_tec3D_surf,format1) ((x_node((i-1)*100+j)-xn_min) * IS3D_unit_conv, j=1,100) ! shift origin to zero
	end do
	if(remainder > 0) then
		write(pw_tec3D_surf,format2) ((x_node(quotient*100+j)-xn_min) * IS3D_unit_conv, j=1,remainder) ! shift origin to zero
	end if
	write(pw_tec3D_surf,*) ! put a line
		
	! write y coordiante
	do i=1,quotient
		write(pw_tec3D_surf,format1) ((y_node((i-1)*100+j)-yn_min) * IS3D_unit_conv, j=1,100) ! shift origin to zero
	end do
	if(remainder > 0) then
		write(pw_tec3D_surf,format2) ((y_node(quotient*100+j)-yn_min) * IS3D_unit_conv, j=1,remainder) ! shift origin to zero
	end if
	write(pw_tec3D_surf,*) ! put a line
		
	! Write variabels 
	do i=1,quotient
		write(pw_tec3D_surf,format1) (eta_node((i-1)*100+j), j=1,100)
	end do
	if(remainder > 0) then
		write(pw_tec3D_surf,format2) (eta_node(quotient*100+j), j=1,remainder)
	end if
	write(pw_tec3D_surf,*) ! put a line
	
	! Write cell connectivity =================================================!
	do i=1,MAXELE
		write(pw_tec3D_surf,'(4I10)') (nodenum_at_cell_tec(k,i),k=1,4)
	end do

	! Write text information at the bottom of each zone =======================!
	if(IS3D_time == 1) then	! second
		time_char = ' sec'
	else if(IS3D_time == 2) then ! minute
		time_char = ' min'
	else if(IS3D_time == 3) then ! hour
		time_char = '  hr'
	else if(IS3D_time == 4) then ! day
		time_char = ' day'
	end if
	write(pw_tec3D_surf,'(A,F15.5,A,A,I5)') &
	&	'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, T="TIME = ', &
	&	julian_day * IS3D_time_conv, time_char, '", ZN =', zone_num_3D_surf

	write(pw_tec3D_surf,*)	! put a line					

end subroutine write_tec3D_surf_head