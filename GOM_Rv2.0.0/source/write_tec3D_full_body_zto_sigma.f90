!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! This subroutine is almost identical to write_tec3D_full_head_zto_sigma.f90
!!  
subroutine write_tec3D_full_body_zto_sigma
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, j, k
	integer :: count1, count2, count3
	integer :: quotient_node, remainder_node, quotient_cell, remainder_cell	
	character(len=15) :: zonetype = 'FEBRICK'
	character(len= 4) :: time_char
	character(len=40) :: format1, format2_cell, format2_node, format3
	
	! jw
	real(dp),dimension(0:maxlayer,maxnod) :: u_sigma, v_sigma, w_sigma 				! jw
	real(dp),dimension(maxlayer,maxele) :: salt_sigma, temp_sigma, rho_sigma 	! jw
	! jw

	! jw
	! jw
	! jw
	! jw
		
	! jw
	call zto_sigma_node_values(u_sigma, v_sigma, w_sigma)

	! jw
	call zto_sigma_cell_values(salt_sigma, temp_sigma, rho_sigma)
	! jw

	! jw
	! jw
	quotient_node = int(maxnod/100)
	remainder_node = mod(maxnod,100)
	quotient_cell = int(maxele/100)
	remainder_cell = mod(maxele,100)

	! jw
! jw
! jw
! jw

! jw
! jw
! jw
! jw

	write(format1,'(A,I0,A)') '(',100,'(F0.5,1x))'
	write(format2_cell,'(A,I0,A)') '(',remainder_cell,'(F0.5,1x))'
	write(format2_node,'(A,I0,A)') '(',remainder_node,'(F0.5,1x))'
	write(format3,'(A,I0,A)') '(',maxlayer,'(F0.5,1x))'
	
	! jw
	! jw
	
	! jw
	! jw
	write(pw_tec3D_full,'(A,F15.5, A,I10, A,I10, A)', advance = 'no') &
	&	'ZONE T = "', julian_day * IS3D_time_conv,	&
	&	'", N =', MAXNOD*(maxlayer+1), &
	&	', E = ', MAXELE*maxlayer, &
	&	', DATAPACKING=BLOCK'

	! jw
	count1 = 0
	do i=1,3
		if(IS3D_variable(i) == 1) then
			count1 = count1+1
		end if
	end do
	
	! jw
	! jw
	! jw
	! jw
	count2 = 0
	do i=4,6
		if(IS3D_variable(i) == 1) then
			count2 = count2 + 1
		end if
	end do
	
	if(count2 > 0) then
		write(pw_tec3D_full,'(A)',advance = 'no') &
		&	', VARLOCATION=(['
		count3 = 0
		do i=1,3
			! jw
			! jw
			if(IS3D_variable(i+3) == 1) then
				count3 = count3 + 1
				write(pw_tec3D_full,'(I1,A)',advance = 'no') (3+count1)+count3, ',' ! jw
			end if
			! jw
		end do
		write(pw_tec3D_full,'(A)',advance = 'no') &
		&	']=CELLCENTERED)'
	end if
	
	write(pw_tec3D_full,'(A,A, A,F15.5, A,F15.5, A)') &
	&	', ZONETYPE=',zonetype, &
	&	', SOLUTIONTIME=', julian_day * IS3D_time_conv, &
	&	', VARSHARELIST=([1-2]=1), CONNECTIVITYSHAREZONE=1, SOLUTIONTIME=', julian_day * IS3D_time_conv, & ! jw
	&	', STRANDID=1' 
	
	! jw
	! jw
	! jw
	do k=0,maxlayer ! jw
		do i=1,quotient_node
			write(pw_tec3D_full,format1) ( &
			&	-h_node((i-1)*100+j) + ((eta_node((i-1)*100+j) + h_node((i-1)*100+j))/maxlayer)*k,  &
			&	j=1,100 )
		end do
		if(remainder_node > 0) then
			write(pw_tec3D_full,format2_node) ( &
			&	-h_node(quotient_node*100+j) + ((eta_node(quotient_node*100+j) + h_node(quotient_node*100+j))/maxlayer)*k,  &
			&	j=1,remainder_node )
		end if
	end do
	write(pw_tec3D_full,*)	! jw

	! jw
	if(IS3D_variable(1) == 1) then
		do k=0,maxlayer
			do i=1,quotient_node
				write(pw_tec3D_full,format1) (u_sigma(k,(i-1)*100+j), j=1,100)
			end do
			if(remainder_node > 0) then
				write(pw_tec3D_full,format2_node) (u_sigma(k,quotient_node*100+j), j=1,remainder_node)
			end if
		end do
		write(pw_tec3D_full,*)	! jw
	end if

	! jw
	if(IS3D_variable(2) == 1) then
		do k=0,maxlayer
			do i=1,quotient_node
				write(pw_tec3D_full,format1) (v_sigma(k,(i-1)*100+j), j=1,100)
			end do
			if(remainder_node > 0) then
				write(pw_tec3D_full,format2_node) (v_sigma(k,quotient_node*100+j), j=1,remainder_node)
			end if
		end do
		write(pw_tec3D_full,*)	! jw
	end if

	! jw
	if(IS3D_variable(3) == 1) then
		do k=0,maxlayer
			do i=1,quotient_node
				write(pw_tec3D_full,format1) (w_sigma(k,(i-1)*100+j), j=1,100)
			end do
			if(remainder_node > 0) then
				write(pw_tec3D_full,format2_node) (w_sigma(k,quotient_node*100+j), j=1,remainder_node)
			end if
		end do
		write(pw_tec3D_full,*)	! jw
	end if
	
	! jw
	if(IS3D_variable(4) == 1) then
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
		do i=1,maxele
			write(pw_tec3D_full,format3) (salt_sigma(k,i), k=1,maxlayer)
		end do
	end if
	
	! jw
	if(IS3D_variable(5) == 1) then
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
		do i=1,maxele
			write(pw_tec3D_full,format3) (temp_sigma(k,i), k=1,maxlayer)
		end do
	end if
	
	! jw
	if(IS3D_variable(6) == 1) then
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
! jw
		do i=1,maxele
			write(pw_tec3D_full,format3) (rho_sigma(k,i), k=1,maxlayer)
		end do
	end if

	! jw
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
	write(pw_tec3D_full,'(A,F15.5,A,A,I5)') &
	&	'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, T="TIME = ', &
	&	julian_day * IS3D_time_conv, time_char, '", ZN =', zone_num_3D_full

	write(pw_tec3D_full,*)	! jw
	
	if(IS3D_binary == 1) then
		! jw
	end if
end subroutine write_tec3D_full_body_zto_sigma