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
	
	! new values for sigma coordinates:
	real(dp),dimension(0:maxlayer,maxnod) :: u_sigma, v_sigma, w_sigma 				! these will be defined at vertical levels
	real(dp),dimension(maxlayer,maxele) :: salt_sigma, temp_sigma, rho_sigma 	! these will be defined at vertical layers
	! End of local variables ==================================================!

	! Data preparation at vertical levels or layers ===========================!	
	! Note:
	! 		u,v,w will be defined at nodes at sigma levels
	! 		salt,temp,density will be defined at cell centers at sigma layer
		
	! interpolate nodal values into sigma levels
	call zto_sigma_node_values(u_sigma, v_sigma, w_sigma)

	! interpolate cell values into sigma layer
	call zto_sigma_cell_values(salt_sigma, temp_sigma, rho_sigma)
	! End of data preparation at vertical levels or layers ====================!

	! Now, let's write tecplot file ===========================================!
	! prepare formats before starts:
	quotient_node = int(maxnod/100)
	remainder_node = mod(maxnod,100)
	quotient_cell = int(maxele/100)
	remainder_cell = mod(maxele,100)

	! prepare the format
! 	write(format1,'(A,I0,A)') '(',100,'E15.5E4)'
! 	write(format2_cell,'(A,I0,A)') '(',remainder_cell,'E15.5E4)'
! 	write(format2_node,'(A,I0,A)') '(',remainder_node,'E15.5E4)'

! 	write(format1,'(A,I0,A)') '(',100,'F15.5)'
! 	write(format2_cell,'(A,I0,A)') '(',remainder_cell,'F15.5)'
! 	write(format2_node,'(A,I0,A)') '(',remainder_node,'F15.5)'
! 	write(format3,'(A,I0,A)') '(',maxlayer,'F15.5)'

	write(format1,'(A,I0,A)') '(',100,'(F0.5,1x))'
	write(format2_cell,'(A,I0,A)') '(',remainder_cell,'(F0.5,1x))'
	write(format2_node,'(A,I0,A)') '(',remainder_node,'(F0.5,1x))'
	write(format3,'(A,I0,A)') '(',maxlayer,'(F0.5,1x))'
	
	! I will use "F" format for x,y,z coordinate and "E" format for other variables.
	! This will allow to show extream values even though we have extreme (wrong) values.
	
	! Write header lines ======================================================!
	! Note: the starting line is different from write_tec3D_full_head_zto_sigma.f90
	write(pw_tec3D_full,'(A,F15.5, A,I10, A,I10, A)', advance = 'no') &
	&	'ZONE T = "', julian_day * IS3D_time_conv,	&
	&	'", N =', MAXNOD*(maxlayer+1), &
	&	', E = ', MAXELE*maxlayer, &
	&	', DATAPACKING=BLOCK'

	! let's count nodal variables (i.e., u,v,w)
	count1 = 0
	do i=1,3
		if(IS3D_variable(i) == 1) then
			count1 = count1+1
		end if
	end do
	
	! then, this is for cell center variables (i.e., salt, temp, rho)
	! thus, if you want to add additional variables in the future, locate nodal variables first.
	! first, I have to count how many cell centered values are asigned.
	! if there is no cell centered values selected, 'VARLOCATION=' term should not be included.
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
			! here (i+3) is for (u,v,w), 
			! thus, if you add additional nodal variables in the future, increase this number, "3"
			if(IS3D_variable(i+3) == 1) then
				count3 = count3 + 1
				write(pw_tec3D_full,'(I1,A)',advance = 'no') (3+count1)+count3, ',' ! here 3 is for (x,y,z), count1 is for u,v,w
			end if
			! note: the last comma will not be a problem, so ignore the last comma here
		end do
		write(pw_tec3D_full,'(A)',advance = 'no') &
		&	']=CELLCENTERED)'
	end if
	
	write(pw_tec3D_full,'(A,A, A,F15.5, A,F15.5, A)') &
	&	', ZONETYPE=',zonetype, &
	&	', SOLUTIONTIME=', julian_day * IS3D_time_conv, &
	&	', VARSHARELIST=([1-2]=1), CONNECTIVITYSHAREZONE=1, SOLUTIONTIME=', julian_day * IS3D_time_conv, & ! note: this is different from write_tec3D_full_head_zto_sigma.f90
	&	', STRANDID=1' 
	
	! skip writing x & y coordinates but z coordinate should be updated =======!
	! write node elevation at each vertical level: ----------------------------!
	! I am creating bricks from bottom to top, thus values also should be written from bottom to top order.
	do k=0,maxlayer ! write (maxlayer + 1) times
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
	write(pw_tec3D_full,*)	! put a line		

	! write u at node ---------------------------------------------------------!
	if(IS3D_variable(1) == 1) then
		do k=0,maxlayer
			do i=1,quotient_node
				write(pw_tec3D_full,format1) (u_sigma(k,(i-1)*100+j), j=1,100)
			end do
			if(remainder_node > 0) then
				write(pw_tec3D_full,format2_node) (u_sigma(k,quotient_node*100+j), j=1,remainder_node)
			end if
		end do
		write(pw_tec3D_full,*)	! put a line
	end if

	! write v at node ---------------------------------------------------------!
	if(IS3D_variable(2) == 1) then
		do k=0,maxlayer
			do i=1,quotient_node
				write(pw_tec3D_full,format1) (v_sigma(k,(i-1)*100+j), j=1,100)
			end do
			if(remainder_node > 0) then
				write(pw_tec3D_full,format2_node) (v_sigma(k,quotient_node*100+j), j=1,remainder_node)
			end if
		end do
		write(pw_tec3D_full,*)	! put a line
	end if

	! write w at node ---------------------------------------------------------!
	if(IS3D_variable(3) == 1) then
		do k=0,maxlayer
			do i=1,quotient_node
				write(pw_tec3D_full,format1) (w_sigma(k,(i-1)*100+j), j=1,100)
			end do
			if(remainder_node > 0) then
				write(pw_tec3D_full,format2_node) (w_sigma(k,quotient_node*100+j), j=1,remainder_node)
			end if
		end do
		write(pw_tec3D_full,*)	! put a line
	end if
	
	! write salinity at cell center -------------------------------------------!
	if(IS3D_variable(4) == 1) then
! 		do k=1,maxlayer
! 			do i=1,quotient_cell
! 				write(pw_tec3D_full,format1) (salt_sigma(k,(i-1)*100+j), j=1,100)
! 			end do
! 			if(remainder_cell > 0) then
! 				write(pw_tec3D_full,format2_cell) (salt_sigma(k,quotient_cell*100+j), j=1,remainder_cell)
! 			end if
! 		end do
! 		write(pw_tec3D_full,*)	! put a line
		do i=1,maxele
			write(pw_tec3D_full,format3) (salt_sigma(k,i), k=1,maxlayer)
		end do
	end if
	
	! write temperature at cell center ----------------------------------------!
	if(IS3D_variable(5) == 1) then
! 		do k=1,maxlayer
! 			do i=1,quotient_cell
! 				write(pw_tec3D_full,format1) (temp_sigma(k,(i-1)*100+j), j=1,100)
! 			end do
! 			if(remainder_cell > 0) then
! 				write(pw_tec3D_full,format2_cell) (temp_sigma(k,quotient_cell*100+j), j=1,remainder_cell)
! 			end if
! 		end do
! 		write(pw_tec3D_full,*)	! put a line
		do i=1,maxele
			write(pw_tec3D_full,format3) (temp_sigma(k,i), k=1,maxlayer)
		end do
	end if
	
	! write density (rho_cell) at cell center ---------------------------------!
	if(IS3D_variable(6) == 1) then
! 		do k=1,maxlayer
! 			do i=1,quotient_cell
! 				write(pw_tec3D_full,format1) (rho_sigma(k,(i-1)*100+j), j=1,100)
! 			end do
! 			if(remainder_cell > 0) then
! 				write(pw_tec3D_full,format2_cell) (rho_sigma(k,quotient_cell*100+j), j=1,remainder_cell)
! 			end if
! 		end do
! 		write(pw_tec3D_full,*)	! put a line
		do i=1,maxele
			write(pw_tec3D_full,format3) (rho_sigma(k,i), k=1,maxlayer)
		end do
	end if

	! skip writing cell connectivity ==========================================!
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
	write(pw_tec3D_full,'(A,F15.5,A,A,I5)') &
	&	'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, T="TIME = ', &
	&	julian_day * IS3D_time_conv, time_char, '", ZN =', zone_num_3D_full

	write(pw_tec3D_full,*)	! put a line
	
	if(IS3D_binary == 1) then
		! not yet activated
	end if
end subroutine write_tec3D_full_body_zto_sigma