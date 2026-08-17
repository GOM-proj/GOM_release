!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_grid_checking_files
	use mod_global_variables
	use mod_file_definition
	
	! jw
	integer :: i, j, k
	integer :: quotient, remainder
	character(len=40) :: format_string1, format_string2
	! jw
	
	! jw
	! jw
	if(check_grid_2D == 1) then
		if(check_grid_format == 1 .or. check_grid_format == 3) then	! jw
			open(pw_check_grid_2D_tec, file = id_check_grid_2D_tec, form = 'formatted', status = 'replace')	! jw
			
			write(pw_check_grid_2D_tec,*) 'Title = "Initial 2D grid - origin shifted version"'
			
			! jw
			! jw
			! jw
			if(check_grid_unit_conv > 0.5) then	! jw
				write(pw_check_grid_2D_tec,'(A)') 'Variables = "X [m]", "Y [m]", "H [m]", "Bottom Elev. [m]"'
			else											! jw
				write(pw_check_grid_2D_tec,'(A)') 'Variables = "X [km]", "Y [km]", "H [m]", "Bottom Elev. [m]"'
			end if
			
			write(pw_check_grid_2D_tec,'(A,F10.2,A,I10,A,I10,A,F10.5,A)') &
			&	'ZONE T = "',0.0,'", N =', MAXNOD, ', E = ', MAXELE, ', DATAPACKING=POINT, &
			&	ZONETYPE=FEQUADRILATERAL, SOLUTIONTIME=', 0.0, ', STRANDID=1' 
		
			do i=1,maxnod
				write(pw_check_grid_2D_tec,'(4F15.5)') 		&
				&	(x_node(i)-xn_min) * check_grid_unit_conv, &	! jw
				&	(y_node(i)-yn_min) * check_grid_unit_conv, &	! jw
				&	h_node(i), 									&	! jw
				&	-h_node(i)										! jw
			end do

			! jw
			do i=1,MAXELE
				write(pw_check_grid_2D_tec,'(4I10)') (nodenum_at_cell_tec(k,i),k=1,4)
			end do
			close(pw_check_grid_2D_tec)
		end if
		
		if(check_grid_format == 2 .or. check_grid_format == 3) then	! jw
			open(pw_check_grid_2D_vtk, file = id_check_grid_2D_vtk, form = 'formatted', status = 'replace')	! jw
			write(pw_check_grid_2D_vtk,'(A)') '# vtk DataFile Version 3.0'
			write(pw_check_grid_2D_vtk,'(A)') 'Title = "Initial 2D grid - origin shifted version"'
			write(pw_check_grid_2D_vtk,'(A)') 'ASCII'
			write(pw_check_grid_2D_vtk,*)
			write(pw_check_grid_2D_vtk,'(A)') 'DATASET UNSTRUCTURED_GRID'
			write(pw_check_grid_2D_vtk,'(A,I10,A)') 'POINTS ', MAXNOD, ' float'
			do i=1,MAXNOD
				write(pw_check_grid_2D_vtk,'(F15.5, F15.5, F15.5)') 	&
				&	(x_node(i)-xn_min) * check_grid_unit_conv, 				& 	! jw
				&	(y_node(i)-yn_min) * check_grid_unit_conv,					&	! jw
				&	0.0
			end do
			write(pw_check_grid_2D_vtk,*)
			
			! jw
			write(pw_check_grid_2D_vtk,'(A,I10,I10)') 'CELLS ', MAXELE, cell_connectivity_size_2D
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_2D_vtk,'(I2, 3I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! jw
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_2D_vtk,'(I2, 4I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! jw
				end if
			end do
			write(pw_check_grid_2D_vtk,*)
			write(pw_check_grid_2D_vtk,'(A,I10)') 'CELL_TYPES ', MAXELE
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_2D_vtk,'(I2)') 5 ! jw
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_2D_vtk,'(I2)') 9 ! jw
				end if
			end do
			write(pw_check_grid_2D_vtk,*)
			write(pw_check_grid_2D_vtk,'(A,I10)') 'POINT_DATA ', MAXNOD
			write(pw_check_grid_2D_vtk,'(A)') 'SCALARS Depth(m) float'
			write(pw_check_grid_2D_vtk,'(A)') 'LOOKUP_TABLE default'
			do i=1,MAXNOD
				write(pw_check_grid_2D_vtk,'(F15.5)') h_node(i)
			end do
			write(pw_check_grid_2D_vtk,*)
			write(pw_check_grid_2D_vtk,'(A)') 'SCALARS Bottom_Elev(m) float'
			write(pw_check_grid_2D_vtk,'(A)') 'LOOKUP_TABLE default'
			do i=1,MAXNOD
				write(pw_check_grid_2D_vtk,'(F15.5)') -h_node(i)
			end do
			
			close(pw_check_grid_2D_vtk)
		end if
	end if
	! jw

	! jw
	! jw
	if(check_grid_2DO == 1) then
		if(check_grid_format == 1 .or. check_grid_format == 3) then ! jw
			open(pw_check_grid_2DO_tec, file = id_check_grid_2DO_tec, form = 'formatted', status = 'replace')	! jw
			
			write(pw_check_grid_2DO_tec,'(A)') 'Title = "Initial 2D grid - original version"'	
			
			if(check_grid_unit_conv > 0.5) then
				write(pw_check_grid_2DO_tec,'(A)') 'Variables = "X [m]", "Y [m]", "H [m]", "Bottom Elev. [m]"'
			else
				write(pw_check_grid_2DO_tec,'(A)') 'Variables = "X [km]", "Y [km]", "H [m]", "Bottom Elev. [m]"'
			end if
			
			write(pw_check_grid_2DO_tec,'(A,F10.2,A,I10,A,I10,A,F10.5,A)') &
				& 'ZONE T = "',0.0,'", N =', MAXNOD, ', E = ', MAXELE, ', DATAPACKING=POINT, &
				&	ZONETYPE=FEQUADRILATERAL, SOLUTIONTIME=', 0.0, ', STRANDID=1' 
		
			do i=1,maxnod
				write(pw_check_grid_2DO_tec,'(4F15.5)') 		&
				&	x_node(i) * check_grid_unit_conv, 			&	! jw
				&	y_node(i) * check_grid_unit_conv, 			&	! jw
				&	h_node(i), 									&	! jw
				&	-h_node(i)										! jw
			end do
			
			! jw
			do i=1,MAXELE
				write(pw_check_grid_2DO_tec,'(4I10)') (nodenum_at_cell_tec(k,i),k=1,4)
			end do
			close(pw_check_grid_2DO_tec)
		end if
		
		if(check_grid_format == 2 .or. check_grid_format == 3) then ! jw
			open(pw_check_grid_2DO_vtk, file = id_check_grid_2DO_vtk, form = 'formatted', status = 'replace')	! jw
			write(pw_check_grid_2DO_vtk,'(A)') '# vtk DataFile Version 3.0'
			write(pw_check_grid_2DO_vtk,'(A)') 'Title = "Initial 2D grid - original version"'	
			write(pw_check_grid_2DO_vtk,'(A)') 'ASCII'
			write(pw_check_grid_2DO_vtk,*)
			write(pw_check_grid_2DO_vtk,'(A)') 'DATASET UNSTRUCTURED_GRID'
			write(pw_check_grid_2DO_vtk,'(A,I10,A)') 'POINTS ', MAXNOD, ' float'
			do i=1,MAXNOD
				write(pw_check_grid_2DO_vtk,'(F15.5, F15.5, F15.5)') 	&
				&	x_node(i) * check_grid_unit_conv, 							&
				&	y_node(i) * check_grid_unit_conv,								&
				&	0.0
			end do
			write(pw_check_grid_2DO_vtk,*)
			
			! jw
			write(pw_check_grid_2DO_vtk,'(A,I10,I10)') 'CELLS ', MAXELE, cell_connectivity_size_2D
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_2DO_vtk,'(I2, 3I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! jw
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_2DO_vtk,'(I2, 4I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! jw
				end if
			end do
			write(pw_check_grid_2DO_vtk,*)
			write(pw_check_grid_2DO_vtk,'(A,I10)') 'CELL_TYPES ', MAXELE
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_2DO_vtk,'(I2)') 5 ! jw
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_2DO_vtk,'(I2)') 9 ! jw
				end if
			end do
			write(pw_check_grid_2DO_vtk,*)
			write(pw_check_grid_2DO_vtk,'(A,I10)') 'POINT_DATA ', MAXNOD
			write(pw_check_grid_2DO_vtk,'(A)') 'SCALARS Depth(m) float'
			write(pw_check_grid_2DO_vtk,'(A)') 'LOOKUP_TABLE default'
			do i=1,MAXNOD
				write(pw_check_grid_2DO_vtk,'(F15.5)') h_node(i)
			end do
			write(pw_check_grid_2DO_vtk,*)
			write(pw_check_grid_2DO_vtk,'(A)') 'SCALARS Bottom_Elev(m) float'
			write(pw_check_grid_2DO_vtk,'(A)') 'LOOKUP_TABLE default'
			do i=1,MAXNOD
				write(pw_check_grid_2DO_vtk,'(F15.5)') -h_node(i)
			end do
			
			close(pw_check_grid_2DO_vtk)
		end if
	end if
	! jw
	
	
	! jw
	if(check_grid_3D == 1) then
		if(check_grid_format == 1 .or. check_grid_format == 3) then ! jw
			quotient = int(maxnod/100)
			remainder = mod(maxnod,100)
			
			! jw
			write(format_string1,'(A,I0,A)') '(',100,'F15.5)' ! jw
			write(format_string2,'(A,I0,A)') '(',remainder,'F15.5)' ! jw
			
			open(pw_check_grid_3D_tec, file=id_check_grid_3D_tec, form='formatted', status='replace')
			
			write(pw_check_grid_3D_tec,*) 'Title = "Initial 3D grid checking"'
			
			if(check_grid_unit_conv > 0.5) then
				write(pw_check_grid_3D_tec,'(A)') 'Variables = "X [m]", "Y [m]", "Z [m]"'
			else
				write(pw_check_grid_3D_tec,'(A)') 'Variables = "X [km]", "Y [km]", "Z [m]"'
			end if
				
			write(pw_check_grid_3D_tec,'(A,F10.2,A,I10,A,I10,A,F10.2,A)') 'ZONE T = "',0.0,'", N =', MAXNOD*2, ', E = ', MAXELE, &
			& ', DATAPACKING=BLOCK, ZONETYPE=FEBRICK, SOLUTIONTIME=', 0.0, ', STRANDID=1' 
		
			! jw
			do k=1,2
				do i=1,quotient
					! jw
					write(pw_check_grid_3D_tec,format_string1) ((x_node((i-1)*100+j)-xn_min) * check_grid_unit_conv, j=1,100)
				end do
				! jw
				write(pw_check_grid_3D_tec,format_string2) ((x_node(quotient*100+j)-xn_min) * check_grid_unit_conv, j=1,remainder)
			end do
			write(pw_check_grid_3D_tec,*) ! jw
				
			! jw
			do k=1,2
				do i=1,quotient
					! jw
					write(pw_check_grid_3D_tec,format_string1) ((y_node((i-1)*100+j)-yn_min) * check_grid_unit_conv, j=1,100)
				end do
				! jw
				write(pw_check_grid_3D_tec,format_string2) ((y_node(quotient*100+j)-yn_min) * check_grid_unit_conv, j=1,remainder)
			end do
			write(pw_check_grid_3D_tec,*) ! jw
							
			! jw
			do i=1,quotient
				write(pw_check_grid_3D_tec,format_string1) (-h_node((i-1)*100+j), j=1,100)	! jw
			end do
			write(pw_check_grid_3D_tec,format_string2) (-h_node(quotient*100+j), j=1,remainder)
			write(pw_check_grid_3D_tec,*) ! jw
			
			! jw
			do i=1,maxnod
				if(h_node(i) > 0.0_dp) then	! jw
					write(pw_check_grid_3D_tec,'(F15.5)',advance = 'no') 0.0
				else
					write(pw_check_grid_3D_tec,'(F15.5)',advance = 'no') -h_node(i)
				end if
				if(mod(i,100) == 0) then ! jw
					write(pw_check_grid_3D_tec,*) ! jw
				end if
			end do
			write(pw_check_grid_3D_tec,*) ! jw
			write(pw_check_grid_3D_tec,*) ! jw
		
			! jw
			do i=1,MAXELE
				write(pw_check_grid_3D_tec,'(8I10)') &
				&	nodenum_at_cell_tec(1,i),nodenum_at_cell_tec(2,i), &
				&	nodenum_at_cell_tec(3,i),nodenum_at_cell_tec(4,i), &
				&	nodenum_at_cell_tec(1,i)+maxnod,nodenum_at_cell_tec(2,i)+maxnod, &
				&	nodenum_at_cell_tec(3,i)+maxnod,nodenum_at_cell_tec(4,i)+maxnod
			end do
		
			! jw
			write(pw_check_grid_3D_tec,'(A,F10.2,A,I5)') 'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, &
			&	T="TIME = ',0.0,' sec", ZN =', 1
			write(pw_check_grid_3D_tec,*)	! jw
			close(pw_check_grid_3D_tec)
		end if
		
		if(check_grid_format == 2 .or. check_grid_format == 3) then ! jw
			open(pw_check_grid_3D_vtk, file = id_check_grid_3D_vtk, form = 'formatted', status = 'replace')	! jw
			write(pw_check_grid_3D_vtk,'(A)') '# vtk DataFile Version 3.0'
			write(pw_check_grid_3D_vtk,'(A)') 'Title = "Initial 3D grid checking"'
			write(pw_check_grid_3D_vtk,'(A)') 'ASCII'
			write(pw_check_grid_3D_vtk,*)
			write(pw_check_grid_3D_vtk,'(A)') 'DATASET UNSTRUCTURED_GRID'
			write(pw_check_grid_3D_vtk,'(A,I10,A)') 'POINTS ', MAXNOD*2, ' float'
			! jw
			do i=1,MAXNOD
				write(pw_check_grid_3D_vtk,'(F15.5, F15.5, F15.5)') 	&
				&	(x_node(i)-xn_min) * check_grid_unit_conv, 				& 	! jw
				&	(y_node(i)-yn_min) * check_grid_unit_conv,					&	! jw
				&	-h_node(i)														! jw
			end do
			! jw
			do i=1,MAXNOD
				if(h_node(i) > 0.0) then ! jw
					write(pw_check_grid_3D_vtk,'(F15.5, F15.5, F15.5)') 	&
					&	(x_node(i)-xn_min) * check_grid_unit_conv, 				& 	! jw
					&	(y_node(i)-yn_min) * check_grid_unit_conv,					&	! jw
					&	0.0																	! jw
				else
					write(pw_check_grid_3D_vtk,'(F15.5, F15.5, F15.5)') 	&
					&	(x_node(i)-xn_min) * check_grid_unit_conv, 				& 	! jw
					&	(y_node(i)-yn_min) * check_grid_unit_conv,					&	! jw
					&	-h_node(i)														! jw
				end if
			end do
			
			write(pw_check_grid_3D_vtk,*)
			! jw
			write(pw_check_grid_3D_vtk,'(A,I10,I10)') 'CELLS ', MAXELE, cell_connectivity_size_3D
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_3D_vtk,'(I2, 6I10)') &
					&	tri_or_quad(i)*2, (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)), &
					&	(nodenum_at_cell(j,i)-1+maxnod, j=1,tri_or_quad(i)) ! jw
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_3D_vtk,'(I2, 8I10)') &
					&	tri_or_quad(i)*2, (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)), &
					&	(nodenum_at_cell(j,i)-1+maxnod, j=1,tri_or_quad(i)) ! jw
				end if
			end do
			write(pw_check_grid_3D_vtk,*)
			write(pw_check_grid_3D_vtk,'(A,I10)') 'CELL_TYPES ', MAXELE
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_3D_vtk,'(I2)') 13 ! jw
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_3D_vtk,'(I2)') 12 ! jw
				end if
			end do
			write(pw_check_grid_3D_vtk,*)
			write(pw_check_grid_3D_vtk,'(A,I10)') 'CELL_DATA ', MAXELE
			write(pw_check_grid_3D_vtk,'(A)') 'SCALARS Depth(m) float'
			write(pw_check_grid_3D_vtk,'(A)') 'LOOKUP_TABLE default'
			do i=1,MAXELE
				write(pw_check_grid_3D_vtk,'(F15.5)') h_cell(i)
			end do
			
			close(pw_check_grid_3D_vtk)			
		end if
	end if
	! jw
	
	
	! jw
	! jw
	if(check_grid_info == 1) then
		! jw
	end if
end subroutine write_grid_checking_files
