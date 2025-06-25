!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_grid_checking_files
	use mod_global_variables
	use mod_file_definition
	
	! local variables ----------------------------------------------------------!
	integer :: i, j, k
	integer :: quotient, remainder
	character(len=40) :: format_string1, format_string2
	! end local variables -----------------------------------------------------!	
	
	! Write grid checking tecplot file (check_grid_2D.dat) ====================!
	! This is the origin shifted version for the nice view.
	if(check_grid_2D == 1) then
		if(check_grid_format == 1 .or. check_grid_format == 3) then	! Tecplot file format
			open(pw_check_grid_2D_tec, file = id_check_grid_2D_tec, form = 'formatted', status = 'replace')	! check_grid_2D.out
			
			write(pw_check_grid_2D_tec,*) 'Title = "Initial 2D grid - origin shifted version"'
			
			! note: do not use:
			!   	if(check_grid_unit_conv == 0.1)
			! 		since, real number annnot be equal to other real number
			if(check_grid_unit_conv > 0.5) then	! 1.0
				write(pw_check_grid_2D_tec,'(A)') 'Variables = "X [m]", "Y [m]", "H [m]", "Bottom Elev. [m]"'
			else											! 0.001
				write(pw_check_grid_2D_tec,'(A)') 'Variables = "X [km]", "Y [km]", "H [m]", "Bottom Elev. [m]"'
			end if
			
			write(pw_check_grid_2D_tec,'(A,F10.2,A,I10,A,I10,A,F10.5,A)') &
			&	'ZONE T = "',0.0,'", N =', MAXNOD, ', E = ', MAXELE, ', DATAPACKING=POINT, &
			&	ZONETYPE=FEQUADRILATERAL, SOLUTIONTIME=', 0.0, ', STRANDID=1' 
		
			do i=1,maxnod
				write(pw_check_grid_2D_tec,'(4F15.5)') 		&
				&	(x_node(i)-xn_min) * check_grid_unit_conv, &	! shift origin to zero
				&	(y_node(i)-yn_min) * check_grid_unit_conv, &	! shift origin to zero
				&	h_node(i), 									&	! positive value is water
				&	-h_node(i)										! bottom elevation from the mean sea level
			end do

			! write cell connectivity
			do i=1,MAXELE
				write(pw_check_grid_2D_tec,'(4I10)') (nodenum_at_cell_tec(k,i),k=1,4)
			end do
			close(pw_check_grid_2D_tec)
		end if
		
		if(check_grid_format == 2 .or. check_grid_format == 3) then	! VTK file format
			open(pw_check_grid_2D_vtk, file = id_check_grid_2D_vtk, form = 'formatted', status = 'replace')	! check_grid_2D.vtk
			write(pw_check_grid_2D_vtk,'(A)') '# vtk DataFile Version 3.0'
			write(pw_check_grid_2D_vtk,'(A)') 'Title = "Initial 2D grid - origin shifted version"'
			write(pw_check_grid_2D_vtk,'(A)') 'ASCII'
			write(pw_check_grid_2D_vtk,*)
			write(pw_check_grid_2D_vtk,'(A)') 'DATASET UNSTRUCTURED_GRID'
			write(pw_check_grid_2D_vtk,'(A,I10,A)') 'POINTS ', MAXNOD, ' float'
			do i=1,MAXNOD
				write(pw_check_grid_2D_vtk,'(F15.5, F15.5, F15.5)') 	&
				&	(x_node(i)-xn_min) * check_grid_unit_conv, 				& 	! shift origin to zero
				&	(y_node(i)-yn_min) * check_grid_unit_conv,					&	! shift origin to zero
				&	0.0
			end do
			write(pw_check_grid_2D_vtk,*)
			
			! write cell connectivity
			write(pw_check_grid_2D_vtk,'(A,I10,I10)') 'CELLS ', MAXELE, cell_connectivity_size_2D
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_2D_vtk,'(I2, 3I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! cell number starts from "0" in VTK
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_2D_vtk,'(I2, 4I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! cell number starts from "0" in VTK
				end if
			end do
			write(pw_check_grid_2D_vtk,*)
			write(pw_check_grid_2D_vtk,'(A,I10)') 'CELL_TYPES ', MAXELE
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_2D_vtk,'(I2)') 5 ! VTK_TRIANGLES
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_2D_vtk,'(I2)') 9 ! VTK_QUAD
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
	! End of check_grid_2D.dat ================================================!

	! Write grid checking tecplot file (check_grid_2DO.dat) ===================!
	! This is the original grid (not origin shifted version)
	if(check_grid_2DO == 1) then
		if(check_grid_format == 1 .or. check_grid_format == 3) then ! Tecplot file format
			open(pw_check_grid_2DO_tec, file = id_check_grid_2DO_tec, form = 'formatted', status = 'replace')	! check_grid_2DO.out
			
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
				&	x_node(i) * check_grid_unit_conv, 			&	! original coordinate
				&	y_node(i) * check_grid_unit_conv, 			&	! original coordinate
				&	h_node(i), 									&	! positive value is water
				&	-h_node(i)										! bottom elevation from the mean sea level
			end do
			
			! write cell connectivity
			do i=1,MAXELE
				write(pw_check_grid_2DO_tec,'(4I10)') (nodenum_at_cell_tec(k,i),k=1,4)
			end do
			close(pw_check_grid_2DO_tec)
		end if
		
		if(check_grid_format == 2 .or. check_grid_format == 3) then ! VTK file format
			open(pw_check_grid_2DO_vtk, file = id_check_grid_2DO_vtk, form = 'formatted', status = 'replace')	! check_grid_2DO.vtk
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
			
			! write cell connectivity
			write(pw_check_grid_2DO_vtk,'(A,I10,I10)') 'CELLS ', MAXELE, cell_connectivity_size_2D
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_2DO_vtk,'(I2, 3I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! cell number starts from "0" in VTK
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_2DO_vtk,'(I2, 4I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! cell number starts from "0" in VTK
				end if
			end do
			write(pw_check_grid_2DO_vtk,*)
			write(pw_check_grid_2DO_vtk,'(A,I10)') 'CELL_TYPES ', MAXELE
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_2DO_vtk,'(I2)') 5 ! VTK_TRIANGLES
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_2DO_vtk,'(I2)') 9 ! VTK_QUAD
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
	! End of check_grid_2DO.dat ===============================================!
	
	
	! Write grid checking tecplot file (check_grid_3D.dat) ====================!
	if(check_grid_3D == 1) then
		if(check_grid_format == 1 .or. check_grid_format == 3) then ! Tecplot file format
			quotient = int(maxnod/100)
			remainder = mod(maxnod,100)
			
			! prepare the format
			write(format_string1,'(A,I0,A)') '(',100,'F15.5)' ! ex, (100F15.5)
			write(format_string2,'(A,I0,A)') '(',remainder,'F15.5)' ! ex, (100F15.5)
			
			open(pw_check_grid_3D_tec, file=id_check_grid_3D_tec, form='formatted', status='replace')
			
			write(pw_check_grid_3D_tec,*) 'Title = "Initial 3D grid checking"'
			
			if(check_grid_unit_conv > 0.5) then
				write(pw_check_grid_3D_tec,'(A)') 'Variables = "X [m]", "Y [m]", "Z [m]"'
			else
				write(pw_check_grid_3D_tec,'(A)') 'Variables = "X [km]", "Y [km]", "Z [m]"'
			end if
				
			write(pw_check_grid_3D_tec,'(A,F10.2,A,I10,A,I10,A,F10.2,A)') 'ZONE T = "',0.0,'", N =', MAXNOD*2, ', E = ', MAXELE, &
			& ', DATAPACKING=BLOCK, ZONETYPE=FEBRICK, SOLUTIONTIME=', 0.0, ', STRANDID=1' 
		
			! write x coordiante ------------------------------------------------!
			do k=1,2
				do i=1,quotient
					! write(pw_check_grid_3D_tec,'(*(F15.5))') ((x_node((i-1)*100+j)-xn_min) * check_grid_unit_conv, j=1,100)
					write(pw_check_grid_3D_tec,format_string1) ((x_node((i-1)*100+j)-xn_min) * check_grid_unit_conv, j=1,100)
				end do
				! write(pw_check_grid_3D_tec,'(*(F15.5))') ((x_node(quotient*100+j)-xn_min) * check_grid_unit_conv, j=1,remainder)
				write(pw_check_grid_3D_tec,format_string2) ((x_node(quotient*100+j)-xn_min) * check_grid_unit_conv, j=1,remainder)
			end do
			write(pw_check_grid_3D_tec,*) ! put a line
				
			! write y coordiante ------------------------------------------------!
			do k=1,2
				do i=1,quotient
					! write(pw_check_grid_3D_tec,'(*(F15.5))') ((y_node((i-1)*100+j)-yn_min) * check_grid_unit_conv, j=1,100)
					write(pw_check_grid_3D_tec,format_string1) ((y_node((i-1)*100+j)-yn_min) * check_grid_unit_conv, j=1,100)
				end do
				! write(pw_check_grid_3D_tec,'(*(F15.5))') ((y_node(quotient*100+j)-yn_min) * check_grid_unit_conv, j=1,remainder)
				write(pw_check_grid_3D_tec,format_string2) ((y_node(quotient*100+j)-yn_min) * check_grid_unit_conv, j=1,remainder)
			end do
			write(pw_check_grid_3D_tec,*) ! put a line
							
			! write bottom elevation and eta at the node points -----------------!
			do i=1,quotient
				write(pw_check_grid_3D_tec,format_string1) (-h_node((i-1)*100+j), j=1,100)	! bottom elevation from MSL
			end do
			write(pw_check_grid_3D_tec,format_string2) (-h_node(quotient*100+j), j=1,remainder)
			write(pw_check_grid_3D_tec,*) ! put a line
			
			! just keed this part as original...
			do i=1,maxnod
				if(h_node(i) > 0.0_dp) then	! located under the MSL
					write(pw_check_grid_3D_tec,'(F15.5)',advance = 'no') 0.0
				else
					write(pw_check_grid_3D_tec,'(F15.5)',advance = 'no') -h_node(i)
				end if
				if(mod(i,100) == 0) then ! change the line every 100 value writing
					write(pw_check_grid_3D_tec,*) ! change the line
				end if
			end do
			write(pw_check_grid_3D_tec,*) ! change the line
			write(pw_check_grid_3D_tec,*) ! put a line
		
			! write cell connectivity ----------------------------------------------!
			do i=1,MAXELE
				write(pw_check_grid_3D_tec,'(8I10)') &
				&	nodenum_at_cell_tec(1,i),nodenum_at_cell_tec(2,i), &
				&	nodenum_at_cell_tec(3,i),nodenum_at_cell_tec(4,i), &
				&	nodenum_at_cell_tec(1,i)+maxnod,nodenum_at_cell_tec(2,i)+maxnod, &
				&	nodenum_at_cell_tec(3,i)+maxnod,nodenum_at_cell_tec(4,i)+maxnod
			end do
		
			! write text information at the bottom of each zone --------------------!
			write(pw_check_grid_3D_tec,'(A,F10.2,A,I5)') 'TEXT CS=FRAME, HU=FRAME, X=50, Y=95, H=2.5, AN=MIDCENTER, &
			&	T="TIME = ',0.0,' sec", ZN =', 1
			write(pw_check_grid_3D_tec,*)	! put a line			
			close(pw_check_grid_3D_tec)
		end if
		
		if(check_grid_format == 2 .or. check_grid_format == 3) then ! VTK file format
			open(pw_check_grid_3D_vtk, file = id_check_grid_3D_vtk, form = 'formatted', status = 'replace')	! check_grid_3D.vtk
			write(pw_check_grid_3D_vtk,'(A)') '# vtk DataFile Version 3.0'
			write(pw_check_grid_3D_vtk,'(A)') 'Title = "Initial 3D grid checking"'
			write(pw_check_grid_3D_vtk,'(A)') 'ASCII'
			write(pw_check_grid_3D_vtk,*)
			write(pw_check_grid_3D_vtk,'(A)') 'DATASET UNSTRUCTURED_GRID'
			write(pw_check_grid_3D_vtk,'(A,I10,A)') 'POINTS ', MAXNOD*2, ' float'
			! bottom grids
			do i=1,MAXNOD
				write(pw_check_grid_3D_vtk,'(F15.5, F15.5, F15.5)') 	&
				&	(x_node(i)-xn_min) * check_grid_unit_conv, 				& 	! shift origin to zero
				&	(y_node(i)-yn_min) * check_grid_unit_conv,					&	! shift origin to zero
				&	-h_node(i)														! bottom elevation from MSL
			end do
			! surface grids
			do i=1,MAXNOD
				if(h_node(i) > 0.0) then ! wet cells
					write(pw_check_grid_3D_vtk,'(F15.5, F15.5, F15.5)') 	&
					&	(x_node(i)-xn_min) * check_grid_unit_conv, 				& 	! shift origin to zero
					&	(y_node(i)-yn_min) * check_grid_unit_conv,					&	! shift origin to zero
					&	0.0																	! set surface at MSL
				else
					write(pw_check_grid_3D_vtk,'(F15.5, F15.5, F15.5)') 	&
					&	(x_node(i)-xn_min) * check_grid_unit_conv, 				& 	! shift origin to zero
					&	(y_node(i)-yn_min) * check_grid_unit_conv,					&	! shift origin to zero
					&	-h_node(i)														! set surface at surface elevation
				end if
			end do
			
			write(pw_check_grid_3D_vtk,*)
			! write cell connectivity
			write(pw_check_grid_3D_vtk,'(A,I10,I10)') 'CELLS ', MAXELE, cell_connectivity_size_3D
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_3D_vtk,'(I2, 6I10)') &
					&	tri_or_quad(i)*2, (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)), &
					&	(nodenum_at_cell(j,i)-1+maxnod, j=1,tri_or_quad(i)) ! cell number starts from "0" in VTK
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_3D_vtk,'(I2, 8I10)') &
					&	tri_or_quad(i)*2, (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)), &
					&	(nodenum_at_cell(j,i)-1+maxnod, j=1,tri_or_quad(i)) ! cell number starts from "0" in VTK
				end if
			end do
			write(pw_check_grid_3D_vtk,*)
			write(pw_check_grid_3D_vtk,'(A,I10)') 'CELL_TYPES ', MAXELE
			do i=1,MAXELE
				if(tri_or_quad(i) == 3) then
					write(pw_check_grid_3D_vtk,'(I2)') 13 ! VTK_WEDGE
				else if(tri_or_quad(i) == 4) then
					write(pw_check_grid_3D_vtk,'(I2)') 12 ! VTK_HEXAHEDRON
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
	! End of check_grid_3D.dat ================================================!
	
	
	! Write additional grid information file: check_grid.dat ==================!	
	! not yet finished... copy this part from REM2D
	if(check_grid_info == 1) then
		! not yet finished
	end if
end subroutine write_grid_checking_files
