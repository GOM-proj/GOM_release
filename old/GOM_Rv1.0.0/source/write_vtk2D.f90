!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_vtk2D
	use mod_global_variables
	use mod_file_definition	
	implicit none
	
	integer :: i, j
	character(len=40) :: format_string
	character(len= 4) :: File_num_buff
	! jw

	format_string = '(I4.4)'	! jw

	write(File_num_buff,format_string) IS2D_vtk_num		
	IS2D_File_name = trim(id_vtk2D)//trim(File_num_buff)//'.vtk'
	open(pw_vtk2D, file = trim(IS2D_File_name), form = 'formatted', status = 'replace')	! jw
	
	write(pw_vtk2D,'(A)') '# vtk DataFile Version 3.0'
	write(pw_vtk2D,'(A)') 'Title = "2D contour plot"'
	write(pw_vtk2D,'(A)') 'ASCII'
	write(pw_vtk2D,*)
	write(pw_vtk2D,'(A)') 'DATASET UNSTRUCTURED_GRID'
	
	! jw
	write(pw_vtk2D,'(A)') 'FIELD FieldData 1'
	write(pw_vtk2D,'(A)') 'TIME 1 1 float'
	write(pw_vtk2D,'(F15.5)') julian_day * IS2D_time_conv
	
	write(pw_vtk2D,'(A,I10,A)') 'POINTS ', MAXNOD, ' float'
	do i=1,MAXNOD
		write(pw_vtk2D,'(F15.5, F15.5, F15.5)')	&
		&	(x_node(i)-xn_min) * IS2D_unit_conv, 	& 	! jw
		&	(y_node(i)-yn_min) * IS2D_unit_conv,	&	! jw
		&	0.0
	end do
	write(pw_vtk2D,*)
		
	! jw
	write(pw_vtk2D,'(A,I10,I10)') 'CELLS ', MAXELE, cell_connectivity_size_2D
	do i=1,MAXELE
		if(tri_or_quad(i) == 3) then
			write(pw_vtk2D,'(I2, 3I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! jw
		else if(tri_or_quad(i) == 4) then
			write(pw_vtk2D,'(I2, 4I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! jw
		end if
	end do
	write(pw_vtk2D,*)
	write(pw_vtk2D,'(A,I10)') 'CELL_TYPES ', MAXELE
	do i=1,MAXELE
		if(tri_or_quad(i) == 3) then
			write(pw_vtk2D,'(I2)') 5 ! jw
		else if(tri_or_quad(i) == 4) then
			write(pw_vtk2D,'(I2)') 9 ! jw
		end if
	end do
	
	! jw
	! jw
	if(maxval(IS2D_variable) > 0) then
		write(pw_vtk2D,*)
		write(pw_vtk2D,'(A,I10)') 'POINT_DATA ', MAXNOD
	end if

	if(IS2D_variable(1) == 1) then
		write(pw_vtk2D,*)
		write(pw_vtk2D,'(A)') 'SCALARS Eta(m) float'
		write(pw_vtk2D,'(A)') 'LOOKUP_TABLE default'
		do i=1,maxnod
			write(pw_vtk2D,'(F15.5)') eta_node(i)
		end do
	end if
	if(IS2D_variable(2) == 1 .or. IS2D_variable(3) == 1) then	
		write(pw_vtk2D,*)
		write(pw_vtk2D,'(A)') 'VECTORS Velocity(m/s) float'
		do i=1,maxnod
			write(pw_vtk2D,'(F15.5 F15.5 F15.5)') &
			&	ubar_node(i), vbar_node(i), 0.0
		end do
	end if
	if(IS2D_variable(4) == 1) then
		write(pw_vtk2D,*)
		write(pw_vtk2D,'(A)') 'SCALARS Salt(psu) float'
		write(pw_vtk2D,'(A)') 'LOOKUP_TABLE default'
		do i=1,maxnod
			write(pw_vtk2D,'(F15.5)') sbar_node(i)
		end do			
	end if
	if(IS2D_variable(5) == 1) then
		write(pw_vtk2D,*)
		write(pw_vtk2D,'(A)') 'SCALARS Temp(C) float'
		write(pw_vtk2D,'(A)') 'LOOKUP_TABLE default'
		do i=1,maxnod
			write(pw_vtk2D,'(F15.5)') tbar_node(i)
		end do		
	end if
	if(IS2D_variable(6) == 1) then
		write(pw_vtk2D,*)
		write(pw_vtk2D,'(A)') 'SCALARS Rho(kg/m3) float'
		write(pw_vtk2D,'(A)') 'LOOKUP_TABLE default'
		do i=1,maxnod
			write(pw_vtk2D,'(F15.5)') rbar_node(i)
		end do		
	end if		

	close(pw_vtk2D)
end subroutine write_vtk2D