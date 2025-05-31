!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_flood_map_vtk
	use mod_global_variables
	use mod_file_definition	
	implicit none
	
	integer :: i, j
	! jw

	write(pw_flood_map_vtk,'(A)') '# vtk DataFile Version 3.0'
	write(pw_flood_map_vtk,'(A)') 'Title = "2D flood map contour plot"'
	write(pw_flood_map_vtk,'(A)') 'ASCII'
	write(pw_flood_map_vtk,*)
	write(pw_flood_map_vtk,'(A)') 'DATASET UNSTRUCTURED_GRID'
	
	! jw
	write(pw_flood_map_vtk,'(A)') 'FIELD FieldData 1'
	write(pw_flood_map_vtk,'(A)') 'TIME 1 1 float'
	write(pw_flood_map_vtk,'(F15.5)') elapsed_time * IS2D_time_conv
	
	write(pw_flood_map_vtk,'(A,I10,A)') 'POINTS ', MAXNOD, ' float'
	do i=1,MAXNOD
		write(pw_flood_map_vtk,'(F15.5, F15.5, F15.5)') 			&
		&	x_node(i) * IS2D_unit_conv, 	&
		&	y_node(i) * IS2D_unit_conv,	&
		&	0.0
	end do
	write(pw_flood_map_vtk,*)
		
	! jw
	write(pw_flood_map_vtk,'(A,I10,I10)') 'CELLS ', MAXELE, cell_connectivity_size_2D
	do i=1,MAXELE
		if(tri_or_quad(i) == 3) then
			write(pw_flood_map_vtk,'(I2, 3I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! jw
		else if(tri_or_quad(i) == 4) then
			write(pw_flood_map_vtk,'(I2, 4I10)') tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! jw
		end if
	end do
	write(pw_flood_map_vtk,*)
	write(pw_flood_map_vtk,'(A,I10)') 'CELL_TYPES ', MAXELE
	do i=1,MAXELE
		if(tri_or_quad(i) == 3) then
			write(pw_flood_map_vtk,'(I2)') 5 ! jw
		else if(tri_or_quad(i) == 4) then
			write(pw_flood_map_vtk,'(I2)') 9 ! jw
		end if
	end do
	write(pw_flood_map_vtk,*)
	write(pw_flood_map_vtk,'(A,I10)') 'POINT_DATA ', MAXNOD
	write(pw_flood_map_vtk,'(A)') 'SCALARS max_eta(m) float'
	write(pw_flood_map_vtk,'(A)') 'LOOKUP_TABLE default'
	do i=1,maxnod
		write(pw_flood_map_vtk,'(F15.5)') max_eta_node(i)
	end do
	write(pw_flood_map_vtk,*)
	write(pw_flood_map_vtk,'(A)') 'SCALARS max_flood_time float'
	write(pw_flood_map_vtk,'(A)') 'LOOKUP_TABLE default'
	do i=1,maxnod
		write(pw_flood_map_vtk,'(F15.5)') max_flood_time(i)
	end do
	write(pw_flood_map_vtk,*)
	write(pw_flood_map_vtk,'(A)') 'SCALARS flood_id integer'
	write(pw_flood_map_vtk,'(A)') 'LOOKUP_TABLE default'
	do i=1,maxnod
		write(pw_flood_map_vtk,'(I3)') flood_id(i)
	end do
	
	close(pw_flood_map_vtk)
end subroutine write_flood_map_vtk