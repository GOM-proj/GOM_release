!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine write_vtk3D
	use mod_global_variables
	use mod_file_definition	
	implicit none
	
	integer :: i, j, k
	character(len=40) :: format_string
	character(len= 4) :: File_num_buff
	
	character(len=10),dimension(maxele) :: format_cell_connectivity_2D, format_cell_connectivity_3D
	integer,dimension(maxele) :: cell_type_2D, cell_type_3D
	
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
	format_string = '(I4.4)'	! jw
	write(File_num_buff,format_string) IS3D_vtk_num		
	
	do i=1,maxele
		if(tri_or_quad(i) == 3) then
			format_cell_connectivity_2D(i) = '(I2, 3I10)'
			format_cell_connectivity_3D(i) = '(I2, 6I10)'
			cell_type_2D(i) = 5 ! jw
			cell_type_3D(i) = 13 ! jw
		else if(tri_or_quad(i) == 4) then
			format_cell_connectivity_2D(i) = '(I2, 4I10)'
			format_cell_connectivity_3D(i) = '(I2, 8I10)'
			cell_type_2D(i) = 9 ! jw
			cell_type_3D(i) = 12 ! jw
		end if
	end do
	
	
	! jw
	if(IS3D_surf_switch == 1) then
		IS3D_File_name_surf = trim(id_vtk3D_surf)//trim(File_num_buff)//'.vtk'
		open(pw_vtk3D_surf, file = trim(IS3D_File_name_surf), form = 'formatted', status = 'replace')	! jw
		write(pw_vtk3D_surf,'(A)') '# vtk DataFile Version 3.0'
		write(pw_vtk3D_surf,'(A)') 'Title = "3D contour plot"'
		write(pw_vtk3D_surf,'(A)') 'ASCII'
		write(pw_vtk3D_surf,*)
		write(pw_vtk3D_surf,'(A)') 'DATASET UNSTRUCTURED_GRID'
		write(pw_vtk3D_surf,'(A,I10,A)') 'POINTS ', maxnod, ' float'
		! jw
		do i=1,maxnod
			write(pw_vtk3D_surf,'(F15.5, F15.5, F15.5)') 	&
			&	(x_node(i)-xn_min) * IS3D_unit_conv, 			& 	! jw
			&	(y_node(i)-yn_min) * IS3D_unit_conv,			&	! jw
			&	eta_cell(i)													! jw
		end do
		write(pw_vtk3D_surf,*)
		! jw
		write(pw_vtk3D_surf,'(A,I10,I10)') 'CELLS ', maxele, cell_connectivity_size_2D
		do i=1,maxele
			! jw
			! jw
			! jw
			! jw
			! jw
			write(pw_vtk3D_surf,format_cell_connectivity_2D(i)) tri_or_quad(i), (nodenum_at_cell(j,i)-1, j=1,tri_or_quad(i)) ! jw
		end do
		write(pw_vtk3D_surf,*)
		write(pw_vtk3D_surf,'(A,I10)') 'CELL_TYPES ', maxele
		do i=1,maxele
			! jw
			! jw
			! jw
			! jw
			! jw
			write(pw_vtk3D_surf,'(I2)') cell_type_2D(i)
		end do
		close(pw_vtk3D_surf)
	end if
	
	
	! jw
	if(IS3D_full_switch == 1) then
		IS3D_File_name_full = trim(id_vtk3D_full)//trim(File_num_buff)//'.vtk'
		open(pw_vtk3D_full, file = trim(IS3D_File_name_full), form = 'formatted', status = 'replace')	! jw
		write(pw_vtk3D_full,'(A)') '# vtk DataFile Version 3.0'
		write(pw_vtk3D_full,'(A)') 'Title = "3D contour plot"'
		write(pw_vtk3D_full,'(A)') 'ASCII'
		write(pw_vtk3D_full,*)
		write(pw_vtk3D_full,'(A)') 'DATASET UNSTRUCTURED_GRID'
		
		! jw
		write(pw_vtk3D_full,'(A,I10,A)') 'POINTS ', maxnod*(maxlayer+1), ' float'
		do k=0,maxlayer
			do i=1,maxnod
				write(pw_vtk3D_full,'(F15.5, F15.5, F15.5)') 	&
				&	(x_node(i)-xn_min) * IS3D_unit_conv, 			& 	! jw
				&	(y_node(i)-yn_min) * IS3D_unit_conv,			&	! jw
				&	-h_node(i) + ((eta_node(i) + h_node(i))/maxlayer)*k
			end do
		end do
		write(pw_vtk3D_full,*)
		
		! jw
		! jw
		! jw
		write(pw_vtk3D_full,'(A,I10,I10)') 'CELLS ', maxele*maxlayer, cell_connectivity_size_3D*maxlayer
		do k=1,maxlayer
			do i=1,maxele
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				write(pw_vtk3D_full,format_cell_connectivity_3D(i)) &
				&	tri_or_quad(i)*2, &
				&	(nodenum_at_cell(j,i)-1 + maxnod*(k-1), j=1,tri_or_quad(i)), &	! jw
				&	(nodenum_at_cell(j,i)-1 + maxnod*k,     j=1,tri_or_quad(i))		! jw
			end do
		end do
		write(pw_vtk3D_full,*)
		
		write(pw_vtk3D_full,'(A,I10)') 'CELL_TYPES ', maxele*maxlayer
		do k=1,maxlayer
			do i=1,maxele
				! jw
				! jw
				! jw
				! jw
				! jw
				write(pw_vtk3D_full,'(I2)') cell_type_3D(i) ! jw
			end do
		end do
		write(pw_vtk3D_full,*)
		
		! jw
		! jw
		if(IS3D_variable(1) == 1 .or. IS3D_variable(2) == 1 .or. IS3D_variable(3) == 1) then
			write(pw_vtk3D_full,'(A,I10)') 'POINT_DATA ', maxnod*(maxlayer+1)
			write(pw_vtk3D_full,'(A)') 'VECTORS Velocity(m/s) float'
			do k=0,maxlayer
				do i=1,maxnod
					write(pw_vtk3D_full,'(F15.5 F15.5 F15.5)') &
					&	u_sigma(k,i), v_sigma(k,i), w_sigma(k,i)
				end do
			end do
		end if
		write(pw_vtk3D_full,*)
		
		! jw
		! jw
		if(IS3D_variable(4) == 1 .or. & 	! jw
		&	IS3D_variable(5) == 1 .or. & 	! jw
		&	IS3D_variable(6) == 1) then	! jw
			write(pw_vtk3D_full,'(A,I10)') 'CELL_DATA ', maxele*maxlayer
		end if
		if(IS3D_variable(4) == 1) then
			write(pw_vtk3D_full,'(A)') 'SCALARS Salt(psu) float'
			write(pw_vtk3D_full,'(A)') 'LOOKUP_TABLE default'
			do k=1,maxlayer
				do i=1,maxele
					write(pw_vtk3D_full,'(F15.5)') salt_sigma(k,i)
				end do
			end do
			write(pw_vtk3D_full,*)
		end if
		
		if(IS3D_variable(5) == 1) then
			write(pw_vtk3D_full,'(A)') 'SCALARS Temp(C) float'
			write(pw_vtk3D_full,'(A)') 'LOOKUP_TABLE default'
			do k=1,maxlayer
				do i=1,maxele
					write(pw_vtk3D_full,'(F15.5)') temp_sigma(k,i)
				end do
			end do
			write(pw_vtk3D_full,*)
		end if
		
		if(IS3D_variable(6) == 1) then
			write(pw_vtk3D_full,'(A)') 'SCALARS Rho(kg/m3) float'
			write(pw_vtk3D_full,'(A)') 'LOOKUP_TABLE default'
			do k=1,maxlayer
				do i=1,maxele
					write(pw_vtk3D_full,'(F15.5)') rho_sigma(k,i)
				end do
			end do
			write(pw_vtk3D_full,*)
		end if
		
		
		close(pw_vtk3D_full)		
	end if
end subroutine write_vtk3D