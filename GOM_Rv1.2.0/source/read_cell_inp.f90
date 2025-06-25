!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! This subroutine reads cell.inp
!! 
subroutine read_cell_inp
	use mod_global_variables
	use mod_file_definition
	use mod_function_library, only : calculate_area
	
	implicit none
	integer :: i, j
	integer :: cell_num
   integer :: num_element, node1, node2, node3, node4
   real(8) :: temp_area1, temp_area2, temp_area3, temp_area4      
	! end of local variables ==================================================!
	
   ! allocate & initialize some variables ====================================!
	allocate(tri_or_quad(maxele), area(maxele))		! (maxele)
   allocate(nodenum_at_cell(4,maxele))					 	! (maxele,4)
   allocate(nodenum_at_cell_tec(4,maxele))
   tri_or_quad = 0
   area = 0.0
   nodenum_at_cell = 0
   nodenum_at_cell_tec = 0
	! end of allocation =======================================================!
	

	write(pw_run_log,*) "	Read cell.inp"
	write(pw_run_log,*) "		Now, you are in 'read_input.f90 -> subroutine read_cell_inp"
	
	! cell.inp should exist, otherwise it will give you an error message
	open(pw_cell_inp, file = id_cell_inp, form='formatted', status = 'old')
	if(cell_mirr == 1) then
		open(pw_cell_mirr, file = id_cell_mirr, form = 'formatted', status = 'replace')
	end if

	! skip header lines
	call skip_header_lines(pw_cell_inp,id_cell_inp)
	
	! read main body ==========================================================!
	! read total_cell_number
	read(pw_cell_inp,*) cell_num
	if(cell_mirr == 1) then
		write(pw_cell_mirr,*) cell_num
	end if
	
	if(cell_num /= MAXELE) then
		write(pw_run_log,*) "Total cell number does not mathch: STOP"
		stop 'read_cell_inp.f90, Error #2'
	end if
	
	! read main body of cell.inp ==============================================!
   do i=1,maxele
      read(pw_cell_inp,*)  num_element, tri_or_quad(i), (nodenum_at_cell(j,i), j = 1, tri_or_quad(i))
      
		! erea calculation for triangle ========================================!
      node1 = nodenum_at_cell(1,i)
      node2 = nodenum_at_cell(2,i)
      node3 = nodenum_at_cell(3,i)
      
      if(tri_or_quad(i) == 3) then ! triangle
         area(i) = calculate_area(x_node(node1), x_node(node2), x_node(node3), &
         &								 y_node(node1), y_node(node2), y_node(node3))
      else	! quadrilateral
         node4=nodenum_at_cell(4,i)
         temp_area1 = calculate_area(x_node(node1), x_node(node2), x_node(node3),&
         &									 y_node(node1), y_node(node2), y_node(node3))
         temp_area2 = calculate_area(x_node(node1), x_node(node3), x_node(node4),& 
         &									 y_node(node1), y_node(node3), y_node(node4))
         temp_area3 = calculate_area(x_node(node1), x_node(node2), x_node(node4),& 
         &									 y_node(node1), y_node(node2), y_node(node4))
         temp_area4 = calculate_area(x_node(node2), x_node(node3), x_node(node4),& 
         &									 y_node(node2), y_node(node3), y_node(node4))

         if(temp_area1 <= 0.0 .or. temp_area2 <= 0.0 .or. temp_area3 <= 0.0 .or. temp_area4 <= 0.0) then
            write(pw_run_log,*)  'concave quadrangle', i, temp_area1, temp_area2, temp_area3, temp_area4
            stop
         end if
         area(i) = temp_area1 + temp_area2	! we need just two triangles 
      end if
      
      if(area(i) <= 0.0) then
         write(pw_run_log,*)'negative area at', i
         stop 'read_cell_inp.f90, Error #3'
      end if
   end do
   
   cell_connectivity_size_2D = maxele + sum(tri_or_quad) 	! total number of data in 'CELLS' for vtk file format, for 2D
   cell_connectivity_size_3D = maxele + sum(tri_or_quad)*2	! total number of data in 'CELLS' for vtk file format, for 3D
   
	close(pw_cell_inp)		! Close 'cell.inp'
	! end of main body of cell.inp ============================================!
	
   ! write mirror image file, cell_mirr.dat
 	if(cell_mirr == 1) then
 		do i=1,MAXELE
 			! write(pw_cell_mirr,*) element(i), (node(k,i),k=1,4), bottom_roughness(i), bed_elev(i), h_old(i)
 			write(pw_cell_mirr,*) i, tri_or_quad(i), (nodenum_at_cell(j,i),j=1, tri_or_quad(i))
 		end do
 		close(pw_cell_mirr)
 	end if

	! Create tecplot type node information
	! If triangle, the fourth node should be equal to the third node
	! This is for the tecplot output file
 	nodenum_at_cell_tec = nodenum_at_cell
 	do i=1,maxele
 		if(nodenum_at_cell(4,i) == 0) then
 			nodenum_at_cell_tec(4,i) = nodenum_at_cell(3,i)
 		end if
 	end do
end subroutine read_cell_inp
