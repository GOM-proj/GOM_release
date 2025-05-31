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
	! jw
	
   ! jw
	allocate(tri_or_quad(maxele), area(maxele))		! jw
   allocate(nodenum_at_cell(4,maxele))					 	! jw
   allocate(nodenum_at_cell_tec(4,maxele))
   tri_or_quad = 0
   area = 0.0
   nodenum_at_cell = 0
   nodenum_at_cell_tec = 0
	! jw
	

	write(pw_run_log,*) "	Read cell.inp"
	write(pw_run_log,*) "		Now, you are in 'read_input.f90 -> subroutine read_cell_inp"
	
	! jw
	open(pw_cell_inp, file = id_cell_inp, form='formatted', status = 'old')
	if(cell_mirr == 1) then
		open(pw_cell_mirr, file = id_cell_mirr, form = 'formatted', status = 'replace')
	end if

	! jw
	call skip_header_lines(pw_cell_inp,id_cell_inp)
	
	! jw
	! jw
	read(pw_cell_inp,*) cell_num
	if(cell_mirr == 1) then
		write(pw_cell_mirr,*) cell_num
	end if
	
	if(cell_num /= MAXELE) then
		write(pw_run_log,*) "Total cell number does not mathch: STOP"
		stop 'read_cell_inp.f90, Error #2'
	end if
	
	! jw
   do i=1,maxele
      read(pw_cell_inp,*)  num_element, tri_or_quad(i), (nodenum_at_cell(j,i), j = 1, tri_or_quad(i))
      
		! jw
      node1 = nodenum_at_cell(1,i)
      node2 = nodenum_at_cell(2,i)
      node3 = nodenum_at_cell(3,i)
      
      if(tri_or_quad(i) == 3) then ! jw
         area(i) = calculate_area(x_node(node1), x_node(node2), x_node(node3), &
         &								 y_node(node1), y_node(node2), y_node(node3))
      else	! jw
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
         area(i) = temp_area1 + temp_area2	! jw
      end if
      
      if(area(i) <= 0.0) then
         write(pw_run_log,*)'negative area at', i
         stop 'read_cell_inp.f90, Error #3'
      end if
   end do
   
   cell_connectivity_size_2D = maxele + sum(tri_or_quad) 	! jw
   cell_connectivity_size_3D = maxele + sum(tri_or_quad)*2	! jw
   
	close(pw_cell_inp)		! jw
	! jw
	
   ! jw
 	if(cell_mirr == 1) then
 		do i=1,MAXELE
 			! jw
 			write(pw_cell_mirr,*) i, tri_or_quad(i), (nodenum_at_cell(j,i),j=1, tri_or_quad(i))
 		end do
 		close(pw_cell_mirr)
 	end if

	! jw
	! jw
	! jw
 	nodenum_at_cell_tec = nodenum_at_cell
 	do i=1,maxele
 		if(nodenum_at_cell(4,i) == 0) then
 			nodenum_at_cell_tec(4,i) = nodenum_at_cell(3,i)
 		end if
 	end do
end subroutine read_cell_inp
