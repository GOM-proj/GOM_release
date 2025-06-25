!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine set_geometry_2
	use mod_global_variables
	use mod_file_definition
	implicit none
	
 	integer :: i, j, k
	
	integer :: ob_node_flag(maxnod)
	integer :: nd1, nd2, el, fc
	! end of local variables ==================================================!

	! (8) determin the open boundary elements and faces =======================!
	! initialize open boundary flags
   ob_node_flag = 0	! ob_node_flag(maxnod,2)
	
	! set open boundary flag to 1 if the node is the open boundary node  
	do i=1,num_ob_cell
		nd1 = ob_nodes(i,1) ! first node at this ob_cell
		nd2 = ob_nodes(i,2) ! second node at this ob_cell
		ob_node_flag(nd1) = 1
		ob_node_flag(nd2) = 1
	end do
	
	write(pw_dia_geometry,*) '8. check_geo(8): ob_node_flag(i)'
	do i=1,maxnod
		write(pw_dia_geometry,'(A9,I10,A1, I10)') &
		&	' Node# = ', i, ',', ob_node_flag(i)
	end do
	write(pw_dia_geometry,*)

	! (9) Find open boundary element and face number (ID) =====================!
	do i=1,num_ob_cell
		! ob_nodes should be counterclockwise
		call find_element_and_face(ob_nodes(i,1),ob_nodes(i,2)) ! it will return element_id & face_id
		
		if(element_id == 0 .or. face_id == 0) then
			stop 'Fails to find the boundary element and face: stop.'
		end if
		
		ob_cell_id(i) = element_id
		ob_face_id(i) = face_id
	end do
	
	write(pw_dia_geometry,*) '9. check_geo(9): ob_cell_id(num_ob_cell), ob_face_id(num_ob_cell)'
	do i=1,num_ob_cell
		write(pw_dia_geometry,'(2I10)') ob_cell_id(i), ob_face_id(i)
	end do
	write(pw_dia_geometry,*)


	! (10) classify boundary elements & faces =================================!
	! initialize variables
	! boundary_type_of_face:
	! 	  -1	: land boundary face, i.e., river boundary face belongs to this
	! 		0	: interior water face
	! 		1~n: ith open boundary face, i.e., not including river boundary
	boundary_type_of_face = 0
   
   do i=1, num_ob_cell
	   el = ob_cell_id(i)
	   fc = ob_face_id(i)
		
		! this is the original method, which sets to 1	   
	   ! ob_element_flag(el) = 1 
	   ! boundary_type_of_face(fc) = 1
	   
	   ! this is the new method, which sets to 1 ~ n
	   ob_element_flag(el) = i
	   boundary_type_of_face(fc) = i
	end do
   
   do j = 1, maxface
   	! if the face does not have the second element (or belongs only to one element), and
   	! if the face is not asigned as one of the open boundary faces,
   	! then, set the face as a land boundary face (boundary between water and land; including river boundary face).
      if(adj_cellnum_at_face(2,j) == 0 .and. boundary_type_of_face(j) == 0) then
         boundary_type_of_face(j) = -1
      end if
   end do
	
	write(pw_dia_geometry,*) '10. check_geo(10)'
	write(pw_dia_geometry,*) '(10.1) ob_element_flag(i): ith = boundary element, 0 = non openboundary element'
	do i=1,maxele
		write(pw_dia_geometry,'(A9,I10,A1, I10)') &
		&	' Cell# = ', i, ',', ob_element_flag(i)
	end do
	write(pw_dia_geometry,*)
	write(pw_dia_geometry,*) '(10.2) boundary_type_of_face(j): -1: land face, 0: inner water face, 1~n: open boundary face'
	do j=1,maxface
		write(pw_dia_geometry,'(A9,I10,A1, I10)') &
		&	' Face# = ', j, ',', boundary_type_of_face(j)
	end do
	write(pw_dia_geometry,*)

	! (11) find the [bottom layer #] at each element/face/node ================!
	write(pw_dia_geometry,*) '11. check_geo(11): bottom layer at element/face/node'
	write(pw_dia_geometry,*) '(11.1) bottom_layer_at_element(i)'
	do i=1,maxele		
	   bottom_layer_at_element(i) = 0
	   do k = 0, maxlayer - 1
	      if((MSL - h_cell(i)) >= z_level(k) .and. &
	      &  (MSL - h_cell(i)) <  z_level(k+1)) then
	         bottom_layer_at_element(i) = k+1  ! bottom actually at level k, not k+1
	         exit
	      end if
	   end do
	   
	   if(bottom_layer_at_element(i) == 0) then
	   	! write an error and stop
	      write(pw_run_log,*)'either MSL < depth at element',i
	      write(pw_run_log,*)'or land height above highest level'
	      write(pw_run_log,*)'check vertical grid information, C2 & C2_1, in main.inp'
	      write(pw_run_log,*)'MSL =', MSL
	      write(pw_run_log,*)'element number = ', i
	      write(pw_run_log,*)'h_cell(i) = ', h_cell(i)

	      write(*,*)'either MSL < depth at element',i
	      write(*,*)'or land height above highest level'
	      write(*,*)'check vertical grid information, C2 & C2_1, in main.inp'
	      write(*,*)'MSL =', MSL
	      write(*,*)'element number = ', i
	      write(*,*)'h_cell(i) = ', h_cell(i)
	      stop 'Vertical grid information error'
	   else
	   	write(pw_dia_geometry,'(A9,I10,A1, I10)') &
	   	&	' Cell# = ', i, ',', bottom_layer_at_element(i)
	   end if
	end do
	write(pw_dia_geometry,*)
	write(pw_dia_geometry,*) '(11.2) bottom_layer_at_face(j)'
	do j=1,maxface
      bottom_layer_at_face(j) = 0
      do k = 0, maxlayer-1
         if(MSL-h_face(j) >= z_level(k) .and. &
         &  MSL-h_face(j) < z_level(k+1)) then
            bottom_layer_at_face(j) = k + 1
            exit
         end if
      end do 

     if(bottom_layer_at_face(j) == 0) then
         write(pw_run_log,*)'either MSL < depth at face', j
         write(pw_run_log,*)'or land height above highest level'
         write(pw_run_log,*)'check vertical grid information, C2 & C2_1, in main.inp'
         write(pw_run_log,*)'MSL =', MSL
         write(pw_run_log,*)'face number = ', j
         write(pw_run_log,*)'h_face(j) = ', h_face(j)

         write(*,*)'either MSL < depth at face',j
         write(*,*)'or land height above highest level'
         write(*,*)'check vertical grid information, C2 & C2_1, in main.inp'
         write(*,*)'MSL =', MSL
         write(*,*)'face number = ', j
         write(*,*)'h_face(j) = ', h_face(j)
         stop 'Vertical grid information error'
      else
      	write(pw_dia_geometry,'(A9,I10,A1, I10)') &
      	&	' Face# = ', j, ',', bottom_layer_at_face(j)
      end if
   end do
	write(pw_dia_geometry,*)
	write(pw_dia_geometry,*) '(11.3) bottom_layer_at_node(i)'
	do i=1,maxnod
      bottom_layer_at_node(i) = 0
      do k = 0, maxlayer-1
         if(MSL-h_node(i) >= z_level( k ) .and. &
         &  MSL-h_node(i) <  z_level(k+1)) then
            bottom_layer_at_node(i) = k+1
            exit
         end if
      end do 

   	if(bottom_layer_at_node(i) == 0) then
         write(pw_run_log,*)'either MSL < depth at side',i
         write(pw_run_log,*)'or land height above highest level'
         write(pw_run_log,*)'check vertical grid information, C2 & C2_1, in main.inp'
         write(pw_run_log,*)'MSL =', MSL
         write(pw_run_log,*)'node number = ', i
         write(pw_run_log,*)'h_node(i) = ', h_node(i)
         write(*,*) 'Vertical grid information error'
         stop
      else
      	write(pw_dia_geometry,'(A9,I10,A1, I10)') &
      	&	' Node# = ', i, ',', bottom_layer_at_node(i)
      end if   ! if(bottom_layer_at_node(i) == 0)
	end do
	write(pw_dia_geometry,*)
	write(pw_dia_geometry,*) 'End of file ====================================='
	close(pw_dia_geometry)
 end subroutine set_geometry_2