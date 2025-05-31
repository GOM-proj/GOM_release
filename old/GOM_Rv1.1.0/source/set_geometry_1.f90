!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! set mesh related geometries
!! General mesh geometries...
!!
subroutine set_geometry_1
	use mod_global_variables
	use mod_file_definition
	implicit none
	
 	integer :: i, j, k, l, l2
 	integer :: ie, nd, nd1, nd2, num_dummy
 	integer :: node_temp
	
	real(dp), allocatable :: xn1_at_face(:), xn2_at_face(:), yn1_at_face(:), yn2_at_face(:)
 	integer :: i_temp_element, index1
 	real(dp):: angle_theta0, angle_theta1, angle_theta2
 	integer :: jsj
 	real(dp):: x1, y1, x2, y2
 	real(dp):: x_sum, y_sum
	! jw
	
	open(pw_dia_geometry,file=id_dia_geometry,form = 'formatted', status = 'replace')
	write(pw_dia_geometry,*) '!============================================================================!'
	write(pw_dia_geometry,*) 'This file, check_geometry.dia, will show general mesh geometric information:'
	write(pw_dia_geometry,*) '1: check_geo(1)'
	write(pw_dia_geometry,*) '   (1.1) adj_cells_at_node'
	write(pw_dia_geometry,*) '   (1.2) adj_cellnum_at_node'
	write(pw_dia_geometry,*) '   (1.3) node_count_each_element'
	write(pw_dia_geometry,*) '2: check_geo(2)'
	write(pw_dia_geometry,*) '   (2.1) adj_nodes_at_node'
	write(pw_dia_geometry,*) '   (2.2) adj_nodenum_at_node'
	write(pw_dia_geometry,*) '3. check_geo(3)'
	write(pw_dia_geometry,*) '   (3.1) adj_cellnum_at_cell'
	write(pw_dia_geometry,*) '4. check_geo(4)'
	write(pw_dia_geometry,*) '   (4.1) maxface'
	write(pw_dia_geometry,*) '   (4.2) face constructing cells, nodes, face length, face center coordinates, depth at face'
	write(pw_dia_geometry,*) '         adj_cellnum_at_face(1,j), adj_cellnum_at_face(2,j),&
									           & nodenum_at_face(1,j), nodenum_at_face(2,j), face_length(j), x_face(j), y_face(j), h_face(j)'
	write(pw_dia_geometry,*) '   (4.3) facenum_at_cell'
	write(pw_dia_geometry,*) '5. check_geo(5)'
	write(pw_dia_geometry,*) '   (5.1) sign_in_outflow'
	write(pw_dia_geometry,*) '6. check_geo(6): Cell center inforamtion'
	write(pw_dia_geometry,*) '7. check_geo(7): unit normal vector at each face'
	write(pw_dia_geometry,*) '8. check_geo(8): ob_node_flag'
	write(pw_dia_geometry,*) '9. check_geo(9): ob_cell_id(num_ob_cell), ob_face_id(num_ob_cell)'
	write(pw_dia_geometry,*) '10. check_geo(11)'
	write(pw_dia_geometry,*) '   (10.1) ob_element_flag'
	write(pw_dia_geometry,*) '   (10.2) boundary_type_of_face: -1: land face, 0: inner water face, 1: open boundary face'
	write(pw_dia_geometry,*) '!============================================================================!'

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
   ! jw
   
   do k = 3, 4		! jw
      do i = 1, k
         do j = 1, k-1
            start_end_node(k,i,j) = i + j
            
            if(start_end_node(k,i,j) > k ) then
               start_end_node(k,i,j) = start_end_node(k,i,j)-k
            end if

            if(start_end_node(k,i,j) < 1 .or.    				&
            &  start_end_node(k,i,j) > k) then
               write(*,*)'set_geometry_1.f90: start_end_node wrong',   				&
               &        i, j, k, start_end_node(k,i,j)
               stop
            end if
            
            ! jw
         end do
      end do
   end do
	

	!==========================================================================!
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	allocate(adj_cells_at_node(maxnod))
   allocate(adj_cellnum_at_node(max_no_neighbor_node,maxnod))
	allocate(node_count_each_element(max_no_neighbor_node,maxnod))
	adj_cells_at_node 		= 0
	adj_cellnum_at_node 		= 0
	node_count_each_element = 0

   do i = 1, maxele
      do l = 1, tri_or_quad(i)
         node_temp = nodenum_at_cell(l,i)
         adj_cells_at_node(node_temp) = adj_cells_at_node(node_temp) + 1
         if(adj_cells_at_node(node_temp) > max_no_neighbor_node) then
            write(pw_run_log,*) 'Too many neighbors at node: ', node_temp
            stop 'set_geometry, Error #1_1'
         end if
         adj_cellnum_at_node(adj_cells_at_node(node_temp),node_temp) = i
         node_count_each_element(adj_cells_at_node(node_temp),node_temp)   = l
      end do
   end do

   ! jw
   do i = 1, maxnod
      if(adj_cells_at_node(i) == 0) then
         write(pw_run_log,*)'hanging node',i
         stop 'set_geometry, Error #1_2'
      end if
   end do

	write(pw_dia_geometry,*) '1: check_geo(1)'
 	write(pw_dia_geometry,*) '(1.1) adj_cells_at_node(i)'
 	do i=1,maxnod
 		write(pw_dia_geometry,'(A9,I10,A1, I10)') &
 		&	' Node# = ', i, ',', adj_cells_at_node(i)
 	end do
 	write(pw_dia_geometry,*)
 	write(pw_dia_geometry,*) '(1.2) adj_cellnum_at_node(j,i), j=1,max_no_neighbor_node'
 	do i=1,maxnod
 		write(pw_dia_geometry,'(A9,I10,A1, *(I10))') &
 		&	' Node# = ', i, ',', (adj_cellnum_at_node(j,i), j=1,max_no_neighbor_node)
 	end do
 	write(pw_dia_geometry,*)
 	write(pw_dia_geometry,*) '(1.3) node_count_each_element(j,i), j=1,max_no_neighbor_node'
 	do i=1,maxnod
 		write(pw_dia_geometry,'(A9,I10,A1, *(I10))') &
 		&	' Node# = ', i, ',', (node_count_each_element(j,i), j=1,max_no_neighbor_node)
 	end do
 	write(pw_dia_geometry,*)
   
   ! jw
	! jw
	! jw
	allocate(adj_nodes_at_node(maxnod))
	allocate(adj_nodenum_at_node(0:max_no_neighbor_node,maxnod))
	adj_nodes_at_node 		= 0
	adj_nodenum_at_node 		= 0

   do i = 1, maxnod
      do j = 1, adj_cells_at_node(i)
         ie = adj_cellnum_at_node(j,i)	! jw
         
         loop_element: do k = 1, tri_or_quad(ie)	! jw
         	nd = nodenum_at_cell(k,ie)	! jw
            do l = 1, adj_nodes_at_node(i)
               if(nd == i .or. nd == adj_nodenum_at_node(l,i)) then	! jw
                  cycle loop_element
               end if
            end do
             
            adj_nodes_at_node(i) = adj_nodes_at_node(i) + 1
             
            if(adj_nodes_at_node(i) > max_no_neighbor_node) then
               write(pw_run_log,*)'Too many surrounding nodes at node:', i
               stop 'set_geometry, Error #2'
            end if
            adj_nodenum_at_node(adj_nodes_at_node(i),i) = nd
   		end do loop_element
      end do
   end do

	write(pw_dia_geometry,*) '2: check_geo(2)'
	write(pw_dia_geometry,*) '(2.1) adj_nodes_at_node(i)'
	do i=1,maxnod
		write(pw_dia_geometry,'(A9,I10,A1, I10)') &
		&	' Node# = ', i, ',', adj_nodes_at_node(i)
	end do
	write(pw_dia_geometry,*)
	write(pw_dia_geometry,*) '(2.2) adj_nodenum_at_node(j,i), j=1,max_no_neighbor_node'
	do i=1,maxnod
		write(pw_dia_geometry,'(A9,I10,A1, *(I10))') &
		&	' Node# = ', i, ', ', (adj_nodenum_at_node(j,i), j=1,max_no_neighbor_node)
	end do
	write(pw_dia_geometry,*)

	! jw
	! jw
	! jw
	! jw
	allocate(adj_cellnum_at_cell(4,maxele))								! jw
	adj_cellnum_at_cell 		= 0

   do i = 1, maxele
      do l = 1, tri_or_quad(i)
         adj_cellnum_at_cell(l,i) = 0
         nd1 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,1),i)
         nd2 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,2),i)
         do k = 1, adj_cells_at_node(nd1)
            num_dummy = adj_cellnum_at_node(k,nd1)
            if(num_dummy /= i .and.	&
            &	(nodenum_at_cell(1,num_dummy) == nd2 .or.	&
            &	 nodenum_at_cell(2,num_dummy) == nd2 .or.	&
            &	 nodenum_at_cell(3,num_dummy) == nd2 .or.	&
            &	(tri_or_quad(num_dummy)  == 4   .and.  &
            &	 nodenum_at_cell(4,num_dummy) == nd2))) then
               adj_cellnum_at_cell(l,i) = num_dummy
            end if            
         end do ! jw
      end do ! jw
   end do ! jw
	
	write(pw_dia_geometry,*) '3. check_geo(3)'
	write(pw_dia_geometry,*) '(3.1) adj_cellnum_at_cell(l,i), l=1,4'
	do i=1,maxele
		write(pw_dia_geometry,'(A9,I10,A1, 4I10)') &
		&	' Cell# = ', i, ', ', (adj_cellnum_at_cell(l,i), l=1,4)
	end do
	write(pw_dia_geometry,*)

	! jw
	! jw
	! jw
	! jw
	! jw
   maxface = 0 ! jw
   do i = 1, maxele
      do l = 1, tri_or_quad(i)
         nd1 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,1),i)
         nd2 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,2),i)

         if(adj_cellnum_at_cell(l,i) == 0 .or. i < adj_cellnum_at_cell(l,i)) then ! jw
            maxface = maxface + 1 
         end if ! jw
      end do ! jw
   end do ! jw
   
	! jw
   if(maxface < maxele .or. maxface < maxnod) then
      write(pw_run_log,*) 'Weird grid with maxface < maxele or maxface < maxnod', maxnod, maxele, maxface
      stop 'set_geometry, Error #4_2'
   end if
   
   write(*,'(A30,I10)') 'maxface = ', maxface
   write(*,*)
	write(*,*) 'Now I am allocating variables & reading input files & preparing GOM...'
	write(*,*) 'Please wait a while...'
 	write(*,*) '.'
 	write(*,*) '.'
 	write(*,*) '.'
 	
	! jw
   allocate(facenum_at_cell(4,maxele))
	facenum_at_cell = 0

	allocate(face_length(maxface), 				&
	&			x_face(maxface), 						&
	&			y_face(maxface), 						&
	&			h_face(maxface),						&
	&			adj_cellnum_at_face(2,maxface),	&
	&			nodenum_at_face(2,maxface))
	face_length 			= 0.0_dp
	x_face					= 0.0_dp
	y_face					= 0.0_dp
	h_face					= 0.0_dp
	adj_cellnum_at_face	= 0
	nodenum_at_face		= 0
	
	! jw
	allocate(xn1_at_face(maxface), 		&
	&			xn2_at_face(maxface), 		&
	&			yn1_at_face(maxface), 		&
	&			yn2_at_face(maxface))
	xn1_at_face = 0.0_dp
	xn2_at_face = 0.0_dp
	yn1_at_face = 0.0_dp
	yn2_at_face = 0.0_dp
	   
   ! jw
   ! jw
	maxface = 0   
   do i = 1, maxele
      do l = 1, tri_or_quad(i)
         nd1 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,1),i)
         nd2 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,2),i)

         if(adj_cellnum_at_cell(l,i) == 0 .or. i < adj_cellnum_at_cell(l,i)) then ! jw
            maxface = maxface + 1 
            facenum_at_cell(l,i)  = maxface
            adj_cellnum_at_face(1,maxface) = i
            nodenum_at_face(1,maxface) = nd1   ! jw
            nodenum_at_face(2,maxface) = nd2   ! jw

            xn1_at_face(maxface) = x_node(nd1)	! jw
            yn1_at_face(maxface) = y_node(nd1)	! jw
            xn2_at_face(maxface) = x_node(nd2)	! jw
            yn2_at_face(maxface) = y_node(nd2)	! jw
            
            x_face(maxface) = (xn1_at_face(maxface) + xn2_at_face(maxface)) * 0.5
            y_face(maxface) = (yn1_at_face(maxface) + yn2_at_face(maxface)) * 0.5
            
            h_face(maxface) = (h_node(nd1) + h_node(nd2)) * 0.5	! jw
            face_length(maxface) = dsqrt((xn2_at_face(maxface) - xn1_at_face(maxface))**2 + &
            &										(yn2_at_face(maxface) - yn1_at_face(maxface))**2) 

            if(face_length(maxface) == 0) then
               write(pw_run_log,*) 'zero face_length', maxface
               stop 'set_geometry, Error #4_3'
            end if
            
            adj_cellnum_at_face(2,maxface) = adj_cellnum_at_cell(l,i)
            if(adj_cellnum_at_cell(l,i) /= 0) then
               i_temp_element = adj_cellnum_at_cell(l,i)
               index1 = 0

               do l2 = 1, tri_or_quad(i_temp_element)
                  if(adj_cellnum_at_cell(l2,i_temp_element) == i) then
                     index1 = l2
                     exit
                  end if
               end do

               if(index1 == 0) then
                  write(pw_run_log,*) 'wrong ball info',i,l		! jw
                  stop 'set_geometry, Error #4_4'
               end if
               facenum_at_cell(index1,i_temp_element) = maxface
            end if ! jw
         end if ! jw
      end do ! jw
   end do ! jw

	
	! jw
	! jw
	write(pw_dia_geometry,*) '4. check_geo(4)'
	write(pw_dia_geometry,*) '(4.1) total face number:'
	write(pw_dia_geometry,*) 'maxface = ', maxface	
	write(pw_dia_geometry,*)	
	write(pw_dia_geometry,*) '(4.2) face constructing cells, nodes, face length, face center coordinates, depth at face:'
	write(pw_dia_geometry,*) 'adj_cellnum_at_face(1,j), adj_cellnum_at_face(2,j), &
	&	nodenum_at_face(1,j), nodenum_at_face(2,j), face_length(j), x_face(j), y_face(j), h_face(j)'
   do j = 1, maxface
      ! jw
      write(pw_dia_geometry,'(A9,I10,A1, 4I10,4E20.10)') &
      &	' Face# = ', j, ',', adj_cellnum_at_face(1,j), adj_cellnum_at_face(2,j), nodenum_at_face(1,j), nodenum_at_face(2,j), &
      &	face_length(j), x_face(j), y_face(j), h_face(j)
   end do
   write(pw_dia_geometry,*)
   write(pw_dia_geometry,*) '(4.3) facenum_at_cell(l,i)'
 	do i=1,maxele
 		write(pw_dia_geometry,'(A9,I10,A1, 4I10)') &
 		&	' Cell# = ', i, ',', (facenum_at_cell(l,i), l=1,4)
 	end do
	write(pw_dia_geometry,*)

	! jw
	! jw
	! jw
   allocate(sign_in_outflow(4,maxele))
   sign_in_outflow = 0
	
   do i = 1, maxele
      do l = 1, tri_or_quad(i)
         jsj = facenum_at_cell(l,i)
         sign_in_outflow(l,i) = (adj_cellnum_at_face(2,jsj) - 2*i + adj_cellnum_at_face(1,jsj)) / &
         &							  (adj_cellnum_at_face(2,jsj) - adj_cellnum_at_face(1,jsj))
      end do
   end do
	
	write(pw_dia_geometry,*) '5. check_geo(5)'
	write(pw_dia_geometry,*) '(5.1) sign_in_outflow(l,i)'
	do i=1,maxele
		write(pw_dia_geometry,'(A9,I10,A1, 4I10)') &
		&	' Cell# = ', i, ', ', (sign_in_outflow(l,i), l=1,4)
	end do
	write(pw_dia_geometry,*)

	! jw
   allocate(x_cell(maxele), &
   &			y_cell(maxele), &
   &			h_cell(maxele))		
   x_cell 	= 0.0_dp
   y_cell 	= 0.0_dp
   h_cell 	= 0.0_dp	! jw
	
   do i = 1, maxele
   	x_sum = 0.0_dp
   	y_sum = 0.0_dp
      
      ! jw
      ! jw
      do l = 1, tri_or_quad(i)
      	x_sum = x_sum + x_node(nodenum_at_cell(l,i))
      	y_sum = y_sum + y_node(nodenum_at_cell(l,i))

         ! jw
         if(h_face(facenum_at_cell(l,i)) > h_cell(i)) then
            h_cell(i) = h_face(facenum_at_cell(l,i))
         end if
      end do
      x_cell(i) = x_sum/tri_or_quad(i)
      y_cell(i) = y_sum/tri_or_quad(i)
   end do

	write(pw_dia_geometry,*) '6. check_geo(6): Cell center inforamtion'
	write(pw_dia_geometry,*) 'Note: if Voronoi == 1, voronoi center information will be used instead of the following information'
	write(pw_dia_geometry,*) 'x_cell(i), y_cell(i), h_cell(i)'
	do i=1,maxele
		! jw
		write(pw_dia_geometry,'(A9,I10,A1, 3E20.10)') &
		&	' Cell# = ', i, ', ', x_cell(i), y_cell(i), h_cell(i)
	end do
	write(pw_dia_geometry,*)

	! jw
	! jw
	! jw
	allocate(delta_j(maxface), 	&
	&			cos_theta(maxface), 	&
	&			sin_theta(maxface), 	&
	&			cos_theta2(maxface),	&
	&			sin_theta2(maxface))
	delta_j 		= 0.0_dp
	cos_theta 	= 0.0_dp
	sin_theta 	= 0.0_dp
	cos_theta2  = 0.0_dp
	sin_theta2  = 0.0_dp
	
	! jw
   if(voronoi == 1)then
   	! jw
		! jw
		! jw
		! jw
		! jw
      call calculate_voronoi_center(xn1_at_face, xn2_at_face, yn1_at_face, yn2_at_face)
   else
		! jw
      do j = 1, maxface
         if(adj_cellnum_at_face(2,j) /= 0) then
				! jw
				! jw
				! jw
				! jw
            delta_j(j) = dsqrt (   &
                    &           (x_cell(adj_cellnum_at_face(2,j))       &
                    &         -  x_cell(adj_cellnum_at_face(1,j)))**2   &
                    &         + (y_cell(adj_cellnum_at_face(2,j))       &
                    &         -  y_cell(adj_cellnum_at_face(1,j)))**2)
            if(delta_j(j) == 0) then
               write(pw_run_log,*) 'zero distance between centers at: ', adj_cellnum_at_face(1,j), adj_cellnum_at_face(2,j)
               stop 'set_geometry.f90, Error #11'
            end if
         end if
         
         if(adj_cellnum_at_face(2,j) == 0) then
         	! jw
         	! jw
         	! jw
         	! jw
         	delta_j(j) = dsqrt( &
         	&	(x_cell(adj_cellnum_at_face(1,j)) - x_face(j))**2 + &
         	&	(y_cell(adj_cellnum_at_face(1,j)) - y_face(j))**2) &
         	&	* 2.0
         	
         	! jw
         	! jw
         	! jw
         	! jw
         end if
      end do         
   end if
   
   ! jw
   ! jw
   ! jw
   ! jw
   ! jw
   do j=1,maxface
		! jw
		! jw
      angle_theta0 = datan2(xn1_at_face(j) - xn2_at_face(j), yn2_at_face(j) - yn1_at_face(j)) ! jw
      cos_theta(j) = dcos(angle_theta0)	! jw
      sin_theta(j) = dsin(angle_theta0)	! jw
      
      ! jw
      ! jw
      if(adj_cellnum_at_face(2,j) == 0) then
      	angle_theta1 = angle_theta0
      else
      	! jw
      	! jw
         x1 = x_cell(adj_cellnum_at_face(1,j))
         y1 = y_cell(adj_cellnum_at_face(1,j))
         
         x2 = x_cell(adj_cellnum_at_face(2,j))
         y2 = y_cell(adj_cellnum_at_face(2,j))
      	
         ! jw
         angle_theta1 = datan2(y2-y1,x2-x1) ! jw
      end if
      
      ! jw
      ! jw
      angle_theta2 = (angle_theta1 - angle_theta0)
      ! jw
      
      ! jw
      ! jw
      ! jw
      
      ! jw
      cos_theta2(j) = dcos(angle_theta2)
      sin_theta2(j) = dsin(angle_theta2)
	end do

   
	write(pw_dia_geometry,*) '7. check_geo(7): unit normal vector at each face'
	write(pw_dia_geometry,*) 'cos_theta(j), cos_theta2(j), sin_theta(j), sin_theta(2), delta_j(j), '
	do j=1,maxface
		! jw
		write(pw_dia_geometry,'(A9,I10,A1, 5E20.10)') &
		&	' Face# = ', j, ',', cos_theta(j), cos_theta2(j), sin_theta(j), sin_theta2(j), delta_j(j)
	end do
	write(pw_dia_geometry,*)
end subroutine set_geometry_1
