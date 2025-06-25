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
	! end of local variables ==================================================!
	
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

   ! define start_end_node ===================================================!
   ! Following do loops as same as below:
   ! start_end_node will have:
   ! 	for triangular cell:
   ! 		start_end_node(3,1,1) = 2
   ! 		start_end_node(3,1,2) = 3
   ! 		start_end_node(3,2,1) = 3
   ! 		start_end_node(3,2,2) = 1
   ! 		start_end_node(3,3,1) = 1
   ! 		start_end_node(3,3,2) = 2
	
	! 	for quadrilateral cell:
   ! 		start_end_node(4,1,1) = 2
   ! 		start_end_node(4,1,2) = 3
   ! 		start_end_node(4,1,3) = 4
   ! 		start_end_node(4,2,1) = 3
   ! 		start_end_node(4,2,2) = 4
   ! 		start_end_node(4,2,3) = 1
   ! 		start_end_node(4,3,1) = 4
   ! 		start_end_node(4,3,2) = 1
   ! 		start_end_node(4,3,3) = 2
   ! 		start_end_node(4,4,1) = 1
   ! 		start_end_node(4,4,2) = 2
   ! 		start_end_node(4,4,3) = 3
   
   do k = 3, 4		! triangle -> quadrilateral
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
            
            ! write(*,*) k,i,j, start_end_node(k,i,j)
         end do
      end do
   end do
	

	!==========================================================================!
	! (1) find adjacent elements at each node =================================!
	! adj_cells_at_node(maxnod):
	! 		total number of adjacent cells at the node
	! adj_cellnum_at_node(1,2,3,...,adj_cells_at_node, maxnod):
	! 		adjacent cell number (ID) at each node
	! node_count_each_element(1,2,3,...,adj_cells_at_node, maxnod)
	! 		the location of this node appears in elements which include this node
	! 		i.e., this node is in the ith position in the adjacent element (i.e., this node appears in the ith position in cell.inp)
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

   ! check if hanging nodes (stand-alone nodes) exist
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
   
   ! (2) find adjacent nodes at each node ====================================!
	! adj_nodes_at_node(maxnod)											: total adjacent node numbers at each node
	! adj_nodenum_at_node(1,2,3,...,adj_cells_at_node, maxnod)	: adjacent node number (ID) at each node
	allocate(adj_nodes_at_node(maxnod))
	allocate(adj_nodenum_at_node(0:max_no_neighbor_node,maxnod))
	adj_nodes_at_node 		= 0
	adj_nodenum_at_node 		= 0

   do i = 1, maxnod
      do j = 1, adj_cells_at_node(i)
         ie = adj_cellnum_at_node(j,i)	! ie = cell ID which includes this node
         
         loop_element: do k = 1, tri_or_quad(ie)	! tri_orquad(ie) = 3 or 4
         	nd = nodenum_at_cell(k,ie)	! nd = nodes in the ie
            do l = 1, adj_nodes_at_node(i)
               if(nd == i .or. nd == adj_nodenum_at_node(l,i)) then	! if nd is the current node skip the rest part
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

	! (3) find adjacent cells at each element =================================!
	! adj_cellnum_at_cell(4,maxele)										: adjacent element number (ID) at each element, 
	! For a triangular cell, maximum 3 adjacent elements should exist. 
	! But for a quadrilateral cell, maximum 4 adjacent elements should exist.
	allocate(adj_cellnum_at_cell(4,maxele))								! (4,maxele), adjacent element number for face number 
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
         end do ! k
      end do ! l
   end do ! i
	
	write(pw_dia_geometry,*) '3. check_geo(3)'
	write(pw_dia_geometry,*) '(3.1) adj_cellnum_at_cell(l,i), l=1,4'
	do i=1,maxele
		write(pw_dia_geometry,'(A9,I10,A1, 4I10)') &
		&	' Cell# = ', i, ', ', (adj_cellnum_at_cell(l,i), l=1,4)
	end do
	write(pw_dia_geometry,*)

	! (4) determin face number of each element & 
	! find maxface (total face number) ========================================!
	! facenum_at_cell(4,maxele)											: face number of each element
	! jw, This method is not very good, but anyway use this method not to use max_face_num
	! First, find the total face number, maxface
   maxface = 0 ! total number of faces
   do i = 1, maxele
      do l = 1, tri_or_quad(i)
         nd1 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,1),i)
         nd2 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,2),i)

         if(adj_cellnum_at_cell(l,i) == 0 .or. i < adj_cellnum_at_cell(l,i)) then ! new sides
            maxface = maxface + 1 
         end if ! adj_cellnum_at_cell(l,i)==0.or.i<adj_cellnum_at_cell(l,i)
      end do ! l=1,tri_or_quad(i)
   end do ! i=1,maxele
   
	! additional error check for maxface calculation
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
 	
	! allocate & initialize global variables ==================================!
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
	
	! allocate & initialize local variables 
	allocate(xn1_at_face(maxface), 		&
	&			xn2_at_face(maxface), 		&
	&			yn1_at_face(maxface), 		&
	&			yn2_at_face(maxface))
	xn1_at_face = 0.0_dp
	xn2_at_face = 0.0_dp
	yn1_at_face = 0.0_dp
	yn2_at_face = 0.0_dp
	   
   ! let's do it again to calculate other variables...
   ! jw, this method is not very good but anyway use it this time...
	maxface = 0   
   do i = 1, maxele
      do l = 1, tri_or_quad(i)
         nd1 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,1),i)
         nd2 = nodenum_at_cell(start_end_node(tri_or_quad(i),l,2),i)

         if(adj_cellnum_at_cell(l,i) == 0 .or. i < adj_cellnum_at_cell(l,i)) then ! new sides
            maxface = maxface + 1 
            facenum_at_cell(l,i)  = maxface
            adj_cellnum_at_face(1,maxface) = i
            nodenum_at_face(1,maxface) = nd1   ! start node number of face of element
            nodenum_at_face(2,maxface) = nd2   ! end   node number of face of element

            xn1_at_face(maxface) = x_node(nd1)	! start node coordinate
            yn1_at_face(maxface) = y_node(nd1)	! start node coordinate
            xn2_at_face(maxface) = x_node(nd2)	! end node coordinate
            yn2_at_face(maxface) = y_node(nd2)	! end node coordinate
            
            x_face(maxface) = (xn1_at_face(maxface) + xn2_at_face(maxface)) * 0.5
            y_face(maxface) = (yn1_at_face(maxface) + yn2_at_face(maxface)) * 0.5
            
            h_face(maxface) = (h_node(nd1) + h_node(nd2)) * 0.5	! we are using water depth at the face center
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
                  write(pw_run_log,*) 'wrong ball info',i,l		! jw, I don't know what this is.
                  stop 'set_geometry, Error #4_4'
               end if
               facenum_at_cell(index1,i_temp_element) = maxface
            end if ! adj_cellnum_at_cell(l,i)/=0
         end if ! adj_cellnum_at_cell(l,i)==0.or.i<adj_cellnum_at_cell(l,i)
      end do ! l=1,tri_or_quad(i)
   end do ! i=1,maxele

	
	! jw, maybe we can put "xn_at_face_center", "y_face", and "h_face" into a single variable: 
	! face_center(x,y,h) which has face_center(maxface,3)
	write(pw_dia_geometry,*) '4. check_geo(4)'
	write(pw_dia_geometry,*) '(4.1) total face number:'
	write(pw_dia_geometry,*) 'maxface = ', maxface	
	write(pw_dia_geometry,*)	
	write(pw_dia_geometry,*) '(4.2) face constructing cells, nodes, face length, face center coordinates, depth at face:'
	write(pw_dia_geometry,*) 'adj_cellnum_at_face(1,j), adj_cellnum_at_face(2,j), &
	&	nodenum_at_face(1,j), nodenum_at_face(2,j), face_length(j), x_face(j), y_face(j), h_face(j)'
   do j = 1, maxface
      ! write(pw_dia_geometry,'(A8, 5I10, 4F15.5)') &
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

	! (5) calculation of sign function associated with the orientation of the normal velocity
	! write(*,*) 'i, j, jsj, adj_cellnum_at_face(1,jsj), adj_cellnum_at_face(2,jsj), sign_in_outflow(j,i)'
	! allocate & initialize global variables ==================================!	
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

	! (6) compute center of each element & h_cell =============================!
   allocate(x_cell(maxele), &
   &			y_cell(maxele), &
   &			h_cell(maxele))		
   x_cell 	= 0.0_dp
   y_cell 	= 0.0_dp
   h_cell 	= -5.0e8	! just setup this initial value as a big negative number; this will be updated later.
	
	
   do i = 1, maxele
   	! The following two forms will give identical results...
   	
   	! this is the more general form ----------------------------------------!
   	x_sum = 0.0_dp
   	y_sum = 0.0_dp
      
      ! x_cell(i) = (xn1 + xn2 + xn3)/3 for triangular element, where xn is the nodal x-point
      ! y_cell(i) = (yn1 + yn2 + ny3)/3 for triangular element, wehre yn is the nodal y-point
      do l = 1, tri_or_quad(i)
      	x_sum = x_sum + x_node(nodenum_at_cell(l,i))
      	y_sum = y_sum + y_node(nodenum_at_cell(l,i))

         ! note: we are using cell depth as a depth of botoom most face depth
         if(h_face(facenum_at_cell(l,i)) > h_cell(i)) then
            h_cell(i) = h_face(facenum_at_cell(l,i))
         end if
      end do
      x_cell(i) = x_sum/tri_or_quad(i)
      y_cell(i) = y_sum/tri_or_quad(i)

		
		! this is the more compact form ----------------------------------------!
! 		do l = 1, tri_or_quad(i)
! 			x_cell(i) = x_cell(i) + x_node(nodenum_at_cell(l,i)) / tri_or_quad(i)
! 			y_cell(i) = y_cell(i) + y_node(nodenum_at_cell(l,i)) / tri_or_quad(i)
! 			if(h_face(facenum_at_cell(l,i)) > h_cell(i)) then
! 				h_cell(i) = h_face(facenum_at_cell(l,i))
! 			end if
! 		end do
   end do
	



	write(pw_dia_geometry,*) '6. check_geo(6): Cell center inforamtion'
	write(pw_dia_geometry,*) 'Note: if Voronoi == 1, voronoi center information will be used instead of the following information'
	write(pw_dia_geometry,*) 'x_cell(i), y_cell(i), h_cell(i)'
	do i=1,maxele
		! write(pw_dia_geometry, '(I8, 3F18.7)') i, x_cell(i), y_cell(i), h_cell(i)
		write(pw_dia_geometry,'(A9,I10,A1, 3E20.10)') &
		&	' Cell# = ', i, ', ', x_cell(i), y_cell(i), h_cell(i)
	end do
	write(pw_dia_geometry,*)

	! (7) calculate voronoi center of elements & face angle ===================!   
	! 		x_cell(i), y_cell(i), x_face(j), y_face(j)									-> update if voronoi is selected
	!  	delta_j(j), cos_theta(j), sin_theta(j), cos_theta2(j), sin_theta2(j) -> new calculation
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
	
	! calculate delta_j(j) depending on (1) voronoi center or (2) geometric center
   if(voronoi == 1)then
   	! use voronoi center as a cell center, otherwise, use centroid of a triangle (or rectangle)
		! calculate_voronoi_center.f90 will update following variables using the Voronoi center information
		! 		x_cell(i), y_cell(i)					-> until now, it was the centroid of an element
		! 		x_face(j), y_face(j)					-> until now, it was at the exact half position at each face
		! 		delta_j(j)								-> new calculation
      call calculate_voronoi_center(xn1_at_face, xn2_at_face, yn1_at_face, yn2_at_face)
   else
		! compute angles (orientation of local x) and perpendicular distances for sides
      do j = 1, maxface
         if(adj_cellnum_at_face(2,j) /= 0) then
				! compute delta_j(j) = distance between the centers of two adjacent elements
				! delta_j for the non Voronoi option has not yet been calculated.
				! So, now calculate it here:
				! note: if "voronoi" is activated, it is calculated in [calculated_voronoi_center.f90]
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
         
         ! I can use the following two approaches here, but the second one will not work since "boundary_type_of_face" has not yet been defined.
         ! so, let's just use the first approach even though this will include land boundary also...
			if(adj_cellnum_at_face(2,j) == 0) then ! this approach will include land boundary also
         ! if(boundary_type_of_face(j) > 0) then ! while this approach will include open boundary only, and that is why this is better approach
         	! compute delta_j(j) at boundary face:
         	! 		calculate distance between the center of the boundary face and the center of the boundary element,
         	! 		then multiply by 2 (for artificial ghost cell center)
         	! 		i.e., I am calculating the distance between the current cell center and the center of the ghost cell at the open boundary face
         	delta_j(j) = dsqrt( &
         	&	(x_cell(adj_cellnum_at_face(1,j)) - x_face(j))**2 + &
         	&	(y_cell(adj_cellnum_at_face(1,j)) - y_face(j))**2) &
         	&	* 2.0         	
         end if
      end do         
   end if
   
   ! now calculate x and y unit normal vector of unit normal direction =======!
   ! I can put following part into the previous do loop, but I didn't.
   ! There is a reason why I separately calculate this loop again.
   ! If I want to pu this part into the previous do loop, I have to keep angle_theta0 and angle_theta1 as an array,
   ! then they will be very big. That is why I separately calculate the following part.
   do j=1,maxface
		! note: this is for the true face normal axis
		! so this angle results will be identical to the values which calculated in calculate_voronoi_center.f90
      angle_theta0 = datan2(xn1_at_face(j) - xn2_at_face(j), yn2_at_face(j) - yn1_at_face(j)) ! angle between the original Cartesian x-axis vs the face-normal x'-axis
      cos_theta(j) = dcos(angle_theta0)	! face normal velocity * cos_theta = true east (u) velocity
      sin_theta(j) = dsin(angle_theta0)	! face normal velocity * sin_theta = true north (v) velocity
      
      ! calculate angle to the fake face normal axis: ---------------------!
      ! this is the additional part from v58
      if(adj_cellnum_at_face(2,j) == 0) then
      	angle_theta1 = angle_theta0
      else
      	! if voronoi center is used, angle_theta1 == angle_theta0,
      	! if geometric ecnter is used, angle_theta1 /= angle_theta0.
         x1 = x_cell(adj_cellnum_at_face(1,j))
         y1 = y_cell(adj_cellnum_at_face(1,j))
         
         x2 = x_cell(adj_cellnum_at_face(2,j))
         y2 = y_cell(adj_cellnum_at_face(2,j))
      	
         ! note these two points are on the x-axis line, thus I don't need to transform, and so use them directly.
         angle_theta1 = datan2(y2-y1,x2-x1) ! this is the angle of the fake face-normal axis
      end if
      
      ! angle from true face-normal to fake face-normal axis in the counterclockwise direction.
      ! thus, if it is positive, the fake face-normal axis is on the countclockwise direction.
      angle_theta2 = (angle_theta1 - angle_theta0)
      ! write(*,*) j, angle_theta0 * 180.0/pi, angle_theta1*180.0/pi, angle_theta2 * 180.0/pi
      
      ! This is for Dr. Sun's approach
      ! cos_theta2(j) = dcos(angle_theta2 + 0.5_dp * pi) ! alpha = 90 + theta2
      ! sin_theta2(j) = dsin(angle_theta2 + 0.5_dp * pi) ! alpha = 90 + theta2
      
      ! This is jw's approach
      cos_theta2(j) = dcos(angle_theta2)
      sin_theta2(j) = dsin(angle_theta2)
	end do

   
	write(pw_dia_geometry,*) '7. check_geo(7): unit normal vector at each face'
	write(pw_dia_geometry,*) 'cos_theta(j), cos_theta2(j), sin_theta(j), sin_theta(2), delta_j(j), '
	do j=1,maxface
		! write(pw_dia_geometry,'(A22, I8, E15.5, E15.5, E15.5, E15.5, F15.5)') &
		write(pw_dia_geometry,'(A9,I10,A1, 5E20.10)') &
		&	' Face# = ', j, ',', cos_theta(j), cos_theta2(j), sin_theta(j), sin_theta2(j), delta_j(j)
	end do
	write(pw_dia_geometry,*)
end subroutine set_geometry_1
