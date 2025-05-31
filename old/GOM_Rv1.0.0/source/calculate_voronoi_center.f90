!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Calculate voronoi center of elements
!! This will calculate Voronoi center, and then update (substitute) following variables 
!! which is already calculated with geometric center information:
!! 	x_cell(i), y_cell(i)
!! 	x_face(j), y_face(j), delta_j(j)
!!
subroutine calculate_voronoi_center(xn1_at_face, xn2_at_face, yn1_at_face, yn2_at_face)
   use mod_global_variables 
   use mod_file_definition

   implicit none

   real(dp),dimension(maxface),intent(in) :: xn1_at_face, xn2_at_face, yn1_at_face, yn2_at_face

   integer :: i, j, l, i_node_no, off_element_voronoi_center
   real(dp):: temp_x, temp_y, xi, xj, xk, xm, yi, yj, yk, ym, xij, yij, xjk, yjk, xkl, ykl
   real(dp):: numerator, denominator, x1, x2, x3, x4, y1, y2, y3, y4, ua, ub
   ! jw

	! jw
   real(dp),dimension(maxele) :: voronoi_center_x, voronoi_center_y
   real(dp),dimension(maxface):: voronoi_face_center_x, voronoi_face_center_y, delta_j_voronoi
   ! jw

	! jw
   open(pw_voronoi_center,file=id_voronoi_center,form='formatted',status='replace')
   write(pw_voronoi_center,'(a)') 'Title = "Voronoi center of element"'
   ! jw
   write(pw_voronoi_center,'(a)') 'i,     x,     y'
   
	! jw
   do i = 1, maxele
      if(tri_or_quad(i) == 3)then
         xi = x_node(nodenum_at_cell(1,i))
         yi = y_node(nodenum_at_cell(1,i))
         xj = x_node(nodenum_at_cell(2,i))
         yj = y_node(nodenum_at_cell(2,i))
         xk = x_node(nodenum_at_cell(3,i))
         yk = y_node(nodenum_at_cell(3,i))
      
         xij = 0.5_dp*(xi+xj)
         yij = 0.5_dp*(yi+yj)
         xjk = 0.5_dp*(xj+xk)
         yjk = 0.5_dp*(yj+yk)

         numerator   = (xij-xjk)*(xj-xi) + (yij-yjk)*(yj-yi)
         denominator = (xj -xi) *(yk-yj) - (xk -xj) *(yj-yi)

			! jw
         if(numerator*denominator >= 0.0_dp) then
            off_element_voronoi_center = off_element_voronoi_center + 1
         end if

         if(denominator /= 0.0) then
            voronoi_center_x(i) = xjk + numerator/denominator*(yk-yj)
            voronoi_center_y(i) = yjk - numerator/denominator*(xk-xj)
         else
            write(pw_run_log,*) 'denominator of voronoi is zero at', i
            stop 'calculate_voronoi_center.f90, Error #1'
         end if
      else if(tri_or_quad(i) == 4) then
         xi = x_node(nodenum_at_cell(1,i))
         yi = y_node(nodenum_at_cell(1,i))
         xj = x_node(nodenum_at_cell(2,i))
         yj = y_node(nodenum_at_cell(2,i))
         xk = x_node(nodenum_at_cell(3,i))
         yk = y_node(nodenum_at_cell(3,i))
         xm = x_node(nodenum_at_cell(4,i))
         ym = y_node(nodenum_at_cell(4,i))

         xij = 0.5_dp*(xi+xj)
         yij = 0.5_dp*(yi+yj)
         xjk = 0.5_dp*(xj+xk)
         yjk = 0.5_dp*(yj+yk)
         xkl = 0.5_dp*(xk+xm)
         ykl = 0.5_dp*(yk+ym)

         numerator   = (xij-xjk)*(xj-xi) + (yij-yjk)*(yj-yi)
         denominator = (xj -xi) *(yk-yj) - (xk -xj) *(yj-yi)


         temp_x = 0.0_dp
         temp_y = 0.0_dp
         do l = 1, tri_or_quad(i)
            i_node_no = nodenum_at_cell(l,i)
            temp_x = temp_x + x_node(i_node_no)
            temp_y = temp_y + y_node(i_node_no)
         end do
         voronoi_center_x(i) = temp_x / 4.0_dp
         voronoi_center_y(i) = temp_y / 4.0_dp
      end if
   end do
   
   ! jw
   do i=1,maxele
   	write(pw_voronoi_center,'(I8, 2f18.7)') i, voronoi_center_x(i), voronoi_center_y(i)
   end do
   close(pw_voronoi_center)
   
   ! jw
	! jw
	! jw
	! jw
	! jw
	open(pw_voronoi_face_center,file=id_voronoi_face_center,form='formatted',status='replace')
   write(pw_voronoi_face_center,'(a)') 'Title = "Voronoi face center of element"'
   write(pw_voronoi_face_center,'(a)') 'Variables = "x", "y", "x_normal_vector", "y_normal_vector", "face number"'

   do j = 1, maxface
      if(adj_cellnum_at_face(2,j) /= 0) then
			! jw
         delta_j_voronoi(j) = dsqrt (   &
                 &           (voronoi_center_x(adj_cellnum_at_face(2,j))       &
                 &         -  voronoi_center_x(adj_cellnum_at_face(1,j)))**2   &
                 &         + (voronoi_center_y(adj_cellnum_at_face(2,j))       &
                 &         -  voronoi_center_y(adj_cellnum_at_face(1,j)))**2)
         if(delta_j_voronoi(j) == 0.0_dp) then
            write(pw_run_log,*) 'zero distance between voronoi centers', adj_cellnum_at_face(1,j), adj_cellnum_at_face(2,j)
            stop 'calculate_voronoi_center.f90, Error #2'
         end if

			! jw
	      ! jw
	      ! jw
	      ! jw
         ! jw
         ! jw
         ! jw

			! jw
         x1 = voronoi_center_x(adj_cellnum_at_face(1,j))
         y1 = voronoi_center_y(adj_cellnum_at_face(1,j))
         x2 = voronoi_center_x(adj_cellnum_at_face(2,j))
         y2 = voronoi_center_y(adj_cellnum_at_face(2,j))

         x3 = xn2_at_face(j)
         y3 = yn2_at_face(j)
         x4 = xn1_at_face(j)
         y4 = yn1_at_face(j)

         denominator = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)

         ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denominator
         ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denominator
  			
  			! jw
  			! jw
         ! jw
         ! jw
         
         ! jw
         voronoi_face_center_x(j) = x_face(j)
         voronoi_face_center_y(j) = y_face(j)
         

         write(pw_voronoi_face_center,'(I8, 2F18.7)') j, voronoi_face_center_x(j), voronoi_face_center_y(j)
      end if      
   end do
   close(pw_voronoi_face_center)
	! jw
	
	
	! jw
	! jw
	! jw
   do i = 1, maxele
      x_cell(i) = voronoi_center_x(i)
      y_cell(i) = voronoi_center_y(i)
   end do
   
   ! jw
   ! jw
   do j = 1, maxface
   	if(adj_cellnum_at_face(2,j) /= 0) then
   		! jw
   		! jw
   		delta_j(j) = delta_j_voronoi(j)
	      x_face(j) = voronoi_face_center_x(j)
	      y_face(j) = voronoi_face_center_y(j)      
   		
      ! jw
      else if(adj_cellnum_at_face(2,j) == 0) then
      	! jw
      	! jw
      	! jw
      	! jw
      	! jw
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
   ! jw
   
end subroutine calculate_voronoi_center 
