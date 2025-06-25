!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! We have already calculated velocities at face center in calculate_velocity_at_face.f90, and now
!! calculate velocities at nodes.
!! jw:
!! 	[i,j,k] order has been updated.
!! 	[division to multiplication] has been updated.
subroutine calculate_velocity_at_node
   use mod_global_variables
   use mod_file_definition
   implicit none
   
   integer :: i, j, k, l
   integer :: icount, isn, nc, ifn, kvt, ite
   real(dp):: weight, factor_weight, total_weight
   real(dp),dimension(0:maxlayer+1) ::	u_face, v_face
	! Note: u_face_level & v_face_level are defined at mod_global_variables since they have [maxface ~] arrays
   ! &		u_face_level,								&	! (0:maxlayer+1,maxface)
   ! &		v_face_level									! (0:maxlayer+1,maxface)
   ! End of local variables ==================================================!
   u_face = 0.0
   v_face = 0.0
   
   	
	! First step:
	! 	calculate face velocities (true east-north velocities) at whole levels (k'th level)
 	!$omp parallel
 	!$omp do private(j,k,u_face,v_face)
	do j = 1, maxface
		if(top_layer_at_face(j) /= 0 ) then	! if this face is wet
			! Note, following to do loops can be simplified as (but I didn't):
			! do k = bottom_layer_at_face(j), top_layer_at_face(j)-1
			! then...
			! do only for top_layer_at_face(j)
			
			! calculate consultane velocities	   
			do k = bottom_layer_at_face(j), top_layer_at_face(j) 
				! Here, u & v are true velocities at true east-north coordinate, and
				!       U & V are face normal velocities (i.e., each face's coordinate)
         	! u = U*cos - V*sin, 
         	! v = U*sin + V*cos
				u_face(k) = un_face_new(k,j) * cos_theta(j) - vn_face_new(k,j) * sin_theta(j) ! u_face(k) is the true east velocity at k'th vertical layer at current face
				v_face(k) = un_face_new(k,j) * sin_theta(j) + vn_face_new(k,j) * cos_theta(j) ! v_face(k) is the true north velocity at k'th vertical layer at current face
			end do
			
			! here, u_face_level & v_face_level are true velocities at k'th level
			u_face_level(bottom_layer_at_face(j)-1,j) = u_face(bottom_layer_at_face(j))
			v_face_level(bottom_layer_at_face(j)-1,j) = v_face(bottom_layer_at_face(j))
			u_face_level(top_layer_at_face(j),j)      = u_face(top_layer_at_face(j))
			v_face_level(top_layer_at_face(j),j)      = v_face(top_layer_at_face(j))
			
			! find true velocities at k'th level using linear interpolation
			do k = bottom_layer_at_face(j), top_layer_at_face(j)-1 !m > m
				u_face_level(k,j) = u_face(k) + dz_face(k,j)*0.5_dp/dzhalf_face(k,j)*(u_face(k+1)-u_face(k))
				v_face_level(k,j) = v_face(k) + dz_face(k,j)*0.5_dp/dzhalf_face(k,j)*(v_face(k+1)-v_face(k))
			end do
		end if
	end do
 	!$omp end do
	
	! Second step:
	! 	calculate horizontal velocities (true east-north velocities) at nodes using u_face_level & v_face_level (face velocities at each vertical level)
	! Here, I am using IDW (Inverse Distance Weighting interpolation) with power = 1.
	! If I have pre-calculated adj_face_at_node information, this step will be much easier.
 	!$omp do private(i,j,k,l,weight,icount,isn,nc,ifn,kvt,factor_weight)
   do i = 1, maxnod
      do k = 0, maxlayer
         u_node(k,i) = 0.0_dp
         v_node(k,i) = 0.0_dp
         weight = 0.0_dp
         icount = 0
         do j = 1, adj_cells_at_node(i)         	
            isn = adj_cellnum_at_node(j,i) ! cell id of neighbor cell
            nc  = node_count_each_element(j,i) ! nc'th node in this neighbor cell

            do l = tri_or_quad(isn)-2, tri_or_quad(isn)-1 ! l = 1,2 for tri_element, l=2,3 for quad_element
            	! ifn is the face_id which connected with this node, and this is what I want to find
               ifn = facenum_at_cell(start_end_node(tri_or_quad(isn),nc,l),isn)
               
               ! if this face has two neighboring cells, this face velocity value will be calculated twice.
               ! But, if this face has only one neighboring cell, this face velocity will be used just once.
               ! And, that is unfare. Thus, it should be counted as twice.
               ! That is what factor_weight means here.
               ! With this approach each face velocity will be used twice, and it is fare now.
               if(adj_cellnum_at_face(2,ifn) == 0)then
                  factor_weight = 2.0_dp
               else
                  factor_weight = 1.0_dp
               end if
               
               ! Here, u_node is the numerator term in the IDW
               ! and weight is the denominator term in the IDW
               if(top_layer_at_face(ifn) /= 0 .and. k >= bottom_layer_at_face(ifn)-1) then
                  kvt = MIN(k,top_layer_at_face(ifn))
                  u_node(k,i) = u_node(k,i)   &
                  &           + u_face_level(kvt,ifn) / face_length(ifn)*factor_weight
                  v_node(k,i) = v_node(k,i)   &
                  &           + v_face_level(kvt,ifn) / face_length(ifn)*factor_weight
                  icount = icount + 1
               end if
               weight = weight + factor_weight/face_length(ifn)
               
               ! write(*,*) 'factor_weight=', factor_weight, 'weight=', weight
            end do
         end do
			
			! This is the final IDW results: IDW_result = numerator/denominator
         if(icount /= 0) then
            u_node(k,i) = u_node(k,i)/weight
            v_node(k,i) = v_node(k,i)/weight
         end if
      end do ! k=0,maxlayer
   end do ! i=1,maxnod
 	!$omp end do nowait
	
	! Third step:
	! 	calculate w at each node by weighted average
 	!$omp do private(i,j,k,total_weight,ite,kvt)
   do i = 1, maxnod
      do k = 0, maxlayer
         w_node(k,i) = 0.0_dp
         total_weight = 0.0_dp  ! sum of weights
         do j = 1, adj_cells_at_node(i)
            ite = adj_cellnum_at_node(j,i)
            if((top_layer_at_element(ite) /= 0) .and. (k >= (bottom_layer_at_element(ite)-1))) then
               kvt  = min(k,top_layer_at_element(ite))
               total_weight = total_weight + area(ite)
               w_node(k,i) = w_node(k,i) + wn_cell(kvt,ite) * area(ite)
            end if
         end do

         if(total_weight /= 0.0_dp) then
         	w_node(k,i) = w_node(k,i)/total_weight
         end if
      end do
   end do
 	!$omp end do
 	!$omp end parallel
   

	! write diagonostic file ==================================================!
	if(dia_node_velocity == 1) then
		write(pw_dia_node_velocity,*) 'it = ', it, ', elapsed_time = ', elapsed_time		
   	write(pw_dia_node_velocity,'(A)') 'u_node(1,i), v_node(1,i), w_node(1,i)'
   	do i=1,maxnod
	   	write(pw_dia_node_velocity,'(A3, I5, 3F10.5)') 'i=', i, u_node(1,i), v_node(1,i), w_node(1,i)
   	end do	
   end if
   
end subroutine calculate_velocity_at_node
