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
	! jw
   ! jw
   ! jw
   ! jw
   u_face = 0.0
   v_face = 0.0
   
   	
	! jw
	! jw
 	!$omp parallel
 	!$omp do private(j,k,u_face,v_face)
	do j = 1, maxface
		if(top_layer_at_face(j) /= 0 ) then	! jw
			! jw
			! jw
			! jw
			! jw
			
			! jw
			do k = bottom_layer_at_face(j), top_layer_at_face(j) 
				! jw
				! jw
         	! jw
         	! jw
				u_face(k) = un_face_new(k,j) * cos_theta(j) - vn_face_new(k,j) * sin_theta(j) ! jw
				v_face(k) = un_face_new(k,j) * sin_theta(j) + vn_face_new(k,j) * cos_theta(j) ! jw
			end do
			
			! jw
			u_face_level(bottom_layer_at_face(j)-1,j) = u_face(bottom_layer_at_face(j))
			v_face_level(bottom_layer_at_face(j)-1,j) = v_face(bottom_layer_at_face(j))
			u_face_level(top_layer_at_face(j),j)      = u_face(top_layer_at_face(j))
			v_face_level(top_layer_at_face(j),j)      = v_face(top_layer_at_face(j))
			
			! jw
			do k = bottom_layer_at_face(j), top_layer_at_face(j)-1 !m > m
				u_face_level(k,j) = u_face(k) + dz_face(k,j)*0.5_dp/dzhalf_face(k,j)*(u_face(k+1)-u_face(k))
				v_face_level(k,j) = v_face(k) + dz_face(k,j)*0.5_dp/dzhalf_face(k,j)*(v_face(k+1)-v_face(k))
			end do
		end if
	end do
 	!$omp end do
	
	! jw
	! jw
	! jw
	! jw
 	!$omp do private(i,j,k,l,weight,icount,isn,nc,ifn,kvt,factor_weight)
   do i = 1, maxnod
      do k = 0, maxlayer
         u_node(k,i) = 0.0_dp
         v_node(k,i) = 0.0_dp
         weight = 0.0_dp
         icount = 0
         do j = 1, adj_cells_at_node(i)         	
            isn = adj_cellnum_at_node(j,i) ! jw
            nc  = node_count_each_element(j,i) ! jw

            do l = tri_or_quad(isn)-2, tri_or_quad(isn)-1 ! jw
            	! jw
               ifn = facenum_at_cell(start_end_node(tri_or_quad(isn),nc,l),isn)
               
               ! jw
               ! jw
               ! jw
               ! jw
               ! jw
               if(adj_cellnum_at_face(2,ifn) == 0)then
                  factor_weight = 2.0_dp
               else
                  factor_weight = 1.0_dp
               end if
               
               ! jw
               ! jw
               if(top_layer_at_face(ifn) /= 0 .and. k >= bottom_layer_at_face(ifn)-1) then
                  kvt = MIN(k,top_layer_at_face(ifn))
                  u_node(k,i) = u_node(k,i)   &
                  &           + u_face_level(kvt,ifn) / face_length(ifn)*factor_weight
                  v_node(k,i) = v_node(k,i)   &
                  &           + v_face_level(kvt,ifn) / face_length(ifn)*factor_weight
                  icount = icount + 1
               end if
               weight = weight + factor_weight/face_length(ifn)
               
               ! jw
            end do
         end do
			
			! jw
         if(icount /= 0) then
            u_node(k,i) = u_node(k,i)/weight
            v_node(k,i) = v_node(k,i)/weight
         end if
      end do ! jw
   end do ! jw
 	!$omp end do nowait
	
	! jw
	! jw
 	!$omp do private(i,j,k,total_weight,ite,kvt)
   do i = 1, maxnod
      do k = 0, maxlayer
         w_node(k,i) = 0.0_dp
         total_weight = 0.0_dp  ! jw
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
   

	! jw
	if(dia_node_velocity == 1) then
		write(pw_dia_node_velocity,*) 'it = ', it, ', elapsed_time = ', elapsed_time		
   	write(pw_dia_node_velocity,'(A)') 'u_node(1,i), v_node(1,i), w_node(1,i)'
   	do i=1,maxnod
	   	write(pw_dia_node_velocity,'(A3, I5, 3F10.5)') 'i=', i, u_node(1,i), v_node(1,i), w_node(1,i)
   	end do	
   end if
   
end subroutine calculate_velocity_at_node
