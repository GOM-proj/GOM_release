!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This will update vertical layer indices and wetting and drying
!! This subroutine is identical to [update_vertical_layers_0.f90] but for new variables (for it > 0).
!! Note: an element will be dry if all nodes or all sides are dry.
!! 
!! [update_vertical_layers_0.f90] vs. [update_vertical_layers_1.f90]:
!! 		dz_cell(k,i) 		-> dz_cell_new(k,i)
!!			dzhalf_cell(k,i)	-> dzhalf_cell_new(k,i)
!!			wn_cell(k,i)		-> wn_cell_new(k,i)
!!			dz_face(k,j)		-> dz_face_new(k,j)
!!			dzhalf_face(k,j)	-> dzhalf_face_new(k,j)
!!			dz_node(k,i)		-> dz_node_new(k,i)
!!			dzhalf_node(k,i)	-> dzhalf_node_new(k,i)
!!
subroutine update_vertical_layers_1
   use mod_global_variables
   use mod_file_definition   
   implicit none
   
   integer :: i, j, k, l 
   integer :: nsum1, nsum2
   real(dp):: total_water_depth, tmp
   real(dp):: eta_max
   integer :: top_layer_at_face_new,top_layer_at_node_new, top_layer_at_element_new
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
   !$omp parallel
   !$omp do private(i,k,total_water_depth,tmp,top_layer_at_element_new)
	do i = 1, maxele
		! jw
		top_layer_at_element_new = 0
		
		! jw
      total_water_depth = h_cell(i) + eta_cell_new(i)

      ! jw
   	! jw
      if(total_water_depth <= dry_depth) then   ! jw
         top_layer_at_element_new = 0			! jw
      else                                      ! jw
         top_layer_at_element_new = 0			! jw
         tmp = MSL + eta_cell_new(i)
         
         do k = maxlayer-1, 0, -1
            if(tmp > z_level(k) .and. tmp <= z_level(k+1)) then
               top_layer_at_element_new = k + 1
               dz_cell_new(top_layer_at_element_new,i) = tmp - z_level(k)
               exit
            end if 
         end do
         
         ! jw
     		if(top_layer_at_element_new < bottom_layer_at_element(i)) then
            if(tmp > z_level(maxlayer)) then
               write(pw_run_log,*)'large elevation at element', i, 'eta_cell_new(i)=', eta_cell_new(i)
               write(pw_run_log,*)'z_level(maxlayer)=', z_level(maxlayer)
               write(*,  *)'large elevation at element', i, 'eta_cell_new(i)=', eta_cell_new(i)
               write(*,  *)'z_level(maxlayer)=', z_level(maxlayer)
               write(*,*) 'press enter to continue'
            else
               write(pw_run_log,*)'z_level greater than layer', i, MSL+eta_cell_new(i), &
                 &         MSL - h_cell(i)
               write(*  ,*)'z_level greater than layer', i, MSL+eta_cell_new(i), &
                 &         MSL - h_cell(i)
               write(*,*) 'press enter to continue' 
            end if   !!! if(tmp > z_level(maxlayer)) then
            stop
         end if   !!! if(top_layer_at_element_new < bottom_layer_at_element(i)) then
         
         if(top_layer_at_element_new == bottom_layer_at_element(i)) then
            dz_cell_new(top_layer_at_element_new,i)     = total_water_depth
            dzhalf_cell_new(top_layer_at_element_new,i) = total_water_depth

				! jw
! jw
! jw
! jw
! jw
! jw
! jw
         else   ! jw
          	! jw
         	! jw
            dz_cell_new(bottom_layer_at_element(i),i) =    &
            &                             z_level(bottom_layer_at_element(i))   &
            &                          - (MSL-h_cell(i))

            do k = bottom_layer_at_element(i)+1, top_layer_at_element_new-1
               dz_cell_new(k,i) = delta_z(k)
            end do

            do k = bottom_layer_at_element(i), top_layer_at_element_new-1
               dzhalf_cell_new(k,i) = (dz_cell_new(k+1,i)+dz_cell_new(k,i))*0.5
            end do
	         dzhalf_cell_new(top_layer_at_element_new,i) = dz_cell_new(top_layer_at_element_new,i)
	         dzhalf_cell_new(bottom_layer_at_element(i)-1,i) = dz_cell_new(bottom_layer_at_element(i),i)            
            
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
         end if
      end if
      
      ! jw
      if(top_layer_at_element(i) == 0 .and. top_layer_at_element_new /= 0) then
         do k = 0, maxlayer
            wn_cell_new(k,i) = 0.0d0
         end do   ! jw
      end if ! jw
      
      ! jw
      top_layer_at_element(i) = top_layer_at_element_new
   end do
 	!$omp end do nowait
   
   
	! jw
 	!$omp do private(j,k,eta_max,total_water_depth,tmp,top_layer_at_face_new)
   do j = 1, maxface
   	! jw
   	top_layer_at_face_new = 0
   	
   	! jw
      if(adj_cellnum_at_face(2,j) /= 0) then
         eta_max = dmax1(eta_cell_new(adj_cellnum_at_face(1,j)), eta_cell_new(adj_cellnum_at_face(2,j)))
      else
         eta_max = eta_cell_new(adj_cellnum_at_face(1,j))
      end if

      total_water_depth = h_face(j) + eta_max

      if(total_water_depth <= dry_depth) then
         top_layer_at_face_new = 0
      else   !!!   not dry
         top_layer_at_face_new = 0
         tmp = MSL + eta_max
         
         do k = maxlayer-1, 0, -1
            if(tmp > z_level(k) .and. tmp <= z_level(k+1)) then
               top_layer_at_face_new = k + 1
               dz_face_new(top_layer_at_face_new,j) = tmp - z_level(k)
               exit
            end if
         end do 
      	
      	if(top_layer_at_face_new < bottom_layer_at_face(j)) then
            if(tmp > z_level(maxlayer)) then
               write(pw_run_log,*) 'large elevation at side', j, eta_max
               write(*,*) 'large elevation at side', j, eta_max
            else
               write(pw_run_log,*)'impossible',j,MSL+eta_max, MSL-h_face(j)
               write(*,*)'impossible',j,MSL+eta_max, MSL-h_face(j)
            end if
            write(*,*) 'press enter to continue'
            stop
          end if
			
			! jw
         if(top_layer_at_face_new == bottom_layer_at_face(j)) then
            dz_face_new(top_layer_at_face_new,j)       = total_water_depth
            dzhalf_face_new(top_layer_at_face_new-1,j) = total_water_depth
            dzhalf_face_new(top_layer_at_face_new,j)   = total_water_depth
         else 
            dz_face_new(bottom_layer_at_face(j),j) =  z_level(bottom_layer_at_face(j))   &
            &                                 - (MSL-h_face(j))
            do k = bottom_layer_at_face(j)+1, top_layer_at_face_new-1
               dz_face_new(k,j) = delta_z(k)
            end do
            do k = bottom_layer_at_face(j), top_layer_at_face_new-1
               dzhalf_face_new(k,j) = (dz_face_new(k+1,j) + dz_face_new(k,j)) * 0.5
            end do
            dzhalf_face_new(top_layer_at_face_new,j)  = dz_face_new(top_layer_at_face_new,j)
            dzhalf_face_new(bottom_layer_at_face(j)-1,j) = dz_face_new(bottom_layer_at_face(j),j)
         end if   ! jw
      end if   ! jw
      
      top_layer_at_face(j) = top_layer_at_face_new
   end do
 	!$omp end do nowait
	

	! jw
 	!$omp do private(i,k,total_water_depth,tmp,top_layer_at_node_new)
   do i = 1, maxnod
   	! jw
   	top_layer_at_node_new = 0
   	
   	! jw
      total_water_depth = h_node(i) + eta_node(i)
      
      if(total_water_depth <= dry_depth) then
         top_layer_at_node_new = 0
      else ! jw
         top_layer_at_node_new = 0
         tmp = MSL + eta_node(i)
         do k = maxlayer-1, 0, -1
            if(tmp > z_level(k) .and. tmp<=z_level(k+1)) then
               top_layer_at_node_new = k+1
               dz_node_new(top_layer_at_node_new,i) = tmp - z_level(k)
               exit
            end if
         end do   ! jw

      	if(top_layer_at_node_new < bottom_layer_at_node(i)) then
            if(tmp > z_level(maxlayer)) then
               write(pw_run_log,*)'large elevation at node', i, eta_node(i)
               write(*  ,*)'large elevation at node', i, eta_node(i)
            else
               write(pw_run_log,*)'impossible', i, MSL+eta_node(i), MSL-h_node(i)
               write(*  ,*)'impossible', i, MSL+eta_node(i), MSL-h_node(i)
            end if
            write(*,*) 'press enter to continue'
            stop
         end if

			! jw
         if(top_layer_at_node_new == bottom_layer_at_node(i)) then
            dz_node_new(top_layer_at_node_new,i)     = total_water_depth
            dzhalf_node_new(top_layer_at_node_new,i) = total_water_depth
         else 
			! jw
            dz_node_new(bottom_layer_at_node(i),i) =  z_level(bottom_layer_at_node(i))   &
            &                           				- (MSL-h_node(i))
            do k = bottom_layer_at_node(i)+1, top_layer_at_node_new-1
               dz_node_new(k,i) = delta_z(k)
            end do
            do k = bottom_layer_at_node(i), top_layer_at_node_new-1
               dzhalf_node_new(k,i) = (dz_node_new(k+1,i)+dz_node_new(k,i)) * 0.5
            end do
         end if
      end if
      
      top_layer_at_node(i) = top_layer_at_node_new
   end do
 	!$omp end do
	
   
	! jw
	! jw
 	!$omp do private(i,l,nsum1,nsum2)
   do i=1,maxele
      nsum1 = 0
      nsum2 = 0
      do l = 1, tri_or_quad(i)
         nsum1 = nsum1 + top_layer_at_node(nodenum_at_cell(l,i))
         nsum2 = nsum2 + top_layer_at_face(facenum_at_cell(l,i))
      end do

      if(nsum1 == 0) then
      	top_layer_at_element(i) = 0 
      end if

      if(nsum2 == 0 .and. top_layer_at_element(i) /= 0) then
         write(pw_run_log,*)'wet element with 3 dry sides',i
         stop
      end if
   end do   ! jw
 	!$omp end do
   
	! jw
	!$omp do private(i,k)
   do i = 1, maxele
      do k = bottom_layer_at_element(i), top_layer_at_element(i)-1 ! jw
         if(dzhalf_cell_new(k,i) < denominator_min_for_matrix) then
            write(pw_run_log,*) 'b',i,k,dzhalf_cell_new(k,i),denominator_min_for_matrix
            stop
         end if
      end do
		
! jw
! jw
! jw
! jw
! jw
   end do
	!$omp end do nowait
	
	!$omp do private(j,k)
   do j=1,maxface
      do k=bottom_layer_at_face(j),top_layer_at_face(j)-1
         if(dzhalf_face_new(k,j) < denominator_min_for_matrix) then
            write(pw_run_log,*) j,k,dzhalf_face_new(k,j),denominator_min_for_matrix
            stop
         end if
      end do

! jw
! jw
! jw
! jw
! jw
   end do
	!$omp end do nowait
   
   !$omp do private(i,k)
   do i = 1, maxnod
      do k = bottom_layer_at_node(i), top_layer_at_node(i)-1
         if(dzhalf_node_new(k,i) < denominator_min_for_matrix) then
            write(pw_run_log,*) i,k,dzhalf_node_new(k,i),denominator_min_for_matrix
            stop
         endif
      end do

! jw
! jw
! jw
! jw
! jw
   end do
	!$omp end do
	!$omp end parallel
   
end subroutine update_vertical_layers_1
