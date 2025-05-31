!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Compute temporary velocitiy at boundary,
!! and calculate interpolated river discharge calling 'calculate_Q.f90'.
!! 
subroutine compute_boundary_velocity
   use mod_global_variables
	use mod_file_definition
   implicit none
   
   integer :: i, j, k
   integer :: n1, n2, nob
   real(dp):: temp_depth
   real(dp),allocatable,dimension(:) :: channel_width, channel_width_wr
   ! jw
   
   ! jw
   u_boundary = 0.0_dp
   v_boundary = 0.0_dp
   
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
	if(num_Qb_cell > 0) then
	   allocate(channel_width(num_Qb_cell))
		channel_width = 0.0_dp
		
		! jw
		call calculate_Q
		
		! jw
		! jw
		! jw
		! jw
		do i=1,num_Qb_cell
			! jw
			channel_width(i) = face_length(Q_boundary(i,2)) ! jw
			
			j = Q_boundary(i,2) ! jw
			n1 = nodenum_at_face(1,j) ! jw
			n2 = nodenum_at_face(2,j) ! jw
			
			temp_depth   = h_face(j) + (eta_node(n1)+eta_node(n2))*0.5_dp
	      if(temp_depth < dry_depth) then
	         write(pw_run_log,*) 'Error: River boundary is dry at Qb_cell#: ', i
	         stop
	      end if

	      ! jw
	      ! jw
	      ! jw
	      ! jw
	      ! jw
	      Qu_boundary(i) = -Q_add(i)/(channel_width(i) * temp_depth) ! jw
	      ! jw
		end do
		deallocate(channel_width)
	end if
	
	! jw
	! jw
	if(num_WR_cell > 0) then
		allocate(channel_width_wr(num_WR_cell))
		channel_width_wr = 0.0_dp
		
		! jw
		call calculate_Q_WR
		
		! jw
		do i=1,num_WR_cell
			! jw
			j = WR_boundary(i,2) ! jw
			n1 = nodenum_at_face(1,j) ! jw
			n2 = nodenum_at_face(2,j) ! jw

			! jw
			channel_width_wr(i) = face_length(j)
			
			! jw
			if(WR_layer(i) == 999) then
				! jw
				k = top_layer_at_face(j)
			else if(WR_layer(i) == 0) then
				! jw
				k = bottom_layer_at_face(j)
			else
				! jw
				if((WR_layer(i) > top_layer_at_face(j)) .or. (WR_layer(i) < bottom_layer_at_face(j))) then
					write(pw_run_log,*) 'Either WR_layer(i) > top_layer_at_face or WR_layer(i) < bottom_layer_at_face'
					write(*,*) 'Either WR_layer(i) > top_layer_at_face or WR_layer(i) < bottom_layer_at_face'
					stop
				end if
				
				k = WR_layer(i)
			end if
			temp_depth = dz_face(k,j)
	      if(temp_depth < dry_depth) then
	         write(pw_run_log,*)  'Error: Writhdraw/Return boundary is dry at WR_cell#: ', i
	         write(*,*)  'Error: Writhdraw/Return boundary is dry at WR_cell#: ', i
	         stop
	      end if

	      ! jw
	      ! jw
	      ! jw
	      ! jw
	      ! jw
	      WRu_boundary(i) = -Q_add_WR(i)/(channel_width_wr(i) * temp_depth) ! jw
	      ! jw
		end do
		deallocate(channel_width_wr)
	end if

	! jw
	! jw
	! jw
	! jw
	! jw
 	if(num_SS_cell > 0) then		
 		! jw
		call calculate_Q_SS
	end if


	! jw
   ! jw
   
end subroutine compute_boundary_velocity
