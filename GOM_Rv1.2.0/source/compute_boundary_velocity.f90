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
   ! End of local variables ==================================================!
   
   ! re-initialize following global variables:   
   u_boundary = 0.0_dp
   v_boundary = 0.0_dp
   
   ! allocate local variables:
!   allocate(channel_width(num_tide_bc)) ! total channel width at the boundary segment 
! 	channel_width = 0.0_dp
   
	
   ! This is for tidal river boundary ========================================!
	! calculate tidal river boundary channel width
!   do i = 1, num_ob_cell
!      if(ob_flow_type(i) /= 0) then
!         channel_width(i) = 0.0
!         do j = 1, num_serial_ob_node(i)
!            channel_width(i) = channel_width(i) + face_length(ob_facenum(i,j)) ! total channel width at the boundary segment 
!         end do
!      end if
!   end do
   
	! calculate velocity at tidal river boundary (this is not a true river)
	! jw, This might be worng, check this again...
!   do i = 1, num_ob_cell
!      if(ob_flow_type(i) /= 0) then
!         nob = ob_facenum(i,1)
!         n1   = nodenum_at_face(1,nob)
!         n2   = nodenum_at_face(2,nob)
!         temp_depth   = h_face(nob) + 0.5*(eta_node(n1)+eta_node(n2)) ! jw, H = h + eta
!         if(temp_depth < dry_depth) then
!            write(pw_run_log,*)'this is dry boundary side at',n1,n2,temp_depth
!            stop 
!         end if
!         u_boundary(i) = tide_Q(i)*spinup_function_tide/(channel_width(i) * temp_depth) * cos_theta(nob) ! (-cos_theta(nob)) ???, jw, check this again
!         v_boundary(i) = tide_Q(i)*spinup_function_tide/(channel_width(i) * temp_depth) * sin_theta(nob)
!      end if
!   end do
!   deallocate(channel_width)
   ! End of tidal river boundary =============================================!
   
   ! This is for real river boundary =========================================!
	! compute temporary velocity at river boundary face
	if(num_Qb_cell > 0) then
	   allocate(channel_width(num_Qb_cell))
		channel_width = 0.0_dp
		
		! This will calculate "Q_add" [m3/s] at each river cell at this time step (using interpolation)
		call calculate_Q
		
		! Now, calculate velocity at each river face at this time step
		! Actually, we don't need to calculate river bounary velocity.
		! We will use just interpolated Q_add [m3/s] in the continuity equation (solve_free_surface_equation.f90).
		! However, I will just keep this...
		do i=1,num_Qb_cell
			! calculate channal width (i.e., face length at river boundary face
			channel_width(i) = face_length(Q_boundary(i,2)) ! Q_boundary(i,1) = element_id, Q_boundary(i,2) = face_id
			
			j = Q_boundary(i,2) ! face number of the river boundary
			n1 = nodenum_at_face(1,j) ! this should equal to Q_node1(i) or Q_node2(i)
			n2 = nodenum_at_face(2,j) ! this should equal to Q_node2(i) or Q_node2(i)
			
			temp_depth   = h_face(j) + (eta_node(n1)+eta_node(n2))*0.5_dp
	      if(temp_depth < dry_depth) then
	         write(pw_run_log,*) 'Error: River boundary is dry at Qb_cell#: ', i
	         stop
	      end if

	      ! we only have to consider Qu_boundary (i.e., face normal velocity) since we are assuming river discharge is perpendicular to the river boundary face
	      ! Here, Qu_boundary is u* value, not u value in the true xy coordinate
	      ! Actually, it is better (clear) to use [sign_in_outflow] to make sure to define the positive direction at this face,
	      ! however, it doesn't required since all boundary face has positive direction toward the outer domain direction.
	      ! That is why I just put (-) sign here.
	      Qu_boundary(i) = -Q_add(i)/(channel_width(i) * temp_depth) ! since positive Q_add should point inside the cell
	      ! write(*,*) 1, Qu_boundary(1)
		end do
		deallocate(channel_width)
	end if
	
	! This is for Withdraw/Return boundary ====================================!
	! compute temporary velocity at WR boundary face
	if(num_WR_cell > 0) then
		allocate(channel_width_wr(num_WR_cell))
		channel_width_wr = 0.0_dp
		
		! This will calculate "Q_add_WR" [m3/s] at each Withdraw/Return cell at this time step (using interpolation)
		call calculate_Q_WR
		
		! Now, calculate velocity at each WR face at this time step
		do i=1,num_WR_cell
			! WR_boundary(i,1) = element_id, WR_boundary(i,2) = face_id
			j = WR_boundary(i,2) ! face number of the WR boundary
			n1 = nodenum_at_face(1,j) ! this should equal to WR_node1(i) or WR_node2(i)
			n2 = nodenum_at_face(2,j) ! this should equal to WR_node2(i) or WR_node2(i)

			! calculate channal width (i.e., face length at river boundary face
			channel_width_wr(i) = face_length(j)
			
			! temp_depth   = h_face(nob) + (eta_node(n1)+eta_node(n2))*0.5_dp
			if(WR_layer(i) == 999) then
				! anytime surface layer
				k = top_layer_at_face(j)
			else if(WR_layer(i) == 0) then
				! anytime bottom layer
				k = bottom_layer_at_face(j)
			else
				! specific vertical layer
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

	      ! we only have to consider WRu_boundary (i.e., face normal velocity) since we are assuming Withdraw/Return discharge is perpendicular to the WR boundary face
	      ! Here, WRu_boundary is u* value, not u value in the true xy coordinate
	      ! Actually, it is better (clear) to use [sign_in_outflow] to make sure to define the positive direction at this face,
	      ! however, it doesn't required since all boundary face has positive direction toward the outer domain direction.
	      ! That is why I just put (-) sign here.
	      WRu_boundary(i) = -Q_add_WR(i)/(channel_width_wr(i) * temp_depth) ! since positive Q_add_WR should point inside the cell
	      ! write(*,*) i, Q_add_WR(i), channel_width_wr(i), temp_depth, WRu_boundary(i)
		end do
		deallocate(channel_width_wr)
	end if

	! This is for Source/Sink boundary ========================================!
	! This was the original approach, but I will set this option as a real source/sink, 
	! that only takes account the total amount of concentration not the flow.
	! So, I will kill all of these approach.
	! compute temporary velocity at SS boundary cell (not at face but at cell...)
 	if(num_SS_cell > 0) then		
 		! This will calculate "Q_add_SS" [m3/s] at each Source/Sink cell at this time step (using interpolation)
		call calculate_Q_SS
	end if


	! stop 'I am here: compute_boundary_velocity'
   ! End of real river boundary ==============================================!
   
end subroutine compute_boundary_velocity
