!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine update_variables
	use mod_global_variables
	implicit none
	
	integer :: i, k
	real(dp):: sum_u, sum_v, sum_salt, sum_temp, sum_rho
	integer :: num_vertical_layer	
	! End of local variables ==================================================!
	! Update old variables to the new variables for the next time step ========!
 	!$omp parallel
 	!$omp workshare
	! at cell
	eta_cell = eta_cell_new
	wn_cell = wn_cell_new
	dz_cell = dz_cell_new
	dzhalf_cell = dzhalf_cell_new
	salt_cell = salt_cell_new
	temp_cell = temp_cell_new
	
	! at face
	un_face = un_face_new
	vn_face = vn_face_new
	dz_face = dz_face_new
	dzhalf_face = dzhalf_face_new
	
	! at node
	dz_node = dz_node_new
	dzhalf_node = dzhalf_node_new
 	!$omp end workshare nowait
	
	! calculate vertically averaged values at node ============================!
	! for u, v, salt, temp, rho (water density)
	! This is for the times series & 2D output
	!$omp do private(i,k,num_vertical_layer,sum_u,sum_v,sum_salt,sum_temp,sum_rho)
   do i=1,maxnod
   	num_vertical_layer = top_layer_at_node(i) - bottom_layer_at_node(i) + 1 

   	! if it this node is dry, set nodal values to 0.0 and move to the next node
   	! if this routine is missed, there will be 'NaN' for dry node since num_vertical_layer is 0.
   	if(num_vertical_layer == 0) then ! dry node
   		ubar_node(i) = 0.0
   		vbar_node(i) = 0.0
   		sbar_node(i) = 0.0
   		tbar_node(i) = 0.0
   		rbar_node(i) = 0.0
   		cycle
   	end if
   	
   	! nodal velocities are defined at vertical level
   	sum_u = 0.0_dp
   	sum_v = 0.0_dp
   	do k=bottom_layer_at_node(i)-1,top_layer_at_node(i) ! since u_node is defined at vertical level
   		sum_u = sum_u + u_node(k,i)
   		sum_v = sum_v + v_node(k,i)
   	end do
   	ubar_node(i) = sum_u/(num_vertical_layer + 1) ! since u_node is defined at vertical level
   	vbar_node(i) = sum_v/(num_vertical_layer + 1) ! since v_node is defined at vertical level
   	
   	! nodal trasportative variables are defined at vertical layer
   	sum_salt = 0.0_dp
   	sum_temp = 0.0_dp
   	sum_rho = 0.0_dp
   	do k=bottom_layer_at_node(i),top_layer_at_node(i)
   		sum_salt = sum_salt + salt_node(k,i)
   		sum_temp = sum_temp + temp_node(k,i)
   		sum_rho = sum_rho + rho_node(k,i)
   	end do
   	sbar_node(i) = sum_salt/num_vertical_layer
   	tbar_node(i) = sum_temp/num_vertical_layer
   	rbar_node(i) = sum_rho/num_vertical_layer
   end do
	!$omp end do
	!$omp end parallel
	! end of calculation ======================================================!
end subroutine update_variables