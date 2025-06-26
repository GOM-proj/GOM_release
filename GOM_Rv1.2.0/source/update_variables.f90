!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine update_variables
	use mod_global_variables
	implicit none
	
	integer :: i, k
	real(dp):: sum_u, sum_v, sum_salt, sum_temp, sum_rho
	integer :: num_vertical_layer	
	! jw
	! jw
 	!$omp parallel
 	!$omp workshare
	! jw
	eta_cell = eta_cell_new
	wn_cell = wn_cell_new
	dz_cell = dz_cell_new
	dzhalf_cell = dzhalf_cell_new
	salt_cell = salt_cell_new
	temp_cell = temp_cell_new
	
	! jw
	un_face = un_face_new
	vn_face = vn_face_new
	dz_face = dz_face_new
	dzhalf_face = dzhalf_face_new
	
	! jw
	dz_node = dz_node_new
	dzhalf_node = dzhalf_node_new
 	!$omp end workshare nowait
	
	! jw
	! jw
	! jw
	!$omp do private(i,k,num_vertical_layer,sum_u,sum_v,sum_salt,sum_temp,sum_rho)
   do i=1,maxnod
   	num_vertical_layer = top_layer_at_node(i) - bottom_layer_at_node(i) + 1 

   	! jw
   	! jw
   	if(num_vertical_layer == 0) then ! jw
   		ubar_node(i) = 0.0
   		vbar_node(i) = 0.0
   		sbar_node(i) = 0.0
   		tbar_node(i) = 0.0
   		rbar_node(i) = 0.0
   		cycle
   	end if
   	
   	! jw
   	sum_u = 0.0_dp
   	sum_v = 0.0_dp
   	do k=bottom_layer_at_node(i)-1,top_layer_at_node(i) ! jw
   		sum_u = sum_u + u_node(k,i)
   		sum_v = sum_v + v_node(k,i)
   	end do
   	ubar_node(i) = sum_u/(num_vertical_layer + 1) ! jw
   	vbar_node(i) = sum_v/(num_vertical_layer + 1) ! jw
   	
   	! jw
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
	! jw
end subroutine update_variables