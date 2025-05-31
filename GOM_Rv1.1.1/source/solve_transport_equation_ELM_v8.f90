!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!!
!! This is the ELM version of transport equation.
!! Backtracking starts from cell center.
!! ===========================================================================! 
subroutine solve_transport_equation_ELM_v8
   use mod_global_variables
   use mod_file_definition
   implicit none


   integer :: i, j, k, l
   integer :: kk, n1, n2, ie, nnel, jlev, num_vertical_layer, &
   &          klev, nd, icount
   integer :: kb, kt
  	integer :: t_layer, b_layer

   integer :: bt_step
   real(dp):: x0, y0, z0, uuint, vvint, wwint, vmag, bt_dt,   &
   &          xt, yt, zt, ttint, ssint
	
	real(dp):: u_temp, v_temp
   real(dp),dimension(maxlayer+1,maxele) :: tc_bt, sc_bt ! jw
   real(dp),dimension(maxlayer+1,maxele) :: trhs, srhs
   real(dp),dimension(maxlayer+1)   		:: a_lower_mat, b_diagonal_mat, c_upper_mat, gam
   real(dp),dimension(maxlayer+1,5) 		:: soln, rrhs

	real(dp):: cp = 4186.0 ! jw
   real(dp):: sum1, sum2, sum0
   real(dp):: rtemp1
   integer :: i_which_backtrack ! jw
   ! jw

	! jw
	! jw
	! jw
	i_which_backtrack = 2 ! jw
   
! jw
	
	! jw
   call find_openboundary_salt_temp_v3
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
	!$omp do private(i,k,l, 										&
	!$omp &			  u_temp,v_temp,bt_step,bt_dt,nd,		&
	!$omp &			  nnel,jlev,x0,y0,z0,xt,yt,zt,			&
	!$omp &			  uuint,vvint,wwint,vmag,ttint,ssint)	
	do i=1,maxele
		! jw
		if(top_layer_at_element(i) == 0) then
			! jw
			do k=1,maxlayer
				tc_bt(k,i) = 0.0_dp ! jw
				sc_bt(k,i) = 0.0_dp ! jw
			end do
			! jw
			! jw
			! jw
		else
			! jw
	      do k=bottom_layer_at_element(i),top_layer_at_element(i)
				! jw
				! jw
	         nnel = i
	         jlev = k
	         x0 = x_cell(i)
	         y0 = y_cell(i)
	         if(k == top_layer_at_element(nnel)) then
	         	z0 = MSL + eta_cell(nnel) - dz_cell(k,nnel)*0.5_dp
	         else
	            z0 = z_level(k) - dz_cell(k,nnel)*0.5_dp
	         end if
	
				! jw
	         u_temp = 0.0
	         v_temp = 0.0
	         bt_step = 0	
				do l=1,tri_or_quad(i)
					nd = nodenum_at_cell(l,i)
					u_temp = u_temp + (u_node(k,nd) + u_node(k-1,nd))*0.5_dp ! jw
					v_temp = v_temp + (v_node(k,nd) + v_node(k-1,nd))*0.5_dp ! jw
					bt_step = bt_step + num_sub_elm_iteration(k,nd) ! jw
				end do
				u_temp = u_temp/tri_or_quad(i) ! jw
				v_temp = v_temp/tri_or_quad(i) ! jw
				bt_step = int(bt_step/tri_or_quad(i)) ! jw
				
				uuint = u_temp ! jw
	         vvint = v_temp ! jw
	         wwint = (wn_cell(k,i) + wn_cell(k-1,i))*0.5_dp ! jw
	         vmag = sqrt(uuint**2 + vvint**2 + wwint**2) ! jw
	
	         if(vmag <= 1.e-4) then
					! jw
	            tc_bt(k,i) = temp_cell(k,i)
	            sc_bt(k,i) = salt_cell(k,i)
	         else
	         	! jw
	            bt_dt = dt/bt_step ! jw
	
	            call ELM_backtrace_v1(i_which_backtrack,i,bt_step,bt_dt,uuint,vvint,wwint, &
	            &							x0,y0,z0,xt,yt,zt,nnel,jlev,ttint,ssint)
	
	            tc_bt(k,i) = ttint
	            sc_bt(k,i) = ssint
	         end if
	         
	         ! jw
	         tc_bt(k,i) = tc_bt(k,i) + dt*phi_sz(k,i)/(rho_o*cp*dz_cell(k,i))
	      end do ! jw
	   end if ! jw
   end do ! jw
 	!$omp end do ! jw
	! jw
	
	! jw
	! jw
	! jw
	! jw
	! jw
	!$omp do private(i,j,k,l,ie, sum1,sum2,rtemp1)
   do i=1,maxele
      do k=bottom_layer_at_element(i),top_layer_at_element(i) ! jw
			! jw
			sum1 = 0.0_dp ! jw
			sum2 = 0.0_dp ! jw
			do l=1,tri_or_quad(i)
				j = facenum_at_cell(l,i)
				ie = adj_cellnum_at_face(2,j)
				if(ie /= 0) then ! jw
					rtemp1 = face_length(j)*dz_face(k,j)*Kh(k,j)/delta_j(j)
					sum1 = sum1 + rtemp1*(salt_cell(k,ie) - salt_cell(k,i))
					sum2 = sum2 + rtemp1*(temp_cell(k,ie) - temp_cell(k,i))
				end if
			end do
			sum1 = sum1*dt/(area(i)*dz_cell(k,i))
			sum2 = sum2*dt/(area(i)*dz_cell(k,i))
			
			! jw
         trhs(k,i) = tc_bt(k,i) + sum1 ! jw
         srhs(k,i) = sc_bt(k,i) + sum2 ! jw
      end do ! jw
   end do ! jw
	!$omp end do
! jw


	! jw
	! jw
	!$omp do private(i,k,kk,num_vertical_layer, &
	!$omp &			  a_lower_mat,b_diagonal_mat,c_upper_mat,rrhs,soln,gam,klev)
   do i=1,maxele
      if(top_layer_at_element(i)==0) then
      	cycle
      end if

      num_vertical_layer = top_layer_at_element(i)-bottom_layer_at_element(i)+1

		! jw
		! jw
		! jw
		! jw
      do k=bottom_layer_at_element(i)+1,top_layer_at_element(i)-1
         kk = top_layer_at_element(i) - k + 1 ! jw
         a_lower_mat(kk) = -dt * Kv(k  ,i)/dzhalf_cell(k  ,i)
         c_upper_mat(kk) = -dt * Kv(k-1,i)/dzhalf_cell(k-1,i)
         b_diagonal_mat(kk) = -a_lower_mat(kk) + dz_cell(k,i) - c_upper_mat(kk)
      end do

		! jw
      if(top_layer_at_element(i)==bottom_layer_at_element(i)) then ! jw
      	! jw
         b_diagonal_mat(1) = dz_cell(top_layer_at_element(i),i)
      else 
			! jw
			! jw
         c_upper_mat(1) = -dt + Kv(top_layer_at_element(i)-1,i) / dzhalf_cell(top_layer_at_element(i)-1,i)
         b_diagonal_mat(1) = dz_cell(top_layer_at_element(i),i) - c_upper_mat(1)

			! jw
			! jw
         a_lower_mat(num_vertical_layer) = -dt * Kv(bottom_layer_at_element(i),i) / dzhalf_cell(bottom_layer_at_element(i),i)	         
         b_diagonal_mat(num_vertical_layer) = - a_lower_mat(num_vertical_layer) + dz_cell(bottom_layer_at_element(i),i)            
     	end if
		! jw
		
		! jw
     	do k=1,num_vertical_layer
     		kk = top_layer_at_element(i) - k + 1 ! jw
        	! jw
        	! jw
        	
        	rrhs(k,1)=dz_cell(kk,i)*trhs(kk,i) ! jw
        	rrhs(k,2)=dz_cell(kk,i)*srhs(kk,i) ! jw
     	end do
		
		! jw
		! jw
     	call tridiagonal_solver(maxlayer+1,num_vertical_layer,2,a_lower_mat,b_diagonal_mat,c_upper_mat,rrhs,soln,gam)

     	do k=1,num_vertical_layer
        	klev = top_layer_at_element(i)+1-k
        	temp_cell_new(klev,i) = soln(k,1)
        	salt_cell_new(klev,i) = soln(k,2)
			
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
     	end do ! jw
	end do ! jw
	!$omp end do
! jw
! jw
					
	! jw
	!$omp do private(i,k,kb,kt)
   do i=1,maxele
      if(top_layer_at_element(i)==0) then
      	cycle
      end if
      
      ! jw
      kb = bottom_layer_at_element(i)
      kt = top_layer_at_element(i)
      do k=1,kb-1
         temp_cell_new(k,i)=temp_cell_new(kb,i)
         salt_cell_new(k,i)=salt_cell_new(kb,i)
      end do
      do k=kt+1,maxlayer
         temp_cell_new(k,i)=temp_cell_new(kt,i)
         salt_cell_new(k,i)=salt_cell_new(kt,i)
      end do
   end do
	!$omp end do

! jw
	! jw
	
	! jw
	!$omp do private(i,k,l,ie,sum0,sum1,sum2,icount)
	do i=1,maxnod
		! jw
		! jw
		! jw
		! jw
		! jw
		do k=top_layer_at_node(i),bottom_layer_at_node(i),-1
			sum0 = 0.0_dp ! jw
			sum1 = 0.0_dp ! jw
			sum2 = 0.0_dp ! jw
			
			icount = 0 ! jw
			do l = 1,adj_cells_at_node(i)
				ie = adj_cellnum_at_node(l,i)
				! jw
				if(k>=bottom_layer_at_element(ie) .and. k<=top_layer_at_element(ie)) then					
					sum0 = sum0 + area(ie)*dz_cell_new(k,ie) ! jw
					sum1 = sum1 + area(ie)*dz_cell_new(k,ie) * temp_cell_new(k,ie) ! jw
					sum2 = sum2 + area(ie)*dz_cell_new(k,ie) * salt_cell_new(k,ie) ! jw
				else
					icount = icount + 1
				end if
			end do

			if(sum0 == 0.0_dp) then
				temp_node(k,i) = 0.0_dp
				salt_node(k,i) = 0.0_dp
			else
				temp_node(k,i) = sum1/sum0 ! jw
				salt_node(k,i) = sum2/sum0 ! jw
			end if
			
			! jw
			! jw
			if(icount == adj_cells_at_node(i)) then
				if(k /= top_layer_at_node(i)) then
					temp_node(k,i) = temp_node(k+1,i)
					salt_node(k,i) = salt_node(k+1,i)
				end if
			end if
		end do		
	end do
	!$omp end do



	
	! jw
	!$omp do private(i,k,kb,kt)
   do i=1,maxnod
      if(top_layer_at_node(i)==0) then
      	cycle
      end if
      
      ! jw
      kb = bottom_layer_at_node(i)
      kt = top_layer_at_node(i)
      do k=1,kb-1
         temp_node(k,i)=temp_node(kb,i)
         salt_node(k,i)=salt_node(kb,i)
      end do
      do k=kt+1,maxlayer
         temp_node(k,i)=temp_node(kt,i)
         salt_node(k,i)=salt_node(kt,i)
      end do
   end do
	!$omp end do
	
	! jw
	!$omp do private(j,k,n1,n2)
	do j=1,maxface
		n1 = nodenum_at_face(1,j)
		n2 = nodenum_at_face(2,j)
		! jw
		do k=1,maxlayer
			salt_face(k,j) = (salt_node(k,n1) + salt_node(k,n2))*0.5
			temp_face(k,j) = (temp_node(k,n1) + temp_node(k,n2))*0.5
		end do
	end do
	!$omp end do
	!$omp end parallel
! jw

	! jw
end subroutine solve_transport_equation_ELM_v8
  