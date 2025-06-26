!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Solve the nonlinear advection term using ELM
!! du/dt = -(udu/dx + vdu/dy + wdu/dz_face)
!! 
!! This will calculate:
!! 	un_ELM - the normal component of the velocity, caused by nonlinear advection, 
!! 				at each side of each vertical layer (at n+1 time step, which is the current time step)
!! 	vn_ELM - the tangential component of the velocity, caused by nonlinear advection, 
!!  				at each side of each vertical layer (at n+1 time step, which is the current time step)
!! 
subroutine solve_nonlinear_advection 
   use mod_global_variables
   implicit none

   integer :: i, j, k, l
   integer :: ie, id, nd, n1, n2, ie0, nnel, jlev, bt_step, iw, idelta
   real(dp):: summ, x0, y0, z0, devm, rl, uvel1, uvel2, vvel1, vvel2, dev, velo_w_node_temp, &
   &          uuint, vvint, wdown, wup, wwint, vmag, bt_dt, xt, yt, zt, ttint, ssint

   real(dp):: adv_turn_off_depth
   integer :: i_which_backtrack ! jw
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
	i_which_backtrack = 1 ! jw
	
	!$omp parallel
	
	! jw
	!$omp workshare
   un_ELM = 0.0_dp
   vn_ELM = 0.0_dp
   !$omp end workshare
   
	! jw
   if(elm_backtrace_flag == 1) then
   	!$omp do private(i,j,k,l, summ,ie,id,devm,nd,rl,uvel1,uvel2,vvel1,vvel2,dev,velo_w_node_temp,iw,idelta)
      do i = 1, maxnod
         do k = 1, maxlayer
            summ = 0.0 ! jw
            do j = 1, adj_cells_at_node(i) ! jw
               ie   = adj_cellnum_at_node(j,i) ! jw
               id   = node_count_each_element(j,i) ! jw
               devm = 0.0
               do l = 1, 2
                  if(l == 1) then
                     nd = nodenum_at_cell(start_end_node(tri_or_quad(ie), id, 1),ie)
                  else
                     nd = nodenum_at_cell(start_end_node(tri_or_quad(ie), id, tri_or_quad(ie)-1),ie)
                  end if
                  rl = SQRT((x_node(i)-x_node(nd))**2 + (y_node(i)-y_node(nd))**2)
                  uvel1 = (u_node(k, i) + u_node(k-1, i))*0.5_dp 	! jw
                  uvel2 = (u_node(k,nd) + u_node(k-1,nd))*0.5_dp	! jw
                  vvel1 = (v_node(k, i) + v_node(k-1, i))*0.5_dp	! jw
                  vvel2 = (v_node(k,nd) + v_node(k-1,nd))*0.5_dp	! jw
                  dev   = dsqrt((uvel1-uvel2)**2+(vvel1-vvel2)**2)/rl*dt ! jw
                  if(dev > devm) then
                  	devm = dev ! jw
                  end if
               end do   ! jw
               
               velo_w_node_temp = (w_node(k,i) + w_node(k-1,i))*0.5_dp ! jw
               iw = INT(2*velo_w_node_temp*dt/delta_z_min) ! jw
               idelta = MAX(int(devm/1.e-1),iw) ! jw
               summ = summ + DBLE(idelta)/adj_cells_at_node(i)
            end do   ! jw
            num_sub_elm_iteration(k,i) = MAX(ELM_min_iter,MIN(ELM_max_iter,INT(summ)))
         end do   ! jw
      end do   ! jw
   	!$omp end do
   end if   ! jw
	
	! jw
	! jw
	
 	!$omp do private(j,k, n1,n2,adv_turn_off_depth,ie0, &
 	!$omp &	nnel,jlev,x0,y0,z0,xt,yt,zt,uuint,vvint,wdown,wup,wwint,vmag,bt_step,bt_dt,ttint,ssint) 
   do j = 1, maxface
      n1 = nodenum_at_face(1,j)
      n2 = nodenum_at_face(2,j)
      
      adv_turn_off_depth = MAX(h_node(n1),h_node(n2))
      
      ! jw
      ! jw
      if(ABS(adv_turn_off_depth) < adv_onoff_depth) then ! jw
	   	do k=bottom_layer_at_face(j),top_layer_at_face(j)
	        	! jw
	         un_ELM(k,j) = un_face(k,j)	
	         vn_ELM(k,j) = vn_face(k,j)
	      end do
         cycle ! jw
      end if
            
      ! jw
      if(top_layer_at_element(adj_cellnum_at_face(1,j)) /= 0) then
      	! jw
         ie0 = adj_cellnum_at_face(1,j)
      else if(adj_cellnum_at_face(2,j) /= 0 .and. top_layer_at_element(adj_cellnum_at_face(2,j)) /= 0) then
      	! jw
      	! jw
      	! jw
         ie0 = adj_cellnum_at_face(2,j)
      else ! jw
      	! jw
      	! jw
      	! jw
      	
      	! jw
      	! jw
         do k = bottom_layer_at_face(j), top_layer_at_face(j)
            un_ELM(k,j) = un_face(k,j)
            vn_ELM(k,j) = vn_face(k,j)            
         end do
         cycle ! jw
      end if
		
		! jw
		! jw
		! jw
      do k = bottom_layer_at_face(j), top_layer_at_face(j) ! jw
			! jw
			! jw
         nnel = ie0
         jlev = k
         
         ! jw
         ! jw
         x0   = x_face(j)
         y0   = y_face(j)
         if(k == top_layer_at_element(nnel)) then
            z0 = MSL+eta_cell(nnel)-dz_cell(k,nnel)*0.5_dp
         else
            z0 = z_level(k)-dz_cell(k,nnel)*0.5_dp
         end if
			
			! jw
         uuint = un_face(k,j)*cos_theta(j)   &	! jw
         &  	- vn_face(k,j)*sin_theta(j)
         vvint = un_face(k,j)*sin_theta(j)   &	! jw
         &   	+ vn_face(k,j)*cos_theta(j)         
         wdown = (w_node(k-1,n1) + w_node(k-1,n2))*0.5_dp ! jw
         wup   = (w_node(k  ,n1) + w_node(k  ,n2))*0.5_dp ! jw
         wwint = (wdown + wup)*0.5_dp ! jw
         vmag  =  SQRT(uuint**2 + vvint**2 + wwint**2) ! jw
         
         ! jw
         ! jw
         ! jw
         if(vmag <= 1.e-4) then 
				! jw
				! jw
            un_ELM(k,j) = un_face(k,j)
            vn_ELM(k,j) = vn_face(k,j)				
         else
        		! jw
            bt_step = INT((num_sub_elm_iteration(k,n1) + num_sub_elm_iteration(k,n2))*0.5_dp) ! jw
            bt_dt   =  dt/bt_step ! jw
            
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
            call ELM_backtrace_v1(i_which_backtrack,j,bt_step,bt_dt,uuint,vvint,wwint, x0,y0,z0,xt,yt,zt,nnel,jlev,ttint,ssint)
				
				! jw
            un_ELM(k,j) =  uuint*cos_theta(j) + vvint*sin_theta(j)	! jw
            vn_ELM(k,j) = -uuint*sin_theta(j) + vvint*cos_theta(j)	! jw
         end if
		end do   ! jw
	end do ! jw
 	!$omp end do
 	!$omp end parallel
end subroutine solve_nonlinear_advection 
