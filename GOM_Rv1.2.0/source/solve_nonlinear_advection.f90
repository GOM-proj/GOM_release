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
   integer :: i_which_backtrack ! transport material definition
   ! end of local variables ==================================================!
	
	! solving nonlinear term(elm) by backtracking =============================!
	! advection flag and elm backtracing method flag --------------------------!
	! ELM_backtrace_flag = 0 :
	! 		a constant value, specified below, is used throughout the domain.
	! 		the next line is: num_sub_elm_iteration = constant number of subdivisions.
	! ELM_backtrace_flag = 1 :
	! 		the number of subdivisions is automatically calculated in the code based on the local
	!     velocity gradient, subject to the max. and min. 
	! 		the next line is: ELM_min_iter, ELM_max_iter : minimum and maximum number of subdivisions allowable
	! -------------------------------------------------------------------------!
	
	! i_which_backtrack:
	! 			1 = momentum transport; defined at face
	! 			2 = salt,temp, and other transport; defined at cell center
	i_which_backtrack = 1 ! since this is for monmentum transport
	
	!$omp parallel
	
	! allocate nonlinear advection velocy by ELM
	!$omp workshare
   un_ELM = 0.0_dp
   vn_ELM = 0.0_dp
   !$omp end workshare
   
	! compute number of subdivisions, num_sub_elm_iteration(k,i) when elm_backtrace_flag=1
   if(elm_backtrace_flag == 1) then
   	!$omp do private(i,j,k,l, summ,ie,id,devm,nd,rl,uvel1,uvel2,vvel1,vvel2,dev,velo_w_node_temp,iw,idelta)
      do i = 1, maxnod
         do k = 1, maxlayer
            summ = 0.0 ! I use 'summ' to avoid confusion with SUM() function.
            do j = 1, adj_cells_at_node(i) ! total number of adjacent elements surrounding this node
               ie   = adj_cellnum_at_node(j,i) ! surrounding cell ID
               id   = node_count_each_element(j,i) ! the location of this node appears in elements which include this node
               devm = 0.0
               do l = 1, 2
                  if(l == 1) then
                     nd = nodenum_at_cell(start_end_node(tri_or_quad(ie), id, 1),ie)
                  else
                     nd = nodenum_at_cell(start_end_node(tri_or_quad(ie), id, tri_or_quad(ie)-1),ie)
                  end if
                  rl = SQRT((x_node(i)-x_node(nd))**2 + (y_node(i)-y_node(nd))**2)
                  uvel1 = (u_node(k, i) + u_node(k-1, i))*0.5_dp 	! u_node at layer
                  uvel2 = (u_node(k,nd) + u_node(k-1,nd))*0.5_dp	! u_node at layer
                  vvel1 = (v_node(k, i) + v_node(k-1, i))*0.5_dp	! v_node at layer
                  vvel2 = (v_node(k,nd) + v_node(k-1,nd))*0.5_dp	! v_node at layer
                  dev   = dsqrt((uvel1-uvel2)**2+(vvel1-vvel2)**2)/rl*dt ! (dv/dx)*dt
                  if(dev > devm) then
                  	devm = dev ! update the velocity gradient
                  end if
               end do   ! do l = 1, 2
               
               velo_w_node_temp = (w_node(k,i) + w_node(k-1,i))*0.5_dp ! w_node at layer
               iw = INT(2*velo_w_node_temp*dt/delta_z_min) ! jw, only interger part
               idelta = MAX(int(devm/1.e-1),iw) ! from now I use MAX() not MAX0, which is an archaic (old form) function of max function
               summ = summ + DBLE(idelta)/adj_cells_at_node(i)
            end do   ! do j = 1, adj_cells_at_node(i)
            num_sub_elm_iteration(k,i) = MAX(ELM_min_iter,MIN(ELM_max_iter,INT(summ)))
         end do   ! do k = 1, maxlayer
      end do   ! do i=1,maxnod
   	!$omp end do
   end if   ! if(elm_backtrace_flag == 1)
	
	! Calculate un_ELM & vn_ELM at face center of each vertical layer =========!
	! from centers of faces
	
 	!$omp do private(j,k, n1,n2,adv_turn_off_depth,ie0, &
 	!$omp &	nnel,jlev,x0,y0,z0,xt,yt,zt,uuint,vvint,wdown,wup,wwint,vmag,bt_step,bt_dt,ttint,ssint) 
   do j = 1, maxface
      n1 = nodenum_at_face(1,j)
      n2 = nodenum_at_face(2,j)
      
      adv_turn_off_depth = MAX(h_node(n1),h_node(n2))
      
      ! if this side depth is less than given advection turn off depth criteria, adv_onoff_depth,
      ! do not calculate nonlinear advection, but just include explicit velocities.
      if(ABS(adv_turn_off_depth) < adv_onoff_depth) then ! advt=0
	   	do k=bottom_layer_at_face(j),top_layer_at_face(j)
	        	! do not update the velocity field from nonlinear advection process but keep the previous value.
	         un_ELM(k,j) = un_face(k,j)	
	         vn_ELM(k,j) = vn_face(k,j)
	      end do
         cycle ! go to the next face
      end if
            
      ! find a wet cell, and if there is no wet cell, there will be no nonlinear advection.
      if(top_layer_at_element(adj_cellnum_at_face(1,j)) /= 0) then
      	! if first cell is not dry (i.e., wet cell)
         ie0 = adj_cellnum_at_face(1,j)
      else if(adj_cellnum_at_face(2,j) /= 0 .and. top_layer_at_element(adj_cellnum_at_face(2,j)) /= 0) then
      	! if first cell is dry &
      	!    second cell exists & 
      	!    the second cell is not dry (i.e., wet)
         ie0 = adj_cellnum_at_face(2,j)
      else ! do not calculate nonlinear advection
      	! if first cell is dry &
      	!    second cell doesn't exist or
      	!    second cell exist but dry
      	
      	! i.e. both cells are dry, same as above
      	! do not update the velocity field from nonlinear advection process but keep the previous value.
         do k = bottom_layer_at_face(j), top_layer_at_face(j)
            un_ELM(k,j) = un_face(k,j)
            vn_ELM(k,j) = vn_face(k,j)            
         end do
         cycle ! skip rest parts of the main do-loop (i.e., go to the next face)
      end if
		
		! Finally, if one of cells is wet, calculate nonlinear advectioin term
		! The backtracking starts from the center of each face, where velocity is specified,
		! i.e., (x0,y0,z0), with true velocities components, i.e., (uuint,vvint,wwint)
      do k = bottom_layer_at_face(j), top_layer_at_face(j) ! wet side
			! initialize (x0,y0,z0),nnel and vel.
			! caution! nnel must be initialized inside this loop as it is updated inside.
         nnel = ie0
         jlev = k
         
         ! set horizontal position at the face center,
         ! and set vertical position, z0, in the middle of the vertical layer
         x0   = x_face(j)
         y0   = y_face(j)
         if(k == top_layer_at_element(nnel)) then
            z0 = MSL+eta_cell(nnel)-dz_cell(k,nnel)*0.5_dp
         else
            z0 = z_level(k)-dz_cell(k,nnel)*0.5_dp
         end if
			
			! scopes for w_node etc. have been extended
         uuint = un_face(k,j)*cos_theta(j)   &	! true u (east/west velocity)
         &  	- vn_face(k,j)*sin_theta(j)
         vvint = un_face(k,j)*sin_theta(j)   &	! true v (nort/south velocity)
         &   	+ vn_face(k,j)*cos_theta(j)         
         wdown = (w_node(k-1,n1) + w_node(k-1,n2))*0.5_dp ! w at top face of this layer
         wup   = (w_node(k  ,n1) + w_node(k  ,n2))*0.5_dp ! w at bottom face of this layer
         wwint = (wdown + wup)*0.5_dp ! vertical velocity magnitude at cell center of the water prism
         vmag  =  SQRT(uuint**2 + vvint**2 + wwint**2) ! true velocity magnitude
         
         ! This is an additional criteria to turn-off nonlinear advection.
         ! If true velocity maginitude is very small, turn off nonlinear advection.
         ! Actually, this criteria should be in main.inp, but I will keep this only in the source code.
         if(vmag <= 1.e-4) then 
				! do not update the velocity field from nonlinear advection process 
				! but keep the previous value, i.e., tendency term only.
            un_ELM(k,j) = un_face(k,j)
            vn_ELM(k,j) = vn_face(k,j)				
         else
        		! Finally, calculate nonlinear advection term 	
            bt_step = INT((num_sub_elm_iteration(k,n1) + num_sub_elm_iteration(k,n2))*0.5_dp) ! number of backtracking steps
            bt_dt   =  dt/bt_step ! sub_time_step [sec], in backtracking
            
            ! The first number should be 1 or 2:
            ! 		1: barotropic backtracking
            !  	2: baroclinic backtracking
            ! See more details in ELM_backtrace.f90
            
            ! Note: if ELM_backtrace_flag is chosen other than 0 and 1, num_sub_elm_iteration will have 0.
            ! Thus, bt_step will be 0, and bt_dt will become infinity.
            ! And, it makes ELM unconditionally stable; anyway it will find a solution
            ! write(*,*) 'ELM_bactrace_flag = ', ELM_backtrace_flag, 'bt_step=', bt_step, 'bt_dt=', bt_dt
            ! where:
            ! 		j				: j'th side (i.e., this side)
            ! 		bt_step		: number of backtracking steps
            ! 		bt_dt			: sub time step [sec] for backtracking
            ! 		uuint,vvint,wwint,ttint,ssint	: true velocities, temperature, and salinity at this side of layer
            ! 		(x0,y0,k0)	: starting position in jth side of kth layer
            ! 		(xt,yt,zt)	: final backtracked position
            ! 		nnel			: adj_cellnum_at_face(1,j) or adj_cellnum_at_face(2,j)
            ! 		jlev			: vertical layer at this side "j"
            call ELM_backtrace_v1(i_which_backtrack,j,bt_step,bt_dt,uuint,vvint,wwint, x0,y0,z0,xt,yt,zt,nnel,jlev,ttint,ssint)
				
				! returning uuint & vvint includes face normal velocity: (un_ELM + un_face)
            un_ELM(k,j) =  uuint*cos_theta(j) + vvint*sin_theta(j)	! u* =  u x cos_theta + v x sin_theta
            vn_ELM(k,j) = -uuint*sin_theta(j) + vvint*cos_theta(j)	! v* = -u x sin_theta + v x cos_theta				
         end if
		end do   ! do k = bottom_layer_at_face(j),top_layer_at_face(j)
	end do ! do j = 1, maxface
 	!$omp end do
 	!$omp end parallel
end subroutine solve_nonlinear_advection 
