!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!!
!! This is the ELM version of transport equation.
!! Backtracking starts from cell center.
!! ===========================================================================! 
subroutine solve_transport_equation_ELM_v7
   use mod_global_variables
   use mod_file_definition
   implicit none


   integer :: i, j, k, l
   integer :: kk, n1, n2, ie, nnel, jlev, num_vertical_layer, &
   &          klev, nd, icount
   integer :: kb, kt
   integer :: bt_step
   real(dp):: x0, y0, z0, uuint, vvint, wwint, vmag, bt_dt,   &
   &          xt, yt, zt, ttint, ssint
	real(dp):: specific_heat_pure_water
	! real(dp):: cp = 4186.0 ! specific_heat of water [J/KgC]

	! real(dp):: zup, zdown, dp1, dp2, sradp1, sradp2
	
	real(dp):: u_temp, v_temp
   real(dp),dimension(maxlayer+1,maxele) :: tc_bt, sc_bt ! backtracted temp & salt at cell
   real(dp),dimension(maxlayer+1,maxele) :: trhs, srhs
   real(dp),dimension(maxlayer+1)   		:: a_lower_mat, b_diagonal_mat, c_upper_mat, gam
   real(dp),dimension(maxlayer+1,5) 		:: soln, rrhs
   
   real(dp):: sum1, sum2, sum0
   real(dp):: rtemp1
   integer :: i_which_backtrack ! transport variable (location)
   ! End of local variables ==================================================!

	! i_which_backtrack:
	! 			1 = momentum transport; defined at face
	! 			2 = salt,temp, and other transport; defined at cell center
	i_which_backtrack = 2 ! since this is for monmentum transport

   specific_heat_pure_water = 4182.0 ! [J/kg.C], specific heat of pure water; originally it was 4184.0
   
! 	write(*,*) 'jww_1'
	
	! before start, let's prepare boundary conditions first ===================!
   call find_openboundary_salt_temp_v3
! 	write(*,*) num_ob_cell
! 	do i=1,num_ob_cell
! 		do k=1,maxlayer
! 			write(*,*) i,k, salt_at_obck(k,i), temp_at_obck(k,i)
! 		end do
! 	end do
! 	stop 'jw'

	! Backtracking from cell center at layer ==================================!
	! write(*,*) 'jw1'
	!$omp parallel
	!$omp do private(i,k,l, 										&
	!$omp &			  u_temp,v_temp,bt_step,bt_dt,nd,		&
	!$omp &			  nnel,jlev,x0,y0,z0,xt,yt,zt,			&
	!$omp &			  uuint,vvint,wwint,vmag,ttint,ssint)	
	do i=1,maxele
		! calculate backtracking only for wet elements
		if(top_layer_at_element(i) == 0) then
			! just use previous temp and salt at this element (i.e., do not backtrack)
			do k=1,maxlayer
				tc_bt(k,i) = 0.0_dp ! temp_cell(k,i)
				sc_bt(k,i) = 0.0_dp ! salt_cell(k,i)
			end do
			! Note: this "cycle" statement will hurt the code if "omp" is used.
			! that is why I used "if ~ else" instead "if ~ cycle"
			! cycle ! skip rest parts of the main do-loop (i.e., go to the next element)
		else
			! if this element is wet, perform backtracking.
	      do k=bottom_layer_at_element(i),top_layer_at_element(i)
				! initialize (x0,y0,z0),nnel and vel.
				! caution ! nnel must be initialized inside this loop as it is updated inside.
	         nnel = i
	         jlev = k
	         x0 = x_cell(i)
	         y0 = y_cell(i)
	         if(k == top_layer_at_element(nnel)) then
	         	z0 = MSL + eta_cell(nnel) - dz_cell(k,nnel)*0.5_dp
	         else
	            z0 = z_level(k) - dz_cell(k,nnel)*0.5_dp
	         end if
	
				! find velocities at the element center
	         u_temp = 0.0
	         v_temp = 0.0
	         bt_step = 0	
				do l=1,tri_or_quad(i)
					nd = nodenum_at_cell(l,i)
					u_temp = u_temp + (u_node(k,nd) + u_node(k-1,nd))*0.5_dp ! u at middle of this vertical layer at node 
					v_temp = v_temp + (v_node(k,nd) + v_node(k-1,nd))*0.5_dp ! v at middle of this vertical layer at node 
					bt_step = bt_step + num_sub_elm_iteration(k,nd) ! backtracking steps
				end do
				u_temp = u_temp/tri_or_quad(i) ! average nodal velocities
				v_temp = v_temp/tri_or_quad(i) ! average nodal velocities
				bt_step = int(bt_step/tri_or_quad(i)) ! average backtracking steps
				
				uuint = u_temp ! u at middle of this vertical layer at cell 
	         vvint = v_temp ! v at middle of this vertical layer at cell 
	         wwint = (wn_cell(k,i) + wn_cell(k-1,i))*0.5_dp ! w at middle of this vertical layer at cell
	         vmag = sqrt(uuint**2 + vvint**2 + wwint**2) ! true velocity magnitude at cell center
	
	         if(vmag <= 1.e-4) then
					! no backtracking.
	            tc_bt(k,i) = temp_cell(k,i)
	            sc_bt(k,i) = salt_cell(k,i)
	         else
	         	! Find backtracted salinity and temperatures for this element and vertical layer
	            bt_dt = dt/bt_step ! sub_time [sec] in backtracking
	
	            call ELM_backtrace_v1(i_which_backtrack,i,bt_step,bt_dt,uuint,vvint,wwint, &
	            &							x0,y0,z0,xt,yt,zt,nnel,jlev,ttint,ssint)
	
	            tc_bt(k,i) = ttint
	            sc_bt(k,i) = ssint
	         end if
	      end do ! k=bottom_layer_at_element(i),top_layer_at_element(i)
	   end if ! if(top_layer_at_element(i) == 0) then
   end do ! i=1,maxele
 	!$omp end do ! no wait
	! stop 'jw2'
	
	! transport equation starts ===============================================!
	! (0) at elements ---------------------------------------------------------!
	! prepare right-hand-side vector: (ELM + horizontal diffusion)
	! actually, I can put this loop into previous backtracking do-loop, 
	! but, I am do this work separately to make it easy to understand.
	!$omp do private(i,j,k,l,ie, sum1,sum2,rtemp1)
   do i=1,maxele
      do k=bottom_layer_at_element(i),top_layer_at_element(i) ! wet element
			! calculate horizontal diffusion term
			sum1 = 0.0_dp ! for salt
			sum2 = 0.0_dp ! for temp
			do l=1,tri_or_quad(i)
				j = facenum_at_cell(l,i)
				ie = adj_cellnum_at_face(2,j)
				if(ie /= 0) then ! if adjacent cell exists, include horizontal diffusion. Otherwise, do not include horizontal diffusion.
					rtemp1 = face_length(j)*dz_face(k,j)*Kh(k,j)/delta_j(j)
					sum1 = sum1 + rtemp1*(salt_cell(k,ie) - salt_cell(k,i))
					sum2 = sum2 + rtemp1*(temp_cell(k,ie) - temp_cell(k,i))
				end if
			end do
			sum1 = sum1*dt/(area(i)*dz_cell(k,i))
			sum2 = sum2*dt/(area(i)*dz_cell(k,i))
			
			! add horizontal advection & diffusion
         trhs(k,i) = tc_bt(k,i) + sum1 ! ELM + horizontal diffusion
         srhs(k,i) = sc_bt(k,i) + sum2 ! ELM + horizontal diffusion

			
			! jw, let's skip this part now...
			! heat conservation
!            if(i_heat_model_flag  /= 0) then
!               if(k == top_layer_at_node(i)) then
!                  zup = MSL + eta_at_node(i)
!               else
!                  zup = z_level(k)
!               end if
!               
!               zdown = zup - dz_node(k,i)
!               dp1 = min(MSL+eta_at_node(i)-zup,1.d2) ! to prevent underflow
!               dp2 = min(MSL+eta_at_node(i)-zdown,1.d2) ! to prevent underflow
!               
!               if(dp1 < 0.0 .or. dp2 < 0.0) then
!               	write(pw_run_log,*) 'depth<0 in nodal transport',i,dp1,dp2
!               end if
!               
!               sradp1 = srad(i)*(0.8*dexp(-dp1/0.9)+0.2*dexp(-dp1/2.1))/ref_water_density/specific_heat_pure_water ! jw, solar radiation? ummm, noway...
!               if(k == bottom_layer_at_node(i)) then
!                  sradp2 = 0
!               else
!                  sradp2 = srad(i)*(0.8*dexp(-dp2/0.9)+0.2*dexp(-dp2/2.1))/ref_water_density/specific_heat_pure_water
!               end if
!               
!               ! trhs(i,j) = trhs(i,j)+trans_ext_iter*dt*(sradp1-sradp2)
!               trhs(i,j) = trhs(i,j) + dt*(sradp1-sradp2)
!               if(k == top_layer_at_node(i)) then
!                  ! trhs(i,j) = trhs(i,j)+trans_ext_iter*dt*sflux(i)/ref_water_density/specific_heat_pure_water
!                  trhs(i,j) = trhs(i,j) + dt*sflux(i)/ref_water_density/specific_heat_pure_water
!               end if
!            end if ! i_heat_model_flag/=0
      end do ! k=bottom_layer_at_node(i),top_layer_at_node(i)
   end do ! i=1,maxele
	!$omp end do
! 	write(*,*) 'jww_4'


	! solve the transport eqs. and find new temp. and salt --------------------!
	! Note, this part is exactly same as the matrix preparation part in "solve_momentum_equation.f90"
	!$omp do private(i,k,kk,num_vertical_layer, &
	!$omp &			  a_lower_mat,b_diagonal_mat,c_upper_mat,rrhs,soln,gam,klev)
   do i=1,maxele
      if(top_layer_at_element(i)==0) then
      	cycle
      end if

      num_vertical_layer = top_layer_at_element(i)-bottom_layer_at_element(i)+1

		! Prepare coefficient matrices, A --------------------------------------!
		! From here, I will use "kk" for matrix ordering, and "k" for corresponding variable's layer ordering
		! see more details in "solve_momentum_equation.f90"
		! at mid layers
      do k=bottom_layer_at_element(i)+1,top_layer_at_element(i)-1
         kk = top_layer_at_element(i) - k + 1 ! row #
         a_lower_mat(kk) = -dt * Kv(k  ,i)/dzhalf_cell(k  ,i)
         c_upper_mat(kk) = -dt * Kv(k-1,i)/dzhalf_cell(k-1,i)
         b_diagonal_mat(kk) = -a_lower_mat(kk) + dz_cell(k,i) - c_upper_mat(kk)
      end do

		! at top & bottom layers
      if(top_layer_at_element(i)==bottom_layer_at_element(i)) then ! if it is 2D
      	! if it is 2D, we only need the diagonal matrix term, b_diagonal_mat
         b_diagonal_mat(1) = dz_cell(top_layer_at_element(i),i)
      else 
			! at the top(surface) layer
			! Only b and c terms are required
         c_upper_mat(1) = -dt + Kv(top_layer_at_element(i)-1,i) / dzhalf_cell(top_layer_at_element(i)-1,i)
         b_diagonal_mat(1) = dz_cell(top_layer_at_element(i),i) - c_upper_mat(1)

			! at the bottom layer
			! Only a and b terms are requried
         a_lower_mat(num_vertical_layer) = -dt * Kv(bottom_layer_at_element(i),i) / dzhalf_cell(bottom_layer_at_element(i),i)	         
         b_diagonal_mat(num_vertical_layer) = - a_lower_mat(num_vertical_layer) + dz_cell(bottom_layer_at_element(i),i)            
     	end if
		! End of Matrix A ------------------------------------------------------!
		
		! right-hand-side vector (rhs) -----------------------------------------!
     	do k=1,num_vertical_layer
     		kk = top_layer_at_element(i) - k + 1 ! row #
        	! rrhs(k,1)=trhs(kk,i) ! for temperature
        	! rrhs(k,2)=srhs(kk,i) ! for salinity
        	
        	rrhs(k,1)=dz_cell(kk,i)*trhs(kk,i) ! for temperature
        	rrhs(k,2)=dz_cell(kk,i)*srhs(kk,i) ! for salinity
     	end do
		
		! Now, we are ready to solve the maxtrix using Thomas Algorithm. 
		! soln for temp. and salt at new level
     	call tridiagonal_solver(maxlayer+1,num_vertical_layer,2,a_lower_mat,b_diagonal_mat,c_upper_mat,rrhs,soln,gam)

     	do k=1,num_vertical_layer
        	klev = top_layer_at_element(i)+1-k
        	temp_cell_new(klev,i) = soln(k,1)
        	salt_cell_new(klev,i) = soln(k,2)
			
			! write(*,*), i,k,soln(k,1),soln(k,2)
			
			
			! check against validity range
			! jw, I will activate this later...
!           	if(salt_node(klev,i) < salinity_min .or. salt_node(klev,i) > salinity_max) then
!              	write(pw_run_log,*) 'reset nodal salinity:',it,i,klev,salt_node(klev,i)
!              	salt_node(klev,i) = max(salinity_min,min(salt_node(klev,i),salinity_max))
!           	end if
!           	if(temp_node(klev,i) < temperature_min .or. temp_node(klev,i) > temperature_max) then
!            	write(pw_run_log,*) 'reset nodal temp:',it,i,klev,temp_node(klev,i)
!              	temp_node(klev,i) = max(temperature_min,min(temp_node(klev,i),temperature_max))
!           end if
     	end do ! k=1,num_vertical_layer
	end do ! i=1,maxele
	!$omp end do
! 		write(*,*) 'jww_5'
! 		write(*,*) 'jww_6'
					
	! extend temp & salt for below bottom & above top layers for wet nodes
	!$omp do private(i,k,kb,kt)
   do i=1,maxele
      if(top_layer_at_element(i)==0) then
      	cycle
      end if
      
      ! if this is a wet element...
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

! 	 write(*,*) 'jww_7'
	! stop 'jw'
	
	! update nodal values using cell values ===================================!
	!$omp do private(i,k,l,ie,sum0,sum1,sum2,icount)
	do i=1,maxnod
		! note, here I am doing from top to bottom, and there is reason
		! previously, there was an unintended node value (i.e., 0.0) could be asigned when all adjacent cells are located at the next vertical layer than this node.
		! To prevent this error, I updated this part as follows:
		! If all adjacent cells are located above this node's vertical level, 
		! just use the top node value for the current nodal value.
		do k=top_layer_at_node(i),bottom_layer_at_node(i),-1
			sum0 = 0.0_dp ! total volume
			sum1 = 0.0_dp ! total mass, temp
			sum2 = 0.0_dp ! total mass, salt
			
			icount = 0 ! count adjacent cells which are in the same vertical layer with this node
			do l = 1,adj_cells_at_node(i)
				ie = adj_cellnum_at_node(l,i)
				! if(dz_cell_new(k,ii) > 0.0_dp) then
				if(k>=bottom_layer_at_element(ie) .and. k<=top_layer_at_element(ie)) then					
					sum0 = sum0 + area(ie)*dz_cell_new(k,ie) ! total volume
					sum1 = sum1 + area(ie)*dz_cell_new(k,ie) * temp_cell_new(k,ie) ! total mass
					sum2 = sum2 + area(ie)*dz_cell_new(k,ie) * salt_cell_new(k,ie) ! total mass
				else
					icount = icount + 1
				end if
			end do

			if(sum0 == 0.0_dp) then
				temp_node(k,i) = 0.0_dp
				salt_node(k,i) = 0.0_dp
			else
				temp_node(k,i) = sum1/sum0 ! [mass]/[m^3], volumetric weighting average
				salt_node(k,i) = sum2/sum0 ! [mass]/[m^3], volumetric weighting average
			end if
			
			! if all adjacent cells are located above this node's vertical level,
			! just use the top node value for the current nodal value.
			if(icount == adj_cells_at_node(i)) then
				if(k /= top_layer_at_node(i)) then
					temp_node(k,i) = temp_node(k+1,i)
					salt_node(k,i) = salt_node(k+1,i)
				end if
			end if
		end do		
	end do
	!$omp end do



	
	! extend temp & salt for below bottom & above top layers for wet nodes
	!$omp do private(i,k,kb,kt)
   do i=1,maxnod
      if(top_layer_at_node(i)==0) then
      	cycle
      end if
      
      ! if this is a wet node...
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
	
	! update face values using nodal values ===================================!
	!$omp do private(j,k,n1,n2)
	do j=1,maxface
		n1 = nodenum_at_face(1,j)
		n2 = nodenum_at_face(2,j)
		! do k=bottom_layer_at_face(j),top_layer_at_face(j)
		do k=1,maxlayer
			salt_face(k,j) = (salt_node(k,n1) + salt_node(k,n2))*0.5
			temp_face(k,j) = (temp_node(k,n1) + temp_node(k,n2))*0.5
		end do
	end do
	!$omp end do
	!$omp end parallel
! 	 write(*,*) 'jww_8'

	! end of solving transport equation =======================================!
end subroutine solve_transport_equation_ELM_v7
  