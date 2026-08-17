!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! ===========================================================================!
!! routine for ELM backtracking.                                       
!! 	input: (x0,y0,z0) and nnel that encloses it, and initial vel.,           
!! 	       initial level (for search_straight_line), and a flag indicating   
!! 	       1st or 2nd tracking.                                              
!! 	output: destination pt (xt,yt,zt), element nnel and level jlev,          
!! 	and vel. there (uuint,vvint,wwint), and (s,t) there.                     
!!   	i_which_backtrack = 1 : barotrophic back tracking
!!    	               = 2 : baroclinic  back tracking
!! ===========================================================================!
!! call ELM_backtrace(1,bt_step,bt_dt,uuint,vvint,wwint, x0,y0,z0,xt,yt,zt,nnel,jlev,ttint,ssint,j)
subroutine ELM_backtrace_v1(i_which_backtrack,j_face,bt_step,bt_dt, & ! jw
&	                      uuint,vvint,wwint,x0,y0,z0,xt,yt,zt, &
&								 nnel,jlev,ttint,ssint) 
	use mod_global_variables
	use mod_file_definition
	use mod_function_library, only : calculate_area  
	implicit none
	
	integer, intent(in) 		:: i_which_backtrack, j_face, bt_step
	real(dp),intent(in) 		:: bt_dt
	integer, intent(inout) 	:: nnel,jlev
	real(dp),intent(inout) 	:: uuint,vvint,wwint,x0,y0,z0
	real(dp),intent(out) 	:: xt,yt,zt,ttint,ssint
	
	real(dp):: vxl(4,2),vyl(4,2),vzl(4,2) ! jw
	real(dp):: vxn(4),vyn(4),vzn(4)
	real(dp):: staint(4),t_xi(4),s_xi(4),sig(4),subrat(4)
	
	integer :: i, j, l, k1, k2, icount, index_1, index_2
	integer :: idt, iflqs1, nd, lev, node1, node2, node3, node4
	integer :: ibnd1, ibnd2
	
	real(dp):: trat, zup, zrat, zrat2, aa, aa1, aa2, aa3, aa4, csi, etta, zanchor, ta, sa, bb, bb1, bb2 
	! jw
	
	do idt = 1, bt_step ! jw
	   trat = dble(idt)/bt_step ! jw
	   
	   ! jw
	   ! jw
	   xt = x0-bt_dt*uuint
	   yt = y0-bt_dt*vvint
	   zt = z0-bt_dt*wwint
	   
	   ! jw
	   ! jw
	   ! jw
	   ! jw
	   ! jw
	   ! jw
	   call search_straight_line(i_which_backtrack,nnel,jlev,bt_dt,x0,y0,z0,xt,yt,zt,iflqs1,idt,j_face)
	
		! jw
		! jw
		do l = 1, tri_or_quad(nnel)
			nd = nodenum_at_cell(l,nnel)
			
			! jw
			! jw
			! jw
			! jw
			do k2 = 1, 2
				lev = jlev + k2 - 2 ! jw
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
				if(i_which_backtrack == 1 .or. i_which_backtrack == 2) then
					vxl(l,k2)=u_node(lev,nd) ! jw
					vyl(l,k2)=v_node(lev,nd) ! jw
					vzl(l,k2)=w_node(lev,nd) ! jw
				else
					! jw
				 	vxl(l,k2)=u_node(lev,nd)*(1.0d0-trat)+velo_u_transport(lev,nd)*trat
				 	vyl(l,k2)=v_node(lev,nd)*(1.0d0-trat)+velo_v_transport(lev,nd)*trat
				 	vzl(l,k2)=w_node(lev,nd)*(1.0d0-trat)+velo_w_transport(lev,nd)*trat
				end if	
			end do
		end do

		! jw
		! jw
		! jw
		if(jlev == top_layer_at_element(nnel)) then
			zup = MSL+eta_cell(nnel)
		else
			zup = z_level(jlev)
		end if
		
		if(zt > zup) then
			write(pw_run_log,*) 	'ELM_backtrace.f90: impossible vertical position: stop'
			write(pw_run_log,*)	'Final vertical position = ', zt, ', but the possible position = ', zup
			write(*,*) 				'ELM backtrace.f90: impossible vertical position: stop'
			write(*,*)				'Final vertical position = ', zt, ', but the possible position = ', zup
			stop
		end if
		
		! jw
		zrat=(zup-zt)/dz_cell(jlev,nnel)
		
		! jw
		do l = 1, tri_or_quad(nnel)
			vxn(l) = vxl(l,2)*(1-zrat)+vxl(l,1)*zrat
			vyn(l) = vyl(l,2)*(1-zrat)+vyl(l,1)*zrat
			vzn(l) = vzl(l,2)*(1-zrat)+vzl(l,1)*zrat
		end do
		
		! jw
		node1 = nodenum_at_cell(1,nnel)
		node2 = nodenum_at_cell(2,nnel)
		node3 = nodenum_at_cell(3,nnel)
		
		if(tri_or_quad(nnel)==3) then ! jw
			! jw
			staint(1) = calculate_area(xt, x_node(node2), x_node(node3), &
			&									yt, y_node(node2), y_node(node3))/area(nnel)
			staint(2) = calculate_area(x_node(node1), xt, x_node(node3), &
			&									y_node(node1), yt, y_node(node3))/area(nnel)
			staint(3) = calculate_area(x_node(node1), x_node(node2), xt, &
			&									y_node(node1), y_node(node2), yt)/area(nnel)
			
			! jw
			! jw
			staint(1) = MAX(0.0d0,MIN(1.0d0,staint(1)))
			staint(2) = MAX(0.0d0,MIN(1.0d0,staint(2)))
			if(staint(1)+staint(2) > 1) then
			   staint(3)=0
			   staint(2)=1.0_dp-staint(1)
			else
			   staint(3)=1.0_dp-staint(1)-staint(2)
			end if
			
		! jw
		else ! jw
			! jw
			node4 = nodenum_at_cell(4,nnel)
			aa1 = calculate_area(x_node(node1),x_node(node2),xt,   &
			&							y_node(node1),y_node(node2),yt)
			aa2 = calculate_area(x_node(node2),x_node(node3),xt,   &
			&							y_node(node2),y_node(node3),yt)
			aa3 = calculate_area(x_node(node3),x_node(node4),xt,   &
			&							y_node(node3),y_node(node4),yt)
			aa4 = calculate_area(x_node(node4),x_node(node1),xt,   &
			&							y_node(node4),y_node(node1),yt)
			aa = abs(aa1)+abs(aa2)+abs(aa3)+abs(aa4)
			
			if(abs(aa-area(nnel)) / area(nnel)>small_06) then
				write(pw_run_log,*)	'ELM_backtrace.f90: tracking is outside element: stop'
				write(pw_run_log,*)	'calculated area = ', aa, ', but element area = ', area(nnel)
				write(*,*)				'ELM_backtrace.f90: tracking is outside element: stop'
				write(*,*)				'calculated area = ', aa, ', but element area = ', area(nnel)
				stop
			end if
			
			call bilinear_interpolation(dble(nnel),   &
			&							x_node(node1), x_node(node2),   &
			&							x_node(node3), x_node(node4),   &
			&							y_node(node1), y_node(node2),   &
			&							y_node(node3), y_node(node4),   &
			&							xt, yt, csi, etta, staint)
      end if ! jw
		
		! jw
		! jw
		! jw
		uuint = 0.0d0
		vvint = 0.0d0
		wwint = 0.0d0		
		do j = 1, tri_or_quad(nnel)
			uuint = uuint+vxn(j)*staint(j)
			vvint = vvint+vyn(j)*staint(j)
			wwint = wwint+vzn(j)*staint(j)
		end do
		
		if(iflqs1 == 1) then
			exit
		end if
		
		! jw
		x0 = xt
		y0 = yt
		z0 = zt
	end do ! jw
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
	! errors may occur due to small errors!!!
	if(i_which_backtrack==2) then
		! jw
		! jw
		! jw
		ibnd1 = ob_element_flag(nnel) ! jw
		ibnd2 = Qb_element_flag(nnel) ! jw
		
		ttint = 0.0_dp
		ssint = 0.0_dp		
		if(ibnd1 > 0) then
			! jw
			ssint = salt_at_obck(jlev,ibnd1)
			ttint = temp_at_obck(jlev,ibnd1)
		else if(ibnd2 > 0) then
			! jw
			ssint = salt_at_Qbc(ibnd2)
			ttint = temp_at_Qbc(ibnd2)
		else
			! jw
			! jw
			! jw
			zanchor = zup - dz_cell(jlev,nnel)*0.5_dp ! jw
			if(jlev == top_layer_at_element(nnel) .and. zt>zanchor) then
				! jw
				k1    = jlev
				k2    = jlev
				zrat2 = 0.5_dp   ! jw
			else if(zt>zanchor) then ! jw
				! jw
				! jw
				k1    = jlev
				k2    = jlev+1
				zrat2 = (zt-zanchor)/dzhalf_cell(jlev,nnel)
			else ! jw
				! jw
				k1    = max(jlev-1,bottom_layer_at_element(nnel))
				k2    = jlev
				zrat2 = 1.0_dp - (zanchor-zt)/dzhalf_cell(k1,nnel)
			end if
			
			zrat2 = max(0.0_dp,min(1.0_dp,zrat2))
			
			icount = 0
			ta     = 0.0_dp
			sa     = 0.0_dp
			
			do l = 1, tri_or_quad(nnel)
				nd = nodenum_at_cell(l,nnel)
				if(top_layer_at_node(nd) /= 0) then
					icount  = icount+1
					t_xi(l) = temp_node(k1,nd)*(1-zrat2)+temp_node(k2,nd)*zrat2
					s_xi(l) = salt_node(k1,nd)*(1-zrat2)+salt_node(k2,nd)*zrat2
					ta = ta + t_xi(l)
					sa = sa + s_xi(l)
				end if
			end do
			! jw
			
			! jw
			if(icount==0) then
				write(pw_run_log,*) 'all dry'
				stop
			else if(icount/=tri_or_quad(nnel)) then
				! jw
				ta = ta/icount
				sa = sa/icount
				ttint = ta
				ssint = sa
			else
				! jw
				ta=ta/icount
				sa=sa/icount
				ttint = 0.0d0
				ssint = 0.0d0
				do l=1,tri_or_quad(nnel)
					ttint=ttint+t_xi(l)*staint(l)
					ssint=ssint+s_xi(l)*staint(l)
				end do
			end if
			
			! jw
			! jw
			index_1 = 0
			index_2 = 0
			
			! jw
			
			if(tri_or_quad(nnel)==3) then
				do l = 1, 4
					if(l <= 3) then
						! jw
						! jw
						! jw
						! jw
						node1=nodenum_at_cell(l,nnel)
						node2=facenum_at_cell(start_end_node(tri_or_quad(nnel),l,2),nnel)
						node3=facenum_at_cell(start_end_node(tri_or_quad(nnel),l,1),nnel)
						
						! jw
						
						! jw
						aa1 = calculate_area(xt,x_face(node2),x_face(node3),   &
						&                yt,y_face(node2),y_face(node3) )
						aa2 = calculate_area(x_node(node1),xt,x_face(node3),   &
						&                y_node(node1),yt,y_face(node3))
						aa3 = calculate_area(x_node(node1),x_face(node2),xt,   &
						&                y_node(node1),y_face(node2),yt)
						aa=abs(aa1)+abs(aa2)+abs(aa3)
						subrat(l)=abs(aa-area(nnel)/4)*4/area(nnel)
						if(subrat(l)<100*small_06) then
							index_1=1
							sig(1)=aa1*4/area(nnel)
							sig(2)=aa2*4/area(nnel)
							sig(1)=max(0.0_dp,min(1.0_dp,sig(1)))
							sig(2)=max(0.0_dp,min(1.0_dp,sig(2)))
							! jw
							if(sig(1)+sig(2)>1) then
								sig(3)=0
								sig(2)=1-sig(1)
							else
								sig(3)=1-sig(1)-sig(2)
							end if
							
							if(top_layer_at_node(node1)/=0 .and. top_layer_at_face(node2)/=0 .and. top_layer_at_face(node3)/=0) then
								index_2=1
								t_xi(1)=temp_node(k1,node1)*(1-zrat2)+temp_node(k2,node1)*zrat2
								t_xi(2)=temp_face(k1,node2)*(1-zrat2)+temp_face(k2,node2)*zrat2
								t_xi(3)=temp_face(k1,node3)*(1-zrat2)+temp_face(k2,node3)*zrat2
								s_xi(1)=salt_node(k1,node1)*(1-zrat2)+salt_node(k2,node1)*zrat2
								s_xi(2)=salt_face(k1,node2)*(1-zrat2)+salt_face(k2,node2)*zrat2
								s_xi(3)=salt_face(k1,node3)*(1-zrat2)+salt_face(k2,node3)*zrat2
								ttint=t_xi(1)*sig(1)+t_xi(2)*sig(2)+t_xi(3)*sig(3)
								ssint=s_xi(1)*sig(1)+s_xi(2)*sig(2)+s_xi(3)*sig(3)
								exit
							end if
						end if
					else ! jw
						node1=facenum_at_cell(1,nnel)
						node2=facenum_at_cell(2,nnel)
						node3=facenum_at_cell(3,nnel)
						aa1 = calculate_area(xt,x_face(node2),x_face(node3),   &
						    &                yt,y_face(node2),y_face(node3))
						aa2 = calculate_area(x_face(node1),xt,x_face(node3),   &
						    &                y_face(node1),yt,y_face(node3))
						aa3 = calculate_area(x_face(node1),x_face(node2),xt,   &
						    &                y_face(node1),y_face(node2),yt)
						aa  = dabs(aa1)+dabs(aa2)+dabs(aa3)
						subrat(l)=abs(aa-area(nnel)/4)*4/area(nnel)
						if(subrat(l) < 100.0d0*small_06) then
							index_1 = 1
							sig(1) = aa1*4/area(nnel)
							sig(2) = aa2*4/area(nnel)
							sig(1) = max(0.0d0,min(1.0d0,sig(1)))
							sig(2) = max(0.0d0,min(1.0d0,sig(2)))
							if(sig(1) + sig(2) > 1.0d0) then
								sig(3) = 0.0d0
								sig(2) = 1.0d0 - sig(1)
							else
								sig(3) = 1.0d0 - sig(1) - sig(2)
							end if
				! jw
							if(top_layer_at_face(node1) /=0 .and.   &
							& top_layer_at_face(node2) /=0 .and.   &
							& top_layer_at_face(node3) /=0 ) then
								index_2=1
								t_xi(1)=temp_face(k1,node1)*(1-zrat2)+temp_face(k2,node1)*zrat2
								t_xi(2)=temp_face(k1,node2)*(1-zrat2)+temp_face(k2,node2)*zrat2
								t_xi(3)=temp_face(k1,node3)*(1-zrat2)+temp_face(k2,node3)*zrat2
								s_xi(1)=salt_face(k1,node1)*(1-zrat2)+salt_face(k2,node1)*zrat2
								s_xi(2)=salt_face(k1,node2)*(1-zrat2)+salt_face(k2,node2)*zrat2
								s_xi(3)=salt_face(k1,node3)*(1-zrat2)+salt_face(k2,node3)*zrat2
								ttint = t_xi(1)*sig(1) + t_xi(2)*sig(2) + t_xi(3)*sig(3)
								ssint = s_xi(1)*sig(1) + s_xi(2)*sig(2) + s_xi(3)*sig(3)
								exit
							end if   ! jw
						end if   ! jw
					end if   ! jw
				end do ! jw
			else ! jw
				do l = 1, 4
					node1=nodenum_at_cell(l,nnel)
					node2=facenum_at_cell(start_end_node(tri_or_quad(nnel),l,3),nnel)
					node4=facenum_at_cell(start_end_node(tri_or_quad(nnel),l,2),nnel)
					aa1 = calculate_area( x_node(node1), x_face(node2), xt,   &
					    &                 y_node(node1), y_face(node2), yt    &
					    &               )
					aa2 = calculate_area( x_face(node2), x_cell(nnel), xt,  &
					    &                 y_face(node2), y_cell(nnel), yt   &
					    &               )
					aa3 = calculate_area( x_cell(nnel),x_face(node4),xt,   &
					    &                 y_cell(nnel),y_face(node4),yt)
					aa4 = calculate_area( x_face(node4),x_node(node1),xt,   &
					    &                 y_face(node4),y_node(node1),yt)
					bb1 = calculate_area(   &
					    &                 x_node(node1),   &
					    &                 x_face(node2),x_face(node4),   &
					    &                 y_node(node1),   &
					    &                 y_face(node2),y_face(node4))
					bb2 = calculate_area(   &
					    &                 x_face(node2),   &
					    &                 x_cell(nnel),x_face(node4),   &
					    &                 y_face(node2),   &
					    &                 y_cell(nnel),y_face(node4))
					
					if(bb1<=0 .or. bb2<=0) then
						write(pw_run_log,*) 'negative sub-element',nnel,bb1,bb2
						stop
					end if
					
					bb = bb1+bb2
					aa = abs(aa1)+abs(aa2)+abs(aa3)+abs(aa4)
				! jw
					subrat(l)=abs(aa-bb)/bb
					if(subrat(l) < 4.0_dp*small_06) then
						index_1=1
						call bilinear_interpolation( nnel+0.1_dp*l,x_node(node1),x_face(node2), &
						&               x_cell(nnel),x_face(node4),         &
						&               y_node(node1), y_face(node2),   &
						&               y_cell(nnel),y_face(node4),         &
						&               xt, yt, csi, etta, staint)
						if(top_layer_at_node(node1) /= 0 .and.   &
						&  top_layer_at_face(node2)   /= 0 .and.   &
						&  top_layer_at_face(node4)   /= 0) then 
							index_2  = 1
							t_xi(1) = temp_node(k1,node1)*(1-zrat2)+temp_node(k2,node1)*zrat2
							t_xi(2) = temp_face(k1,node2)*(1-zrat2)+temp_face(k2,node2)*zrat2
							t_xi(3) = ta ! jw
							t_xi(4) = temp_face(k1,node4)*(1-zrat2)+temp_face(k2,node4)*zrat2
							s_xi(1) = salt_node(k1,node1)*(1-zrat2)+salt_node(k2,node1)*zrat2
							s_xi(2) = salt_face(k1,node2)*(1-zrat2)+salt_face(k2,node2)*zrat2
							s_xi(3) = sa
							s_xi(4) = salt_face(k1,node4)*(1-zrat2)+salt_face(k2,node4)*zrat2
							
							ttint   = 0.0_dp
							ssint   = 0.0_dp
							
							do j = 1, 4
								ttint = ttint + t_xi(j)*staint(j)
								ssint = ssint + s_xi(j)*staint(j)
							end do   ! jw
						end if   ! jw
					end if ! jw
				end do ! jw
			end if
			
			! jw
			
			if(index_1==0) then
				write(pw_run_log,*) 'not in any sub-element',nnel,tri_or_quad(nnel),(subrat(i),i=1,4)
				write(pw_run_log,'(e26.16, e26.16)') xt,yt
				stop
			end if
			
			if(index_2==0 .and. icount==tri_or_quad(nnel)) then
				write(pw_run_log,*) 'all sides are not wet',nnel
			end if
		end if
	end if ! jw

end subroutine ELM_backtrace_v1

