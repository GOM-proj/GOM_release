!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!!
!! straightline search algorithm. 
!! Initially nnel is an element that encompasses (x0,y0).
!!     iloc=0: do not nudge initial pointt
!!     iloc=1: nudge.	 
!!     input : iloc,nnel,x0,y0,z0,xt,yt,zt,jlev,time,and vn_face,w_node for
!!             abnormal cases;
!!     output: the updated end point (xt,yt,zt) (if so), nnel, jlev, a flag nfl.	
!!     exit btrack if a bnd or dry element is hit and vel. there is small, or death trap is reached.
!! ===========================================================================!
!! call    search_straight_line(1,   nnel,jlev,bt_dt,x0,y0,z0,xt,yt,zt,iflqs1,idt,j_face) ! this is the corresponding "call" in ELM_bactrace.f90
subroutine search_straight_line(i_which_backtrack,nnel,jlev, time,x0,y0,z0,xt,yt,zt,   nfl,idt,j_face)
	use mod_global_variables
	use mod_file_definition
	use mod_function_library, only : calculate_area
	implicit none
	
	integer, intent(in) 		:: i_which_backtrack,idt,j_face
	real(dp),intent(in) 		:: time,x0,y0,z0
	integer, intent(out) 	:: nfl
	integer, intent(inout) 	:: nnel,jlev 	! jw
	real(dp),intent(inout) 	:: xt,yt,zt 	! jw
	
	integer :: k, l, k1, k2
	integer :: nel, node1, node2, jd1, jd2, iflag, nel_j, &
	& 			  md1, md2, iter_temp, isd, iteration
	
	real(dp):: trm, area1, area2, ae, xcg, ycg, pathl, xin, yin, &
	& 				zin, dist, xvel, yvel, zvel, hvel 
	! jw
	
	! jw
	! jw
	if(top_layer_at_element(nnel) == 0 .or. &
	&	MSL-h_cell(nnel) > MSL+eta_cell(nnel)) then 
		! jw
		! jw
		! jw
		! jw
	   write(pw_run_log,*) 	'search_straight_line.f90: starting element is dry: stop'
	   write(*,*) 				'search_straight_line.f90: starting element is dry: stop'
	   stop
	end if
	
	
	! jw
	! jw
	nfl = 0
	trm = time 	! jw
	
	! jw
	! jw
	! jw
	! jw
	! jw
	nel = nnel 	! jw
	area1 = 0.0
	area2 = 0.0
	do l = 1, tri_or_quad(nel)
		! jw
		! jw
		node1 = nodenum_at_cell(l,nel)
		node2 = nodenum_at_cell(start_end_node(tri_or_quad(nel),l,1),nel)
		
		! jw
		! jw
		! jw
		! jw
		area1 = area1 + abs( calculate_area(x_node(node1), x_node(node2), x0,   &
		&                                   y_node(node1), y_node(node2), y0))
		area2 = area2 + abs( calculate_area(x_node(node1), x_node(node2), xt,   &
		&                                   y_node(node1), y_node(node2), yt))
	end do

	! jw
	! jw
	! jw
	! jw
	ae = abs(area1 - area(nel))/area(nel) ! jw
	if(ae > small_06) then
		write(pw_run_log,*) 	'search_straight_line.f90: (x0,y0) not in the current element (nnel) initially: stop'
		write(pw_run_log,*) 	'calculaated area with (x0,y0) = ', area1, ', current element area = ', area(nnel)
		write(*,*) 				'search_straight_line.f90: (x0,y0) not in the current element (nnel) initially: stop'
		write(*,*) 				'calculaated area with (x0,y0) = ', area1, ', current element area = ', area(nnel)
		stop
	end if
	
	! jw
	! jw
	! jw
	! jw
	! jw
	ae = abs(area2 - area(nel))/area(nel)
	if(ae < small_06) then
		! jw
		! jw
		nnel = nel
	else
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
 		xcg = (1.0-1.0d-4)*x0 + 1.0d-4*x_cell(nel)
		ycg = (1.0-1.0d-4)*y0 + 1.0d-4*y_cell(nel)
		
		! jw
		pathl = sqrt((xt-xcg)**2+(yt-ycg)**2) ! jw
		! jw
		! jw
		if((xcg == xt .and. ycg == yt) .or. pathl == 0.0) then
			write(pw_run_log,*) 	'search_straight_line.f90: backtracked path has zero length: stop:'
			write(pw_run_log,*) 	'(x0,y0):   ', x0, y0
			write(pw_run_log,*) 	'(xcg,ycg): ', xcg, ycg
			write(pw_run_log,*) 	'(xt,yt):   ', xt, yt
			write(*,*) 				'search_straight_line.f90: backtracked path has zero length: stop:'
			write(*,*) 				'(x0,y0):   ', x0, y0
			write(*,*) 				'(xcg,ycg): ', xcg, ycg
			write(*,*) 				'(xt,yt):   ', xt, yt
			stop
		end if
		
		! jw
		! jw
		! jw
		do l = 1, tri_or_quad(nel)
			jd1 = nodenum_at_cell(start_end_node(tri_or_quad(nel),l,1),nel) ! jw
			jd2 = nodenum_at_cell(start_end_node(tri_or_quad(nel),l,2),nel) ! jw
			
			! jw
			! jw
			! jw
			call check_intersection(xcg,xt,x_node(jd1),x_node(jd2), &
			&								ycg,yt,y_node(jd1),y_node(jd2), &
			&								iflag,xin,yin)
			
			if(iflag == 1) then
				! jw
				nel_j = l
				exit ! jw
			end if
		end do
		
		if(iflag == 0) then
			! jw
			! jw
			write(pw_run_log,*) 	'search_straigth_line.f90: no intersecting edge was found: stop'
			write(*,*) 				'search_straigth_line.f90: no intersecting edge was found: stop'
			! jw
			stop
		end if
		! jw
		
		! jw
		! jw
		! jw
		! jw
		
		! jw
		! jw
		zin = z0 ! jw
		iteration = 0		
		search_loop: do
			iteration = iteration + 1
			
			if(iteration > 1000) then
				write(pw_run_log,*) 'search_straight_line.f90: death trap reached', idt, j_face
				nfl = 1
				xt  = xin
				yt  = yin
				zt  = zin
				nnel= nel
				exit search_loop
			end if
			
			! jw
			md1 = nodenum_at_cell(start_end_node(tri_or_quad(nel),nel_j,1),nel) ! jw
			md2 = nodenum_at_cell(start_end_node(tri_or_quad(nel),nel_j,2),nel) ! jw
		
			dist = sqrt((xin-xt)**2+(yin-yt)**2) ! jw
			if(dist/pathl > 1.0+1.0d-4) then
				! jw
				! jw
				! jw
				! jw
				write(pw_run_log,*) 	'search_straight_line.f90: path overshot: stop'
				write(*,*) 				'search_straight_line.f90: path overshot: stop'
				stop
			end if
			
			! jw
			! jw
			! jw
			! jw
			zin = zt - dist/pathl*(zt-zin)
			
			! jw
			trm = trm*dist/pathl ! jw
			pathl = sqrt((xin-xt)**2+(yin-yt)**2) ! jw
			if(pathl==0.0 .or. trm==0.0) then
				write(pw_run_log,*) 	'target reached'
				write(*,*)				'target reached'
				stop
			end if
			
			iter_temp = 0 ! jw
			
			! jw
			! jw
			! jw
			! jw
			! jw
			! jw
			if(adj_cellnum_at_cell(nel_j,nel)==0  .or.  & 							! jw
			& 	top_layer_at_element(adj_cellnum_at_cell(nel_j,nel))==0) then 	! jw
				iter_temp = 1 ! jw
				isd = facenum_at_cell(nel_j,nel) ! jw
				
				if(nodenum_at_face(1,isd)+nodenum_at_face(2,isd)/=md1+md2) then ! jw
				! jw
					write(pw_run_log,*)	'search_straight_line.f90: wrong side: stop'
					write(*,*)				'search_straight_line.f90: wrong side: stop'
					stop
				end if
				
				! jw
				! jw
				xin =(1.0-1.0d-4)*xin+1.0d-4*x_cell(nel)
				yin =(1.0-1.0d-4)*yin+1.0d-4*y_cell(nel)
				xcg = xin
				ycg = yin
				
				xvel = -vn_face(jlev,isd)*sin_theta(isd) ! jw
				yvel =  vn_face(jlev,isd)*cos_theta(isd) ! jw
				zvel = (w_node(jlev,md1)+w_node(jlev,md2))*0.5_dp
				xt   = xin-xvel*trm
				yt   = yin-yvel*trm
				zt   = zin-zvel*trm
				hvel = sqrt(xvel**2+yvel**2)
				
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
			 	nfl=1 
			 	xt=xin
			 	yt=yin
			 	zt=zin
			 	nnel=nel
			 	exit search_loop				
			end if 
			! jw
	
			! jw
			! jw
			! jw
			! jw
			! jw
			if(iter_temp == 0) then ! jw
				! jw
				nel = adj_cellnum_at_cell(nel_j,nel) ! jw
			end if
			
			! jw
			! jw
			area1 = 0.0
			do l=1,tri_or_quad(nel)
				k1=nodenum_at_cell(l,nel)
				k2=nodenum_at_cell(start_end_node(tri_or_quad(nel),l,1),nel)
				area1 = area1 + dabs(calculate_area(x_node(k1),x_node(k2),xt, &
				&												y_node(k1),y_node(k2),yt))
			end do
			
			! jw
			ae = abs(area1-area(nel))/area(nel)
	
			! jw
			! jw
			! jw
			if(ae < small_06) then
				! jw
				! jw
				! jw
				! jw
				! jw
				! jw
				nnel = nel
				exit search_loop
			else
				! jw
				! jw
				! jw
				! jw
				! jw
				do l=1,tri_or_quad(nel)
					jd1 = nodenum_at_cell(start_end_node(tri_or_quad(nel),l,1),nel)
					jd2 = nodenum_at_cell(start_end_node(tri_or_quad(nel),l,2),nel)
					if(jd1==md1 .and. jd2==md2 .or. jd2==md1 .and. jd1==md2) then
						cycle
					end if
					call check_intersection(xcg,xt,x_node(jd1),x_node(jd2), &
					&								ycg,yt,y_node(jd1),y_node(jd2), &
					&								iflag,xin,yin)
				   if(iflag == 1) then
				   	nel_j = l ! jw
				   	cycle search_loop
				   end if
				end do
			end if
			
		 	if(iflag == 0) then
				! jw
				! jw
				write(pw_run_log,*)	'search_straight_line.f90: failed to find the next edge: stop', &
				&							iter_temp,xin,yin,xt,yt,nel,  &
				&							md1,md2,idt,  &
				&							nodenum_at_face(1,j_face), &
				&							nodenum_at_face(2,j_face)
				write(*,*)				'search_straight_line.f90: failed to find the next edge: stop'
				stop
			end if
		end do search_loop
	end if ! jw
	! jw
	
	
	! jw
	! jw
	! jw
	! jw
	if(top_layer_at_element(nnel) == 0 .or. &
	& 	MSL-h_cell(nnel) > MSL+eta_cell(nnel)) then
		! jw
		! jw
		! jw
		! jw
	   write(pw_run_log,*) 	'search_straight_line.f90: ending element is dry: stop'
	   write(*,*) 				'search_straight_line.f90: ending element is dry: stop'
		stop
	end if
	
	! jw
	! jw
	zt = MIN(MAX(zt, MSL-h_cell(nnel)), MSL+eta_cell(nnel))
	
	! jw
	if(zt <= z_level(0)) then
		write(pw_run_log,*) 'search_straight_line.f90: ilegal vertical position: exit 1!'
		stop
	else if(zt > z_level(maxlayer)) then
		write(pw_run_log,*) 'search_straight_line.f90: ilegal vertical position: exit 2!'
		stop
	else
		! jw
		do k=0,maxlayer-1
			if(zt > z_level(k) .and. zt <= z_level(k+1)) then
				jlev = k+1 ! jw
				exit
			end if
		end do
	end if
	
	! jw
	jlev = MIN(MAX(jlev,bottom_layer_at_element(nnel)), top_layer_at_element(nnel))
	
	! jw
end subroutine search_straight_line

