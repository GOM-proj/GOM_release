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
	integer, intent(inout) 	:: nnel,jlev 	! initial element and vertical level, then output with destination element and vertical level 
	real(dp),intent(inout) 	:: xt,yt,zt 	! final destination coordinate
	
	integer :: k, l, k1, k2
	integer :: nel, node1, node2, jd1, jd2, iflag, nel_j, &
	& 			  md1, md2, iter_temp, isd, iteration
	
	real(dp):: trm, area1, area2, ae, xcg, ycg, pathl, xin, yin, &
	& 				zin, dist, xvel, yvel, zvel, hvel 
	! End of local variables ==================================================!
	
	! (0) before start, check if the current element is dry or not. ===========!
	! we cannot perform backtracking from the dry element!
	if(top_layer_at_element(nnel) == 0 .or. &
	&	MSL-h_cell(nnel) > MSL+eta_cell(nnel)) then 
		! Here, MSL-h_cell(nnel) = bottom elevation
		!       MSL+eta_cell(nnel) = water surface elevation at current time
		! i.e., if bottom elevation is above the water surface elevation, 
		! that means the element turns into a dry cell.
	   write(pw_run_log,*) 	'search_straight_line.f90: starting element is dry: stop'
	   write(*,*) 				'search_straight_line.f90: starting element is dry: stop'
	   stop
	end if
	
	
	! (1) Horizontal element finding ==========================================!
	! Here, I will find the element that encompass the final destination point (xt,yt,zt)
	nfl = 0
	trm = time 	! (= bt_dt = dt/bt_step) i.e., sub time in [sec]
	
	! (1.1) Check if the given initial point is really located in the given initial element
	! calculate sub triangular areas with each node and given (x0,y0) or (xt,yt)
	! if (x0,y0) locates in this element (including on the face), the calculated area, area1, must be the same as the element area.
	! if (xt,yt) locates outside the element, the calculated area must differ from the current element area.
	! So, this process is to find if the (x0,y0) and (xt,yt) locate in the initial element.
	nel = nnel 	! nel = starting element; from now: nel = staring element, nnel = final element
	area1 = 0.0
	area2 = 0.0
	do l = 1, tri_or_quad(nel)
		! node1 starts from the second counterclockwise direction node from (x0,y0)
		! node2 is the one more counterclockwise direction node
		node1 = nodenum_at_cell(l,nel)
		node2 = nodenum_at_cell(start_end_node(tri_or_quad(nel),l,1),nel)
		
		! function calculate_area(x1,x2,x3,y1,y2,y3) will return a triangular area with three node information: 
		! i.e., (x1,y1), (x2,y2), (x3,y3)
		! area1: area with the initial point
		! area2: area with the final point
		area1 = area1 + abs( calculate_area(x_node(node1), x_node(node2), x0,   &
		&                                   y_node(node1), y_node(node2), y0))
		area2 = area2 + abs( calculate_area(x_node(node1), x_node(node2), xt,   &
		&                                   y_node(node1), y_node(node2), yt))
	end do

	! (1.2) Check if the given terminal point is located in the given initial element		
	! if calculated area is different from the current element's area,
	! that means the (x0,y0) is not in the current element.
	! Note, here we have to use "small" like value since we are working with real number. 
	ae = abs(area1 - area(nel))/area(nel) ! calculate area difference respect to the current element's area
	if(ae > small_06) then
		write(pw_run_log,*) 	'search_straight_line.f90: (x0,y0) not in the current element (nnel) initially: stop'
		write(pw_run_log,*) 	'calculaated area with (x0,y0) = ', area1, ', current element area = ', area(nnel)
		write(*,*) 				'search_straight_line.f90: (x0,y0) not in the current element (nnel) initially: stop'
		write(*,*) 				'calculaated area with (x0,y0) = ', area1, ', current element area = ', area(nnel)
		stop
	end if
	
	! same as above, but with (xt,yt), i.e., backtracted point
	! if the area difference is small enough, it means the point (xt,yt) locates in the current element,
	! and that means we found the ending element, and that is the current element.
	! if the area difference is big engough, it means the point (xt,yt) locates outside the current element,
	! and that means we have to find the ending element that encompass (xt,yt)
	ae = abs(area2 - area(nel))/area(nel)
	if(ae < small_06) then
		! update the endding element (nnel) with the starting element (nel)
		! now, we found the terminal element, so go to the vertical position (layer) finding section
		nnel = nel
	else
		! Below is for the case of (ae > small_06):
		! the terminal point is not located in the starting element, 
		! thus we should find the terminal element.
		! (xt,yt) is not in the starting element, i.e. nel, and thus (x0,y0) and (xt,yt) are distinctive
		! If needed, let's slightly modify the initial point to prevent underflow for iloc >=1.
		
		! This was the original approach for nudging...
 		! if(i_which_backtrack == 1) then ! nudge initial point since the variable is located at face
			! nudge the initial point slightly toward the cell center
			! (xcg, ycg) will be the new starting point:
			! 		if x0 < x_cell(nel), shift a bit to right (i.e., toward x_cell(nel))
			! 		if x0 > x_cell(nel), shift a bit to left (i.e., toward x_cell(nel)
 		! 	xcg = (1.0-1.0d-4)*x0 + 1.0d-4*x_cell(nel)
		! 	ycg = (1.0-1.0d-4)*y0 + 1.0d-4*y_cell(nel)
		! else if(i_which_backtrack == 2) then ! do not nudge initial point since we are using velocity at cell center
		! 	xcg = x0
		! 	ycg = y0			
		! end if
		
		! This is the new approach (just turn on nudging every case...
 		xcg = (1.0-1.0d-4)*x0 + 1.0d-4*x_cell(nel)
		ycg = (1.0-1.0d-4)*y0 + 1.0d-4*y_cell(nel)
		
		! calculate the distance between the (updated) starting point (xcg,ycg) and the backtracked point (xt,yt)
		pathl = sqrt((xt-xcg)**2+(yt-ycg)**2) ! horizontal distance between the staring point (xcg,ycg) and terminal point (xt,yt)
		! actually, we don't need to check this again since they are already in different elements.
		! but anyway, let's just keep this for the double check. 
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
		
		! if these two points are in different elements, 
		! then the straight line must pass one of the initial element's side (face).
		! So, let's find the initial crossing face, nel_j. 
		do l = 1, tri_or_quad(nel)
			jd1 = nodenum_at_cell(start_end_node(tri_or_quad(nel),l,1),nel) ! current element's lth face's starting node
			jd2 = nodenum_at_cell(start_end_node(tri_or_quad(nel),l,2),nel) ! current element's lth face's ending node
			
			! check_intersection will return:
			! 		iflag: whether intersection point exists (= 1) or not (= 0)
			! 		(xin,yin): intersection point
			call check_intersection(xcg,xt,x_node(jd1),x_node(jd2), &
			&								ycg,yt,y_node(jd1),y_node(jd2), &
			&								iflag,xin,yin)
			
			if(iflag == 1) then
				! intersection is found between the starting/ending line and the current side, l
				nel_j = l
				exit ! exit do loop
			end if
		end do
		
		if(iflag == 0) then
			! if there is no intersection, that means something is wrong; so stop.
			! actually, this is also redundant, but let's keep it to make sure.
			write(pw_run_log,*) 	'search_straigth_line.f90: no intersecting edge was found: stop'
			write(*,*) 				'search_straigth_line.f90: no intersecting edge was found: stop'
			! write(*,*) it, nel, xcg, xt, ycg, yt, xin, yin
			stop
		end if
		! ======================================================================!
		
		! Finding new element for (xt,yt) ======================================!
		! until now, we just checked if the terminal point really locates outside this element.
		! from this point, we are sure the terminal point really locates outside this element.
		! So, we will find the new element information.
		
		! Following loop will find the horizontal element, nnel, which encompass the final backtracted point (xt,yt)
		! if the element is found, exit the loop and go to the vertical position finding section. 
		zin = z0 ! set the initial vertical position to the initial z0 (which locates at the middle of the vertical layer)
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
			
			! backtrack from the intersecting point, (xin,yin), and face, nel_j
			md1 = nodenum_at_cell(start_end_node(tri_or_quad(nel),nel_j,1),nel) ! starting element's intersecting face's starting node
			md2 = nodenum_at_cell(start_end_node(tri_or_quad(nel),nel_j,2),nel) ! starting element's intersecting face's ending node
		
			dist = sqrt((xin-xt)**2+(yin-yt)**2) ! horizontal distance between the intersecting point and the terminal point
			if(dist/pathl > 1.0+1.0d-4) then
				! i.e., dist cannot bigger than path1
				! dist should be slightly shorter than (if nudge was used)
				! or equal to the total path line length (path1)
				! (if nudge was not used, i.e., when (x0,y0) is used, and (x0,y0) starts from the face as for the velocity backtracking) 
				write(pw_run_log,*) 	'search_straight_line.f90: path overshot: stop'
				write(*,*) 				'search_straight_line.f90: path overshot: stop'
				stop
			end if
			
			! now, let's compute z position
			! new zin position (i.e., vertical intersection position) with linear interpolation, i.e.,:
			! (dist/path1) = (zt - zin_new)/(zt - zin_old)
			! Thus, zin_new = zt - dist/path1 * (zt - zin_old)
			zin = zt - dist/pathl*(zt-zin)
			
			! update remained time and path1
			trm = trm*dist/pathl ! time remaining [sec] from the intersection to the final point
			pathl = sqrt((xin-xt)**2+(yin-yt)**2) ! new path1: from intersection to the final point
			if(pathl==0.0 .or. trm==0.0) then
				write(pw_run_log,*) 	'target reached'
				write(*,*)				'target reached'
				stop
			end if
			
			iter_temp = 0 ! flag
			
			! abnormal cases control --------------------------------------------!
			! Now, let's check if the terminal point locates in a dry element or off the boundary.
			! if so, we have to relocate the terminal point to the intersecting boundary element face instead the original terminal point.
			! i.e., when the face is located at the boundary or the adjacent element is dry:
			! for horizontal exit and dry elements, compute tangential vel.,
			! update target (xt,yt,zt) and continue.
			if(adj_cellnum_at_cell(nel_j,nel)==0  .or.  & 							! if terminal point is in a land cell or over the boundary element
			& 	top_layer_at_element(adj_cellnum_at_cell(nel_j,nel))==0) then 	! if the terminal point is in dry element
				iter_temp = 1 ! we found the terminal element; set this to "1" to exit loop
				isd = facenum_at_cell(nel_j,nel) ! terminal point is on the face "nel_j" of the element "nel", i.e., the face ID is "isd"
				
				if(nodenum_at_face(1,isd)+nodenum_at_face(2,isd)/=md1+md2) then ! jw, I think this criteria is not clear enough; this is the original approach.
				! if(nodenum_at_face(1,isd) /= md1 .or. nodenum_at_face(2,isd) /= md2) then ! jw, I think this criteria is more clear, check this again.
					write(pw_run_log,*)	'search_straight_line.f90: wrong side: stop'
					write(*,*)				'search_straight_line.f90: wrong side: stop'
					stop
				end if
				
				! nudge intersect (xin,yin) point backward, and update starting pt
				! This is to make sure this final point is in the element.
				xin =(1.0-1.0d-4)*xin+1.0d-4*x_cell(nel)
				yin =(1.0-1.0d-4)*yin+1.0d-4*y_cell(nel)
				xcg = xin
				ycg = yin
				
				xvel = -vn_face(jlev,isd)*sin_theta(isd) ! u = u*cos - v*sin, then u* = 0.0 at the land boundary face
				yvel =  vn_face(jlev,isd)*cos_theta(isd) ! v = u*sin + v*cos, then u* = 0.0 at the land boundary face
				zvel = (w_node(jlev,md1)+w_node(jlev,md2))*0.5_dp
				xt   = xin-xvel*trm
				yt   = yin-yvel*trm
				zt   = zin-zvel*trm
				hvel = sqrt(xvel**2+yvel**2)
				
				! this was the original approach, but I am not sure about this.
				! so, I update this part...
				! from here...
				! if(hvel<1.e-4) then ! jw, I am not sure about this criteria... (shouldn't we check un_face? not the the horizontal velocity?)
					! we found the new horizontal terminal position, i.e., not original (xt,yt,zt) but with (xin,yin,zin),
					! and the terminal element is "nel", i.e, the open boundary element or dry land element's neighbor.
				! 	nfl=1 
				! 	xt=xin
				! 	yt=yin
				! 	zt=zin
				! 	nnel=nel
				! 	exit search_loop
				! end if
				! pathl=hvel*trm
				! to here...
				
				! this is the new approach...
				! i.e., we just found the terminal element, so just exit the loop
			 	nfl=1 
			 	xt=xin
			 	yt=yin
			 	zt=zin
			 	nnel=nel
			 	exit search_loop				
			end if 
			! end of abnormal cases control -------------------------------------!
	
			! search for nel's neighbor with edge nel_j, or in abnormal cases, the same
			! update nel to the next element, which is in the backtracking direction
			! Then, find if the final backtracted point (xt,yt) locates in this new element.
			! the following routin is exactly same as the routine in the very beginning.
			! so, do this until final (xt,yt) is located in the new "nel"
			if(iter_temp == 0) then ! i.e., we still couldn't find the terminal element, so we should keep doing backtracking.
				! update "nel" to the new element that is in backtracking direction.
				nel = adj_cellnum_at_cell(nel_j,nel) ! next front element
			end if
			
			! claculate area with the terminal points (xt,yt) and the nodes in a new "nel"
			! again, if (xt,yt) is located in the new "nel", that is the final element we are looking for.
			area1 = 0.0
			do l=1,tri_or_quad(nel)
				k1=nodenum_at_cell(l,nel)
				k2=nodenum_at_cell(start_end_node(tri_or_quad(nel),l,1),nel)
				area1 = area1 + dabs(calculate_area(x_node(k1),x_node(k2),xt, &
				&												y_node(k1),y_node(k2),yt))
			end do
			
			! calculate area difference respect to the current element's area
			ae = abs(area1-area(nel))/area(nel)
	
			! if calculated area is different from the current element's area,
			! that means the (xt,yt) is not in the current element.
			! Note, here we have to use "small" like value since we are working with real number. 			
			if(ae < small_06) then
				! update the ending element (nnel) to the current element (nel)
				! now, we found the terminal element, so exit this element searching loop
				! and go to the vertical position (layer) finding section			
				! what we found here:
				! 		final element: "nel"
				! 		termianl point: (xt,yt,zt) in normal water cell
				nnel = nel
				exit search_loop
			else
				! we haven't yet reach to the terminal element.
				! i.e., the terminal point locates outside this element.
				! Thus, we have to check this search loop again.			
				! To do so, let's find the next intersecting edge
				! then, we will start this search loop again the final point reach to the either "abnormal cases" or "ae < small_06".
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
				   	nel_j = l ! next front edge          
				   	cycle search_loop
				   end if
				end do
			end if
			
		 	if(iflag == 0) then
				! if there is no intersection, that means something is wrong; so stop.
				! actually, this is also redundant, but let's keep it to make sure.		 		
				write(pw_run_log,*)	'search_straight_line.f90: failed to find the next edge: stop', &
				&							iter_temp,xin,yin,xt,yt,nel,  &
				&							md1,md2,idt,  &
				&							nodenum_at_face(1,j_face), &
				&							nodenum_at_face(2,j_face)
				write(*,*)				'search_straight_line.f90: failed to find the next edge: stop'
				stop
			end if
		end do search_loop
	end if ! if(ae < small_06) then
	! End of horizontal element finding =======================================!
	
	
	! =========================================================================!
	! Until now, I have found the horizontal element, nnel, which encompass the destination point (xt,yt).
	! Now, I will find the vertical level for the final destination point (xt,yt,zt); (after this sub elm step).
	! note: nnel is the terminal (ending) element
	if(top_layer_at_element(nnel) == 0 .or. &
	& 	MSL-h_cell(nnel) > MSL+eta_cell(nnel)) then
		! Here, MSL-h_cell(nnel) = bottom elevation
		!       MSL+eta_cell(nnel) = water surface elevation at current time
		! i.e., if bottom elevation is above the water surface elevation, 
		! that means the element turns into a dry cell.
	   write(pw_run_log,*) 	'search_straight_line.f90: ending element is dry: stop'
	   write(*,*) 				'search_straight_line.f90: ending element is dry: stop'
		stop
	end if
	
	! find the vertical position
	! zt must locate between bottom elevation and the surface elevation, i.e., (bottom_elevation <= zt <= surface_elevation)
	zt = MIN(MAX(zt, MSL-h_cell(nnel)), MSL+eta_cell(nnel))
	
	! find the vertical layer
	if(zt <= z_level(0)) then
		write(pw_run_log,*) 'search_straight_line.f90: ilegal vertical position: exit 1!'
		stop
	else if(zt > z_level(maxlayer)) then
		write(pw_run_log,*) 'search_straight_line.f90: ilegal vertical position: exit 2!'
		stop
	else
		! if zt locates in the valid vertical water column, find the vertical layer.
		do k=0,maxlayer-1
			if(zt > z_level(k) .and. zt <= z_level(k+1)) then
				jlev = k+1 ! this is vertical layer not level
				exit
			end if
		end do
	end if
	
	! bring jlev into local range (due to small inconsistency in defining top_layer_at_element and bottom_layer_at_element)
	jlev = MIN(MAX(jlev,bottom_layer_at_element(nnel)), top_layer_at_element(nnel))
	
	! Finally, we found (xt,yt,zt), nnel, and jlev ============================!
end subroutine search_straight_line

