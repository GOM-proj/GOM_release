!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
subroutine create_regular_grid
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, j, k, ii
	real(dp):: xmin, xmax, ymin, ymax, dx, dy, x_origin, y_origin
	real(dp),allocatable :: x_node2(:), y_node2(:)
	integer :: max_regular_grid_count, total_count
	real(dp):: minx, maxx, miny, maxy
	integer :: x_ghost_cell, y_ghost_cell
	! jw
	! jw
	x_ghost_cell = 5
	y_ghost_cell = 5
	
	allocate(x_node2(maxnod), y_node2(maxnod))
	x_node2 = 0.0
	y_node2 = 0.0
	
	! jw
	xmin = minval(x_node)
	xmax = maxval(x_node)
	ymin = minval(y_node)
	ymax = maxval(y_node)
	
	! jw
	dx = 20 * 1000.0 ! jw
	dy = 20 * 1000.0 ! jw
	regular_grid_half_dx = dx * 0.5
	regular_grid_half_dy = dy * 0.5
	
	! jw
! jw
! jw
	
	! jw
	regular_grid_xi = int((xmax - xmin)/dx) + 1 + x_ghost_cell * 2 ! jw
	regular_grid_yi = int((ymax - ymin)/dy) + 1 + y_ghost_cell * 2 ! jw
	
	! jw
	allocate(regular_grid_xc(regular_grid_xi), regular_grid_yc(regular_grid_yi))
	regular_grid_xc = 0.0
	regular_grid_yc = 0.0
	do i=1,regular_grid_xi
		! jw
		regular_grid_xc(i) = xmin + dx*(i-1) + regular_grid_half_dx - x_ghost_cell*dx
	end do
	do j=1,regular_grid_yi
		! jw
		regular_grid_yc(j) = ymin + dy*(j-1) + regular_grid_half_dy - y_ghost_cell*dy
	end do
	
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
	allocate(regular_grid_node_count(regular_grid_yi,regular_grid_xi))
	total_count = 0
	do i=1,regular_grid_xi
		do j=1,regular_grid_yi
			! jw
			minx = regular_grid_xc(i) - regular_grid_half_dx
			maxx = regular_grid_xc(i) + regular_grid_half_dx
			miny = regular_grid_yc(j) - regular_grid_half_dy
			maxy = regular_grid_yc(j) + regular_grid_half_dy
			
			! jw
			k = 0
			do ii=1,maxnod
				if(x_node(ii) >= minx .and. x_node(ii) < maxx) then
					if(y_node(ii) >= miny .and. y_node(ii) < maxy) then
						k = k+1
					end if
				end if
			end do
			total_count = total_count + k
			regular_grid_node_count(j,i) = k
		end do
	end do

	! jw
	max_regular_grid_count = maxval(regular_grid_node_count)
	
! jw
! jw
! jw
! jw

	if(total_count /= maxnod) then
		write(pw_run_log,*) 'Error: create_regular_grid.f90: '
		write(pw_run_log,*) 'Total node numbers in a regular grid should be equal to maxnod'
		write(*,*) 'Error: create_regular_grid.f90: '
		write(*,*) 'Total node numbers in a regular grid should be equal to maxnod'
		stop
	end if
!	stop 'I am here'
	
	! jw
	allocate(regular_grid(regular_grid_yi,regular_grid_xi,max_regular_grid_count))
	regular_grid = -999
	
	! jw
	do i=1,regular_grid_xi
		do j=1,regular_grid_yi
			! jw
			minx = regular_grid_xc(i) - regular_grid_half_dx
			maxx = regular_grid_xc(i) + regular_grid_half_dx
			miny = regular_grid_yc(j) - regular_grid_half_dy
			maxy = regular_grid_yc(j) + regular_grid_half_dy
			
			k = 0
			do ii=1,maxnod
				if(x_node(ii) >= minx .and. x_node(ii) < maxx) then
					if(y_node(ii) >= miny .and. y_node(ii) < maxy) then
						k = k+1
						regular_grid(j,i,k) = ii
					end if
				end if
			end do
		end do
	end do
	
end subroutine create_regular_grid