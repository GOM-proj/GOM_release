!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine will create a rectangular regular grid for space interpolation.
!! This approach looks less stable...
!! 
subroutine create_regular_grid
	use mod_global_variables
	use mod_file_definition
	implicit none
	
	integer :: i, j, k, ii
	real(dp):: xmin, xmax, ymin, ymax, dx, dy
	! real(dp):: x_origin, y_origin
	real(dp),dimension(maxnod) :: x_node2, y_node2
	integer :: max_regular_grid_count, total_count
	real(dp):: minx, maxx, miny, maxy
	integer :: ghost_cell
	! End of local variables ==================================================!
	! include 5 ghost cells both at the beginning and end
	ghost_cell = 5
	
	x_node2 = 0.0
	y_node2 = 0.0
	
	! find the minimun and maximum values
	xmin = minval(x_node)
	xmax = maxval(x_node)
	ymin = minval(y_node)
	ymax = maxval(y_node)
	
	! define grid size in a new regular grid
	dx = 20 * 1000.0 ! 20 km
	dy = 20 * 1000.0 ! 20 km
	regular_grid_half_dx = dx * 0.5
	regular_grid_half_dy = dy * 0.5
	
	! shift original grid to (0,0)
	x_node2 = x_node - xmin
	y_node2 = y_node - ymin
	
	! total number of grids in a new regular grid
	regular_grid_xi = int((xmax - xmin)/dx) + 1 + ghost_cell * 2 ! 5 left and 5 right ghost cells
	regular_grid_yi = int((ymax - ymin)/dy) + 1 + ghost_cell * 2 ! 5 bottom and 5 top ghost cells
	
	! calculate center coordinates for a new regular grid
	allocate(regular_grid_xc(regular_grid_xi), regular_grid_yc(regular_grid_yi))
	regular_grid_xc = 0.0
	regular_grid_yc = 0.0
	do i=1,regular_grid_xi
		! regular_grid_xc(i) = dx*i - regular_grid_half_dx
		regular_grid_xc(i) = xmin + dx*(i-1) + regular_grid_half_dx - ghost_cell*dx
	end do
	do j=1,regular_grid_yi
		! regular_grid_yc(j) = dy*j - regular_grid_half_dy
		regular_grid_yc(j) = ymin + dy*(j-1) + regular_grid_half_dy - ghost_cell*dy
	end do
	
! 	write(*,*) regular_grid_xi, regular_grid_yi
! 	do i=1,regular_grid_xi
! 		write(*,*) regular_grid_xc(i)
! 	end do
! 	write(*,*)
! 	do i=1,regular_grid_yi
! 		write(*,*) regular_grid_yc(i)
! 	end do	
! 	stop 'IIIII'
	
	! find maximum node counts in each regular grid
	allocate(regular_grid_node_count(regular_grid_yi,regular_grid_xi))
	total_count = 0
	do i=1,regular_grid_xi
		do j=1,regular_grid_yi
			! grid boundary for this regular grid
			minx = regular_grid_xc(i) - regular_grid_half_dx
			maxx = regular_grid_xc(i) + regular_grid_half_dx
			miny = regular_grid_yc(j) - regular_grid_half_dy
			maxy = regular_grid_yc(j) + regular_grid_half_dy
			
			! find how many nodes belongs to this grid cell			
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

	! find maximum node count
	max_regular_grid_count = maxval(regular_grid_node_count)

! 	write(*,*) 'total count = ', total_count, 'max_count = ', max_regular_grid_count
! 	do j=1,regular_grid_yi
! 		write(*,'(*(I5))') (regular_grid_node_count(j,i), i=1,regular_grid_xi)
! 	end do	

	if(total_count /= maxnod) then
		write(pw_run_log,*) 'Error: create_regular_grid.f90: '
		write(pw_run_log,*) 'Total node numbers in a regular grid should be equal to maxnod'
		write(*,*) 'Error: create_regular_grid.f90: '
		write(*,*) 'Total node numbers in a regular grid should be equal to maxnod'
		stop
	end if
! 	stop 'I am here'
	
	! Now, allocate regular_grid
	allocate(regular_grid(regular_grid_yi,regular_grid_xi,max_regular_grid_count))
	regular_grid = -999
	
	! Now, store node information which belongs each regular grid
	do i=1,regular_grid_xi
		do j=1,regular_grid_yi
			! grid boundary for this regular grid
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