!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Read initial temperature at horizontally nodal and vertically each center of vertical layers.
!! Then, calculate temperature at face and cell center.
!! 
subroutine read_temp_init
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, j, k, l
	integer :: n1, n2
	real(dp):: sum1
	integer :: serial_num ! node_id
	! End of local variables ==================================================!

	open(pw_temp_init,file=id_temp_init,form='formatted',status='old')

	! skip header lines
	call skip_header_lines(pw_temp_init,id_temp_init)
	
	! read main body ==========================================================!
	do i=1,maxnod
		read(pw_temp_init,*) serial_num, (temp_node(k,i), k=1,maxlayer)
	end do
	close(pw_temp_init)
	
	
	! Calculate temperature at face ===========================================!
	do j=1,maxface
		n1 = nodenum_at_face(1,j)
		n2 = nodenum_at_face(2,j)
		do k=1,maxlayer
			temp_face(k,j) = (temp_node(k,n1) + temp_node(k,n2)) * 0.5
		end do
	end do
	
	! Calculate temperature at cell ===========================================!
	! note: here, I simply use arithmetic mean, but later I may need to update this to volume weighted mean (not really important for the initial condition...)
	do i=1,maxele
		do k=1,maxlayer
			sum1 = 0.0
			do l=1,tri_or_quad(i)
				n1 = nodenum_at_cell(l,i)
				sum1 = sum1 + temp_node(k,n1)
			end do
			temp_cell(k,i) = sum1/tri_or_quad(i)
		end do
	end do
	! initially set temp_cell_new to old value
	temp_cell_new = temp_cell
	
	! check...
	! write(*,*) 'I am here'
! 	do i=1,3
! 	 	write(*,*) (temp_node(k,i), k=1,maxlayer)
! 	 	write(*,*) (temp_cell(k,i), k=1,maxlayer)
! 	 	write(*,*) (temp_face(k,i), k=1,maxlayer)
! 	 	write(*,*)
! 	 	write(*,*) (temp_cell_new(k,i), k=1,maxlayer)
! 	end do
end subroutine read_temp_init