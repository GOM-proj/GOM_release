!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Read initial salinity at horizontally nodal and vertically each center of vertical layers.
!! Then, calculate salinity at face and cell center.
!! 
subroutine read_salt_init
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, j, k, l
	integer :: n1, n2
	real(dp):: sum1
	integer :: serial_num ! jw
	! jw

	open(pw_salt_init,file=id_salt_init,form='formatted',status='old')

	! jw
	call skip_header_lines(pw_salt_init,id_salt_init)
	
	! jw
	do i=1,maxnod
		read(pw_salt_init,*) serial_num, (salt_node(k,i), k=1,maxlayer)
	end do
	close(pw_salt_init)
	
	
	! jw
	do j=1,maxface
		n1 = nodenum_at_face(1,j)
		n2 = nodenum_at_face(2,j)
		do k=1,maxlayer
			salt_face(k,j) = (salt_node(k,n1) + salt_node(k,n2)) * 0.5
		end do
	end do
	
	! jw
	! jw
	do i=1,maxele
		do k=1,maxlayer
			sum1 = 0.0
			do l=1,tri_or_quad(i)
				n1 = nodenum_at_cell(l,i)
				sum1 = sum1 + salt_node(k,n1)
			end do
			salt_cell(k,i) = sum1/tri_or_quad(i)
		end do
	end do
	! jw
	salt_cell_new = salt_cell	
end subroutine read_salt_init