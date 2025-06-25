!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Set initial density analytically
!! 
subroutine ana_initial_density
	use mod_global_variables
	implicit none
	
	integer :: i, j, k, l, n1, n2
	real(dp):: sum1
	! End of local variables ==================================================!

	! (1) calculate density at nodes ==========================================!	
	do i=1, 100 ! maxnod
		do k = 1,maxlayer
			rho_node(k,i) = 1.03
		end do
	end do

	do i=101, 202 ! maxnod
		do k=1,maxlayer
			rho_node(k,i) = 1.00
		end do
	end do
	
	! (2) calculate density at face ===========================================!
	do j = 1, maxface
		n1 = nodenum_at_face(1,j)
		n2 = nodenum_at_face(2,j)
		do k=1,maxlayer
			rho_face(k,j) = (n1+n2)*0.5
		end do
	end do
	
	! (3) calculate density at cell ===========================================!
	do i=1,maxele
		do k=1,maxlayer
			sum1 = 0.0
			do l=1,tri_or_quad(i)
				n1 = nodenum_at_cell(l,i)
				sum1 = sum1 + rho_node(k,n1)
			end do
			rho_cell(k,i) = sum1/tri_or_quad(i)
		end do
	end do
end subroutine ana_initial_density