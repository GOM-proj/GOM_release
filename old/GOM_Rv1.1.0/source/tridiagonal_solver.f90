!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! Tridiagonal matrix solver by thomas algorithm:
!!		nmax  		: dimension of a etc in the driving routine.     
!!		n     		: actual rank of the system. actual number of vertical layers (= active water column layers at jth side)         
!!		nc    		: input; actual # of columns of rhs. We are using 3 in solve_momentum_equation.f90.
!!		mat_a 		: lower matrix
!!		mab_b 		: diagonal matrix
!!		mat_c 		: upper matrix
!!		rrhs  		: right hand side                                
!!		solution    : output with nc columns (depending on input nc).
!!		gam   		: a working array.                               
! jw
subroutine tridiagonal_solver(nmax,n,nc,mat_a,mat_b,mat_c,rhs,solution,gam)
	use mod_global_variables, only : dp
   use mod_file_definition

   implicit none   

   integer,intent(in) :: nmax, n, nc
   real(dp),dimension(nmax)  ,intent(in)  :: mat_a, mat_b, mat_c
   real(dp),dimension(nmax,nc),intent(in)  :: rhs 						! jw
   real(dp),dimension(nmax,nc),intent(out) :: solution 				! jw
   real(dp),dimension(nmax)  ,intent(out) :: gam

   integer :: i, j
   real(dp) :: denominator 
   ! jw
	
	! jw
   if(mat_b(1) == 0.0) then
      write(pw_run_log,*) 'Diagonal matrix mat_b(1) = 0.0, should not be zero'
      stop
   end if
	
	! jw
	! jw
   denominator = mat_b(1)
   do i = 1, nc
      solution(1,i) = rhs(1,i) / denominator
   end do

   do j = 2, n
      gam(j) = mat_c(j-1) / denominator
      denominator = mat_b(j) - mat_a(j)*gam(j)
      if(denominator == 0.0_dp) then
         write(pw_run_log,*) 'Tridiagonal solver failed'
         stop
      end if
      do i = 1, nc
         solution(j,i) = (rhs(j,i)-mat_a(j)*solution(j-1,i)) / denominator
      end do !i
   end do !j

	! jw
   do j = n-1, 1, -1
   	do i = 1, nc
        solution(j,i) = solution(j,i)-gam(j+1)*solution(j+1,i)
     end do
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

end subroutine tridiagonal_solver
