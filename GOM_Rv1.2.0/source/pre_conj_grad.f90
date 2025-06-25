!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This is the Jacobi preconditioned conjugate gradient method (jcg)
!! This version is based on the original pre_conj_grad_v0.f90
!! This version is much faster than the original version since I removed (reduced) omp setup overhead as much as possible.
!! 
subroutine pre_conj_grad(na,nz,A,colidx,rowstr,b,x)
	use mod_global_variables, only : dp, it, max_iteration_pcg, error_tolerance, pcg_result_show, ishow_frequency
	use mod_file_definition, only : pw_run_log
	use omp_lib
   implicit none
	
	integer,intent(in) :: na, nz
	real(dp),dimension(nz),intent(in) :: A
	integer, dimension(nz),intent(in) :: colidx
	integer, dimension(na+1),intent(in) :: rowstr
	real(dp),dimension(na),intent(in) :: b
	real(dp),dimension(na),intent(inout) :: x
	
   integer :: i,j,k
   integer :: cgit, cgitmax, cgit2
   real(dp):: dq, delta_new, delta_old, delta_0, alpha, beta
   
   real(dp),dimension(na):: q,r,d,Ax
   
   ! for precondioner matrix, M
   real(dp),dimension(na):: M_inverse, s ! M is the preconditioner
   integer :: na2, nz2
   integer :: colidx2(na),rowstr2(na+1)
   
   real(dp):: sum1
   ! integer :: num_threads, myrank
   ! End of local variables ==================================================!
      
	cgit2 = 0
   cgitmax = max_iteration_pcg ! maximum iteration   
	
	! Find the Jacobi preconditioner, M, and its inverse matrix, M_inverse.
	! Jacobi preconditioner, M, is the diagonal matrix of the original matrix.
	! And, the inverse matrix of M is the 1/M(i)
	! And, Calculate: A.x
	! Ax = multmv(na,nz,colidx,rowstr,A,x) ! Ax	
	! r = b-Ax
	!$omp parallel private(cgit)
	!$omp do private(i,k,sum1)
		do i=1,na
			sum1 = 0.0_dp
			do k=rowstr(i),rowstr(i+1)-1
				if(colidx(k)==i) then
					M_inverse(i) = 1.0/A(k)
				end if
				sum1 = sum1 + A(k)*x(colidx(k)) ! row-wise matrix multiplication
			end do
			colidx2(i) = i
			rowstr2(i) = i
			Ax(i) = sum1 ! = A.v
			r(i) = b(i) - Ax(i) ! r = b-Ax
		end do
	!$omp end do nowait
	
	!$omp single
		rowstr2(na+1) = na+1
	!$omp end single	

	! d = M_inverse.r
	na2 = na ! matrix size
	nz2 = na ! total indices in the diagonal matrix
	! d = multmv(na2,nz2,colidx2,rowstr2,M_inverse,r)
	!$omp do private(j,k, sum1)
		do j=1,na2
	   	sum1 = 0.0_dp
	      do k=rowstr2(j),rowstr2(j+1)-1
	        	sum1 = sum1 + M_inverse(k)*r(colidx2(k)) ! row-wise matrix multiplication
			end do
			d(j) = sum1 ! = A.v
		end do
	!$omp end do nowait
	
	
	! delta_new = dotprod(r,d) ! delta_new = r^T.d
	! initialize for reduction
	!$omp single
		delta_new = 0.0_dp
	!$omp end single
	
	!$omp do private(i) reduction(+:delta_new)
	do i=1,size(r)
	    delta_new = delta_new + r(i) * d(i)
	enddo
	!$omp end do
	
	!$omp single
		delta_0 = delta_new		! size of residual norm
	!$omp end single
	!!$omp end parallel

! 	write(*,*) 'delta_new, delta_0 = ', delta_new, delta_0
! 	write(*,*) 'r, d = :'
! 	do i=1,na
! 		write(*,*) r(i), d(i)
! 	end do
! 	stop 'jjj'

	
	! Note: I tried keep omp for the following do loop also.
	! It works, but the final 'cgit' printing values was wrong even though the results are ok. 
	! So, I set !$omp paralle in each do loop
	! Now, it works correctly with including cgit2. 
	do cgit = 1, cgitmax
		!!$omp parallel
		
		! write(*,*) cgit 
		! q = multmv(A,d) 						! q = A.d, matrix & vector multiplication
		! q = multmv(na,nz,colidx,rowstr,A,d) ! q = A.d, matrix & vector multiplication
		!$omp do private(j,k, sum1)
			do j=1,na
		   	sum1 = 0.0_dp
		      do k=rowstr(j),rowstr(j+1)-1
		        	sum1 = sum1 + A(k)*d(colidx(k)) ! row-wise matrix multiplication
				end do
				q(j) = sum1 ! = A.d
			end do
		!$omp end do nowait
		
		! dq =  dotprod(d,q) 						! d^T.q
		!$omp single
			dq = 0.0_dp
		!$omp end single
		!$omp do private(i) reduction(+:dq)
			do i=1,size(d)
			    dq = dq + d(i) * q(i)
			enddo
		!$omp end do

! 		write(*,*) 'q:'
! 		do i=1,na
! 			write(*,*) q(i)
! 		end do
! 		write(*,*) 'dq = ', dq
! 		stop 'jjj'


		! alpha = delta_new/dq   					! delta_new/d^T.q		
		!$omp single
			alpha = delta_new/dq
		!$omp end single

		!$omp do private(i)
			do i=1,na
				x(i) = x(i) + alpha*d(i)
				r(i) = r(i) - alpha*q(i) 		! if i is in-divisible by 50		            
         end do
		!$omp end do
		
		
		! s = multmv(na2,nz2,colidx2,rowstr2,M_inverse,r)
		!$omp do private(j,k, sum1)
			do j=1,na2
		   	sum1 = 0.0_dp
		      do k=rowstr2(j),rowstr2(j+1)-1
		        	sum1 = sum1 + M_inverse(k)*r(colidx2(k)) ! row-wise matrix multiplication
				end do
				s(j) = sum1 ! = A.v
			end do
		!$omp end do nowait
		
		
		! delta_old = delta_new
		! delta_new = dotprod(r,s) 				! delta_new = r^T.s
		!$omp single
			delta_old = delta_new

			! initialize for reduction
			delta_new = 0.0_dp
		!$omp end single

		!$omp do private(i) reduction(+:delta_new)
			do i=1,size(r)
			    delta_new = delta_new + r(i) * s(i)
			enddo
		!$omp end do
		
		! beta = delta_new / delta_old			! beta = delta_new/delta_old
		!$omp single
			beta = delta_new / delta_old
		!$omp end single

		! d = s + beta*d
		!$omp do private(i)
			do i=1,na
				d(i) = s(i) + beta*d(i)
			end do
		!$omp end do		
		!!$omp end parallel
 		
 		! This is the first approach, and it is faster than the second approach
 	   if(sqrt(delta_new**2/delta_0**2) < error_tolerance) then ! This is the original CG method
	   ! if(sqrt(delta_new/na) < error_tolerance) then ! kinds of RMSE is less than 1 mm
	   	! All threads will be reached to identical cgit,
	   	! and all threads will update cgit2 (global value) to cgit (its own private value)
	   	! thus, the final cgit2 will have correct cgit value.
 	   	cgit2 = cgit
 	   	exit
	   end if
		
		! This is the second approach, and I think this is more stable approach.
		! However, it is slower than the first approach.
! 		!$omp single
! 			if(sqrt(delta_new**2/delta_0**2) < error_tolerance) then ! This is the original CG method
! 				cgit2 = cgit
! 			end if
! 		!$omp end single
! 		
! 		if(cgit2 > 0) then
! 			if(myrank < num_threads) then
! 				exit
! 			end if
! 		end if
 	end do
	!$omp end parallel
		
	if(pcg_result_show == 1) then
		if(mod(it,ishow_frequency)==0) then
			write(*,*) 'Current simulation time step = ', it, ', pre_conj_grad.f90 converged at: ', cgit2
			write(pw_run_log,'(A32,I10,A33,I10)') ' Current simulation time step = ', it, ', pre_conj_grad.f90 converged at:', cgit2
		end if
	end if
end subroutine pre_conj_grad
