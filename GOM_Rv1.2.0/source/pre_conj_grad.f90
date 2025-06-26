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
   
   ! jw
   real(dp),dimension(na):: M_inverse, s ! jw
   integer :: na2, nz2
   integer :: colidx2(na),rowstr2(na+1)
   
   real(dp):: sum1
   ! jw
   ! jw
      
	cgit2 = 0
   cgitmax = max_iteration_pcg ! jw
	
	! jw
	! jw
	! jw
	! jw
	! jw
	! jw
	!$omp parallel private(cgit)
	!$omp do private(i,k,sum1)
		do i=1,na
			sum1 = 0.0_dp
			do k=rowstr(i),rowstr(i+1)-1
				if(colidx(k)==i) then
					M_inverse(i) = 1.0/A(k)
				end if
				sum1 = sum1 + A(k)*x(colidx(k)) ! jw
			end do
			colidx2(i) = i
			rowstr2(i) = i
			Ax(i) = sum1 ! jw
			r(i) = b(i) - Ax(i) ! jw
		end do
	!$omp end do nowait
	
	!$omp single
		rowstr2(na+1) = na+1
	!$omp end single	

	! jw
	na2 = na ! jw
	nz2 = na ! jw
	! jw
	!$omp do private(j,k, sum1)
		do j=1,na2
	   	sum1 = 0.0_dp
	      do k=rowstr2(j),rowstr2(j+1)-1
	        	sum1 = sum1 + M_inverse(k)*r(colidx2(k)) ! jw
			end do
			d(j) = sum1 ! jw
		end do
	!$omp end do nowait
	
	
	! jw
	! jw
	!$omp single
		delta_new = 0.0_dp
	!$omp end single
	
	!$omp do private(i) reduction(+:delta_new)
	do i=1,size(r)
	    delta_new = delta_new + r(i) * d(i)
	enddo
	!$omp end do
	
	!$omp single
		delta_0 = delta_new		! jw
	!$omp end single
	!!$omp end parallel

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
	do cgit = 1, cgitmax
		!!$omp parallel
		
		! jw
		! jw
		! jw
		!$omp do private(j,k, sum1)
			do j=1,na
		   	sum1 = 0.0_dp
		      do k=rowstr(j),rowstr(j+1)-1
		        	sum1 = sum1 + A(k)*d(colidx(k)) ! jw
				end do
				q(j) = sum1 ! jw
			end do
		!$omp end do nowait
		
		! jw
		!$omp single
			dq = 0.0_dp
		!$omp end single
		!$omp do private(i) reduction(+:dq)
			do i=1,size(d)
			    dq = dq + d(i) * q(i)
			enddo
		!$omp end do

! jw
! jw
! jw
! jw
! jw
! jw


		! jw
		!$omp single
			alpha = delta_new/dq
		!$omp end single

		!$omp do private(i)
			do i=1,na
				x(i) = x(i) + alpha*d(i)
				r(i) = r(i) - alpha*q(i) 		! jw
         end do
		!$omp end do
		
		
		! jw
		!$omp do private(j,k, sum1)
			do j=1,na2
		   	sum1 = 0.0_dp
		      do k=rowstr2(j),rowstr2(j+1)-1
		        	sum1 = sum1 + M_inverse(k)*r(colidx2(k)) ! jw
				end do
				s(j) = sum1 ! jw
			end do
		!$omp end do nowait
		
		
		! jw
		! jw
		!$omp single
			delta_old = delta_new

			! jw
			delta_new = 0.0_dp
		!$omp end single

		!$omp do private(i) reduction(+:delta_new)
			do i=1,size(r)
			    delta_new = delta_new + r(i) * s(i)
			enddo
		!$omp end do
		
		! jw
		!$omp single
			beta = delta_new / delta_old
		!$omp end single

		! jw
		!$omp do private(i)
			do i=1,na
				d(i) = s(i) + beta*d(i)
			end do
		!$omp end do		
		!!$omp end parallel
 		
 		! jw
 	   if(sqrt(delta_new**2/delta_0**2) < error_tolerance) then ! jw
	   ! jw
	   	! jw
	   	! jw
	   	! jw
 	   	cgit2 = cgit
 	   	exit
	   end if
		
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
 	end do
	!$omp end parallel
		
	if(pcg_result_show == 1) then
		if(mod(it,ishow_frequency)==0) then
			write(*,*) 'Current simulation time step = ', it, ', pre_conj_grad.f90 converged at: ', cgit2
			write(pw_run_log,'(A32,I10,A33,I10)') ' Current simulation time step = ', it, ', pre_conj_grad.f90 converged at:', cgit2
		end if
	end if
end subroutine pre_conj_grad
