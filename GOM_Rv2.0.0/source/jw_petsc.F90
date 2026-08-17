!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
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

!=============================================================================!
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
!=============================================================================!

subroutine jw_petsc(sparsem,rhs2)
! jw
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscviewer.h>
	use petsc
	implicit none

	! jw
	Mat					:: A 		! jw
	Vec					:: x, b 	! jw
	KSP					:: ksp	! jw
	PC						:: pc		! jw
	MPI_Comm 			:: comm
 	PetscErrorCode 	:: ierr
	PetscInt				:: MPI_rank, MPI_size ! jw
	PetscInt				:: m, n 	! jw
	
	PetscInt,allocatable 	:: row(:), col(:) ! jw
	PetscScalar,allocatable	:: value(:)  ! jw
	PetscViewer 		:: view_out
	PetscReal			:: rtol
	KSPConvergedReason:: converge_reason 
	PetscInt				:: its 	! jw
	! jw
	
	! jw
	real(dp),intent(in) :: sparsem, rhs2
	integer :: i, j
	character(len=100) :: textbuff
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
	
	call PetscInitialize(PETSC_NULL_CHARACTER,ierr);  CHKERRA(ierr) ! jw
	comm = PETSC_COMM_WORLD
	
	! jw
	! jw
	! jw
	! jw
	! jw
	call MPI_Comm_size(MPI_COMM_WORLD,MPI_size,ierr);  CHKERRA(ierr)		! jw

	if (MPI_size /= 1) then
		call MPI_Comm_rank(PETSC_COMM_WORLD,MPI_rank,ierr)
		if (MPI_rank == 0) then
			write(*,*) 'This is a uniprocessor example only!'
		end if
      ! jw
      write(*,'(A,I5,A)') 'We are using ', MPI_size, ' nodes.'
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
	m = maxele
	n = maxele
	call MatCreate(PETSC_COMM_WORLD,A,ierr)	
	
	! jw
	call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n,ierr) 
	
	! jw
	! jw
	call MatSetFromOptions(A,ierr)
	call MatSetUp(A,ierr) ! jw
		
	! jw
	! jw
	! jw
	! jw

   ! jw
   do i=1,maxele
   	do j=1,maxele
   		row(i) = i-1 ! jw
   		col(j) = j-1 ! jw
   		if(i == j) then
   			value(j) = sparsem(i,0) ! jw
   		else
   			if()
   			value(j) = sparsem(i,ii) ! jw
 			call MatSetValues(A,1,row,maxele,col,value,INSERT_VALUES,ierr)
   		
   	end do
   end do




	! jw
	do i=1,3
		row(i) = 0 ! jw
		col(i) = i-1 ! jw
	end do
	value(1) = 5
	value(2) = -2
	value(3) = 0
 	call MatSetValues(A,1,row,3,col,value,INSERT_VALUES,ierr)
	
	! jw
	do i=1,3
		row(i) = 1 ! jw
		col(i) = i-1 ! jw
	end do
	value(1) = -2
	value(2) = 5
	value(3) = 1
 	call MatSetValues(A,1,row,3,col,value,INSERT_VALUES,ierr)
	
	! jw
	do i=1,3
		row(i) = 2 ! jw
		col(i) = i-1 ! jw
	end do
	value(1) = 0
	value(2) = 1
	value(3) = 5	
 	call MatSetValues(A,1,row,3,col,value,INSERT_VALUES,ierr)
	
	! jw
	! jw
	! jw
	! jw
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
	
	! jw
	! jw
	call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr) ! jw
	
	! jw
	view_out = PETSC_VIEWER_STDOUT_WORLD
	call MatView(A,view_out,ierr) ! jw
	
	! jw
	call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./mat.output",view_out,ierr) ! jw
	call MatView(A,view_out,ierr) ! jw
	! jw
	
	! jw
	call VecCreate(PETSC_COMM_WORLD,x,ierr) ! jw
	call VecSetSizes(x,PETSC_DECIDE,n,ierr) ! jw
	call VecSetFromOptions(x,ierr)
	call VecDuplicate(x,b,ierr) ! jw

	! jw
	! jw
	! jw
	! jw
	do i=1,3
		row(i) = i-1 ! jw
	end do
	value(1) = 20
	value(2) = 10
	value(3) = -10
	call VecSetValues(b,3,row,value,INSERT_VALUES,ierr) ! jw
	
	! jw
	call VecAssemblyBegin(b,ierr)
	call VecAssemblyEnd(b,ierr)
	
	! jw
	! jw
	call VecView(b,view_out,ierr) ! jw
	! jw
	
	! jw
	! jw
	! jw
	call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
	! jw
	call KSPSetOperators(ksp,A,A,ierr) ! jw

	! jw
	! jw
	! jw
	! jw
	call KSPSetFromOptions(ksp,ierr)

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
	call KSPSetType(ksp,KSPCG,ierr) ! jw
	
	! jw
	call KSPGetPC(ksp,pc,ierr) ! jw
	! jw
	! jw
	call PCSetType(pc,PCJACOBI,ierr) ! jw
	
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
	rtol = 1.d-7
	call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
	
	! jw
	! jw
	call KSPSolve(ksp,b,x,ierr) ! jw
	
	! jw
	call KSPGetConvergedReason(ksp,converge_reason,ierr)
	if(converge_reason < 0) then
		write(textbuff,'(A,I3,A)') 'Failure to converge', converge_reason,'.\n'
		call PetscPrintf(PETSC_COMM_WORLD,textbuff,ierr) ! jw
		call PetscViewerASCIIPrintf(view_out,textbuff,ierr) ! jw
	else
		call KSPGetIterationNumber(ksp,its,ierr)
		write(textbuff,'(A,I5,A)') 'Converged in', its, ' iterations.\n'
		call PetscPrintf(PETSC_COMM_WORLD,textbuff,ierr) ! jw
		call PetscViewerASCIIPrintf(view_out,textbuff,ierr) ! jw
	end if
	
	! jw
	call VecView(x,view_out,ierr) ! jw
	
	 	
	! jw
	! jw
	! jw
	! jw
	! jw
	call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)	
	! jw
	
	
	! jw
	! jw
	call VecDestroy(x,ierr)
	call VecDestroy(b,ierr)
	call MatDestroy(A,ierr)
	call KSPDestroy(ksp,ierr)
	call PetscFinalize(ierr)
	
end subroutine jw_petsc