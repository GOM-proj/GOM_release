!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
! Basic PETSc instruction ====================================================!
!                    Include files
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 
!  This program uses CPP for preprocessing, as indicated by the use of
!  PETSc include files in the directory petsc/include/finclude.  This
!  convention enables use of the CPP preprocessor, which allows the use
!  of the #include statements that define PETSc objects and variables.
! 
!  Use of the conventional Fortran include statements is also supported
!  In this case, the PETsc include files are located in the directory
!  petsc/include/foldinclude.
! 
!  Since one must be very careful to include each file no more than once
!  in a Fortran routine, application programmers must exlicitly list
!  each file needed for the various PETSc components within their
!  program (unlike the C/C++ interface).
! 
!  See the Fortran section of the PETSc users manual for details.
! 
!  The following include statements are required for KSP Fortran programs:
!     petscsys.h    - base PETSc routines
!     petscvec.h    - vectors
!     petscmat.h    - matrices
!     petscksp.h    - Krylov subspace methods
!     petscpc.h     - preconditioners
!  Other include statements may be needed if using additional PETSc
!  routines in a Fortran program, e.g.,
!     petscviewer.h - viewers
!     petscis.h     - index sets
! 
! 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                   Variable declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 
!  Variables:
!     ksp      - linear solver context
!     ksp      - Krylov subspace method context
!     pc       - preconditioner context
!     x, b, u  - approx solution, right-hand-side, exact solution vectors
!     A        - matrix that defines linear system
!     its      - iterations for convergence
!     norm     - norm of error in solution
! 
! End of PETSc instruction ===================================================!

!=============================================================================!
! jw's note:
! In this example, I will:
! 		(1) create a sparse matrix
! 		(2) insert/add values in the matrix
! 		(3) solve Ax = b with 
! 		(4) print the outcoms.
!     (5) destroy the matrix
! all by myself.
! Adding more to ex5.F90
! This might be good enough for my purpose in GOM...
!=============================================================================!

subroutine jw_petsc(sparsem,rhs2)
! note that this include statement should be included and without any space before #
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscviewer.h>
	use petsc
	implicit none

	! declare PETSc variables -------------------------------------------------!
	Mat					:: A 		! 3 by 3 matrix
	Vec					:: x, b 	! x = approx solution, b = RHS (known vector)
	KSP					:: ksp	! Krylov subspace method context
	PC						:: pc		! preconditioner context
	MPI_Comm 			:: comm
 	PetscErrorCode 	:: ierr
	PetscInt				:: MPI_rank, MPI_size ! integer for Petsc & MPI
	PetscInt				:: m, n 	! global m by n matrix
	
	PetscInt,allocatable 	:: row(:), col(:) ! indexes of row and column
	PetscScalar,allocatable	:: value(:)  ! values to be inserted.
	PetscViewer 		:: view_out
	PetscReal			:: rtol
	KSPConvergedReason:: converge_reason 
	PetscInt				:: its 	! total iterations
	! End of petsc variables --------------------------------------------------!
	
	! local variables ---------------------------------------------------------!
	real(dp),intent(in) :: sparsem, rhs2
	integer :: i, j
	character(len=100) :: textbuff
	! PetscScalar, pointer :: xx_v(:) ! to read result vector
	! End of local variables ==================================================!


	! Required field to use Petsc =============================================!
	! The general form of a PETSc program includes the following structure:
	! It should start calling 'PetscInitialize()':
	! 		call PetscInitialize(character(*) file, integer ierr)
	! PetscInitialize() calls MPI_Init() if MPI has not been previously initialized.
	! If MPI needs to be initialized directly (or is initialized by some other library), MPI_Init() can be called before PetscInitialize().
	! By default, PetscInitialize() sets the PETSc world communitor PETSC_COMM_WORLD to MPI_COMM_WORLD.
	! In Fortran, you MUST use PETSC_NULL_CHARACTER to indicate a null character string. 
	
	call PetscInitialize(PETSC_NULL_CHARACTER,ierr);  CHKERRA(ierr) ! It is a good idea to check this error code after every call to a PETSc routine!
	comm = PETSC_COMM_WORLD
	
	! As I mentioned above, PetscInitialize() already sets PETSC_COMM_WORLD to MPI_COMM_WORLD.
	! Thus, I can use MPI_COMM_WORLD directly in MPI functions/subroutines.
	! These are for MPI setting.
	! I will put this later...
	! call MPI_Comm_rank(MPI_COMM_WORLD,MPI_rank,ierr);  CHKERRA(ierr) 	! determin the rank of the calling process in the communicator
	call MPI_Comm_size(MPI_COMM_WORLD,MPI_size,ierr);  CHKERRA(ierr)		! returns the size of the group associated with a communicator

	if (MPI_size /= 1) then
		call MPI_Comm_rank(PETSC_COMM_WORLD,MPI_rank,ierr)
		if (MPI_rank == 0) then
			write(*,*) 'This is a uniprocessor example only!'
		end if
      ! SETERRQ(PETSC_COMM_WORLD,1,' ',ierr) ! jw,
      write(*,'(A,I5,A)') 'We are using ', MPI_size, ' nodes.'
	end if
	! End of required field ===================================================!
	
	! Now, we can use Petsc ===================================================!	
	! (1) create matrix A:
	! 		A = [5 -2  0
	!         -2  5  1
	!          0  1  5]
	! The answer should be:
	! 		b = [6 5 -3]^T (T: transpose)
	! (1.1) Create an empty (no size yet) matrix A, actually A is a pointer on the matrix not a real matrix as in Fortran. 
	! By default, MatCreate use sparse AIJ format
	m = maxele
	n = maxele
	call MatCreate(PETSC_COMM_WORLD,A,ierr)	
	
	! (1.2) set matrix A = [m by n] matrix
	call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n,ierr) 
	
	! (1.3) The default matrix type is AIJ. So, actually it is not required if you just want to use AIJ type matrix.
	! But, this may enforce the matrix format before runtime.
	call MatSetFromOptions(A,ierr)
	call MatSetUp(A,ierr) ! I don't know wheter I really need this or not.
		
	! now, let's include values in A:
	! MatSetValues() inserts a [m by n] block of values in the matrix.
	! This is kind of a weird process.
	! Tha matrix index should starts from "0", but the row() and col() should starts from "1".

   ! construct A matrix
   do i=1,maxele
   	do j=1,maxele
   		row(i) = i-1 ! row index: [0:maxele-1]
   		col(j) = j-1 ! coloumn index: [0:maxele-1]
   		if(i == j) then
   			value(j) = sparsem(i,0) ! diagonal value
   		else
   			if()
   			value(j) = sparsem(i,ii) ! something like this
 			call MatSetValues(A,1,row,maxele,col,value,INSERT_VALUES,ierr)
   		
   	end do
   end do




	! 1st row
	do i=1,3
		row(i) = 0 ! row index
		col(i) = i-1 ! column index
	end do
	value(1) = 5
	value(2) = -2
	value(3) = 0
 	call MatSetValues(A,1,row,3,col,value,INSERT_VALUES,ierr)
	
	! 2nd row
	do i=1,3
		row(i) = 1 ! row index
		col(i) = i-1 ! column index
	end do
	value(1) = -2
	value(2) = 5
	value(3) = 1
 	call MatSetValues(A,1,row,3,col,value,INSERT_VALUES,ierr)
	
	! 3rd row
	do i=1,3
		row(i) = 2 ! row index
		col(i) = i-1 ! column index
	end do
	value(1) = 0
	value(2) = 1
	value(3) = 5	
 	call MatSetValues(A,1,row,3,col,value,INSERT_VALUES,ierr)
	
	! After all elements have been inserted into the matrix,
	! it must be processed with the pair of following commands.
	! This is necessary because values may be buffered; so this pair of calls
	! ensures all relevant buffers are flushed.
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
	
	! check the matrix:
	! MatView() requires <petsc/finclude/petscviewer.h>
	call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr) ! display on the standard output device
	
	! this is another way:
	view_out = PETSC_VIEWER_STDOUT_WORLD
	call MatView(A,view_out,ierr) ! display on the standard output device
	
	! print out to a file:
	call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./mat.output",view_out,ierr) ! open an ASCII file
	call MatView(A,view_out,ierr) ! then, display on the output file
	! End of matrix A =========================================================!
	
	! Now, lets create vectors x & b: =========================================!
	call VecCreate(PETSC_COMM_WORLD,x,ierr) ! jw, create empty vector x (no size yet)
	call VecSetSizes(x,PETSC_DECIDE,n,ierr) ! jw, define x vector to be [n by 1] vector for every MPI nodes???
	call VecSetFromOptions(x,ierr)
	call VecDuplicate(x,b,ierr) ! jw, duplicate x vector to b

	! feed vector b, which is the known values
	! There are two methods: 
	! (1) using VecSet() - this will feed a vector with just one identical values
	! (2) using VecSetValues() - this will feed ad vector with different values
	do i=1,3
		row(i) = i-1 ! row index (actually, it is not a row index but just an index since this is not a matrix but a vector)
	end do
	value(1) = 20
	value(2) = 10
	value(3) = -10
	call VecSetValues(b,3,row,value,INSERT_VALUES,ierr) ! 3 is the length of the vector
	
	! After all elements have been inserted into the vector, you have to do this:
	call VecAssemblyBegin(b,ierr)
	call VecAssemblyEnd(b,ierr)
	
	! check the matrix:
	! VecView() requires <petsc/finclude/petscvec.h>
	call VecView(b,view_out,ierr) ! display on the output file: mat.output since view_out is alrealy defined as this file previously
	! End of vector creation ==================================================!
	
	! Solve linear/nonlinear system ===========================================!
	! Create linear solver context
	! Note the following two subroutines work as a pair.
	call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
	! Before actually solving a linear system with KSP, the user must call the following routine:
	call KSPSetOperators(ksp,A,A,ierr) ! normally this setup would work

	! Set runtime options, e.g.,
	! -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
	! These options will override those specified above as long as
	! KSPSetFromOptions() is called _after_ any other customization routines.
	call KSPSetFromOptions(ksp,ierr)

	! Set linear solver defaults for this problem (optional).
	! - By extracting the KSP and PC contexts from the KSP context,
	!   we can then directly call any KSP and PC routines to set various options.
	! - The following four statements are optional; all of these
	!   parameters could alternatively be specified at runtime via
	!   KSPSetFromOptions();
	! The command KSPSetFromOptions() enables customization of the linear solution method at runtime
	! via the options database, e.g., choice of iterative method, preconditioner, convergence tolerance, etc.
	
	! Default solver is the block Jacobi method with ILU(0): with one block per process, each of which is solved with ILU(0).
	! However, we can retrieve ksp and pc context and then select those method with other ones.
	
	! first, let's change the solver
	call KSPSetType(ksp,KSPCG,ierr) ! let's change the default GMRES to CG
	
	! second, let's change the preconditioner
	call KSPGetPC(ksp,pc,ierr) ! jw, ksp: input, pc: output; this will return 'pc', which is the preconditioner context, now we can choose other pc
	! Select PC for a particular preconditioner type rather than just using the default pc (i.e., ILU)
	! PCSetType(PC pc, PCType type)
	call PCSetType(pc,PCJACOBI,ierr) ! PCJACOBI is the preconditioner type we want to use
	
	! set the convergence criteria: e.g., tolerance and maximum iteration
	! There are three types of tolerances:
	! 		rtol: the residual norm relative to the norm of the right hand side
	! 				default value = 1e-5
	! 		atol: the absolute size of the residual norm
	! 				default value = 1e-50
	! 		dtol: the relative increase in the residual
	! 				default value = 1e5
	! 	The maximum iteration (maxits):
	! 				default value = 1e4
	! See more details in the manual: 4.3.2
	! This function has the following form:
	! 		KSPSetTolerances(ksp,rtol,atol,dtol,maxits)
	! Let's just reduce the rtol here, so other values are same as the default values
	rtol = 1.d-7
	call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
	
	! Now, you are ready to solve the linear system ---------------------------!
	! Solve the linear system
	call KSPSolve(ksp,b,x,ierr) ! b = right hand side vector, x = solution
	
	! analyze the solving process
	call KSPGetConvergedReason(ksp,converge_reason,ierr)
	if(converge_reason < 0) then
		write(textbuff,'(A,I3,A)') 'Failure to converge', converge_reason,'.\n'
		call PetscPrintf(PETSC_COMM_WORLD,textbuff,ierr) ! print on the screen
		call PetscViewerASCIIPrintf(view_out,textbuff,ierr) ! print on the file
	else
		call KSPGetIterationNumber(ksp,its,ierr)
		write(textbuff,'(A,I5,A)') 'Converged in', its, ' iterations.\n'
		call PetscPrintf(PETSC_COMM_WORLD,textbuff,ierr) ! print on the screen
		call PetscViewerASCIIPrintf(view_out,textbuff,ierr) ! print on the file
	end if
	
	! check the results
	call VecView(x,view_out,ierr) ! display on the output file: mat.output since view_out is alrealy defined as this file previously
	
	 	
	! View solver info; we could instead use the option -ksp_view	
	! There are two options:
	!  	PETSC_VIEWER_STDOUT_SELF; 	standard output (default)
	!  	PETSC_VIEWER_STDOUT_WORLD; synchronized standard output where only the first processor opens the file
	!                                All other porcessors send their data to the first processor to print.
	call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)	
	! End of solving ==========================================================!
	
	
	! Destroy the matrix ======================================================!	
	! Free work space. All PETSc objects should be destroyed when they are no longer needed
	call VecDestroy(x,ierr)
	call VecDestroy(b,ierr)
	call MatDestroy(A,ierr)
	call KSPDestroy(ksp,ierr)
	call PetscFinalize(ierr)
	
end subroutine jw_petsc