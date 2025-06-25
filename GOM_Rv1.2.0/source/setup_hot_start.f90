!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Read some essential variables from restart.inp
!! These variables are required for the momentum and transport equation. 
!! 
subroutine setup_hot_start
   use mod_global_variables
   use mod_file_definition
   implicit none

   integer :: i,j,k
   ! End of local variables ==================================================!
   
   open(pw_restart_inp,file=trim(id_restart_inp),form='unformatted',status='old')

	! eta_cell(maxele), eta_cell_new(maxele)
	do i=1,maxele
		read(pw_restart_inp) eta_cell(i), eta_cell_new(i)
	end do
	
	! eta_node(maxnod)
	do i=1,maxnod
		read(pw_restart_inp) eta_node(i)
	end do
	
	! un_face(maxlayer,maxface), vn_face(maxlayer,maxface)
	! Note, I wrote this as (in write_restart.f90):
	! 		(maxface,maxlayer) format
	do j=1,maxface
		read(pw_restart_inp) (un_face(k,j), k=1,maxlayer), (vn_face(k,j), k=1,maxlayer)
	end do
	
	! u_node(0:maxlayer,maxnod), v_node(0:maxlayer,maxnod)
	do i=1,maxnod
		read(pw_restart_inp) (u_node(k,i), k=0,maxlayer), (v_node(k,i), k=0,maxlayer)
	end do

	! wn_cell(0:maxlayer,maxele)
	! Note, I wrote this as (in write_restart.f90):
	! 		(0:maxlayer,maxele) format
	do i=1,maxele
		read(pw_restart_inp) (wn_cell(k,i), k=0,maxlayer)
	end do
	
	! w_node(0:maxlayer,maxnod)
	do i=1,maxnod
		read(pw_restart_inp) (w_node(k,i), k=0,maxlayer)
	end do
	
	! Transport variables at cell, node, and face =============================!	
	! at cell
	do i=1,maxele
		read(pw_restart_inp) (salt_cell(k,i), k=1,maxlayer), (temp_cell(k,i), k=1,maxlayer), (rho_cell(k,i), k=1,maxlayer)
	end do
	
	! at node
	do i=1,maxnod
		read(pw_restart_inp) (salt_node(k,i), k=1,maxlayer), (temp_node(k,i), k=1,maxlayer), (rho_node(k,i), k=1,maxlayer)
	end do
	
	! at face
	do j=1,maxface
		read(pw_restart_inp) (salt_face(k,j), k=1,maxlayer), (temp_face(k,j), k=1,maxlayer), (rho_face(k,j), k=1,maxlayer)
	end do

	close(pw_restart_inp)
end subroutine setup_hot_start