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
   ! jw
   
   open(pw_restart_inp,file=trim(id_restart_inp),form='unformatted',status='old')

	! jw
	do i=1,maxele
		read(pw_restart_inp) eta_cell(i), eta_cell_new(i)
	end do
	
	! jw
	do i=1,maxnod
		read(pw_restart_inp) eta_node(i)
	end do
	
	! jw
	! jw
	! jw
	do j=1,maxface
		read(pw_restart_inp) (un_face(k,j), k=1,maxlayer), (vn_face(k,j), k=1,maxlayer)
	end do
	
	! jw
	do i=1,maxnod
		read(pw_restart_inp) (u_node(k,i), k=0,maxlayer), (v_node(k,i), k=0,maxlayer)
	end do

	! jw
	! jw
	! jw
	do i=1,maxele
		read(pw_restart_inp) (wn_cell(k,i), k=0,maxlayer)
	end do
	
	! jw
	do i=1,maxnod
		read(pw_restart_inp) (w_node(k,i), k=0,maxlayer)
	end do
	
	! jw
	! jw
	do i=1,maxele
		read(pw_restart_inp) (salt_cell(k,i), k=1,maxlayer), (temp_cell(k,i), k=1,maxlayer), (rho_cell(k,i), k=1,maxlayer)
	end do
	
	! jw
	do i=1,maxnod
		read(pw_restart_inp) (salt_node(k,i), k=1,maxlayer), (temp_node(k,i), k=1,maxlayer), (rho_node(k,i), k=1,maxlayer)
	end do
	
	! jw
	do j=1,maxface
		read(pw_restart_inp) (salt_face(k,j), k=1,maxlayer), (temp_face(k,j), k=1,maxlayer), (rho_face(k,j), k=1,maxlayer)
	end do

	close(pw_restart_inp)
end subroutine setup_hot_start