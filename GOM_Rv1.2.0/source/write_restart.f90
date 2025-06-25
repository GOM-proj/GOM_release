!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! Wrtie some essential variables, which are required for momentum and transport equations.
!! 
subroutine write_restart
   use mod_global_variables 
   use mod_file_definition
   implicit none

   integer :: i, j, k
   character(len=100):: restart_file_name
	character(len=40) :: format_string
	character(len= 7) :: File_num_buff   
   ! End of local variables ==================================================!
   
   format_string = '(I7.7)'	! an integer of width 7 with zeros at the left. For example: 0000001 ~ 9999999
	write(File_num_buff,format_string) it
	restart_file_name = trim(id_restart_out)//trim(File_num_buff)//'.out'
		
	open(pw_restart_out,file=trim(restart_file_name),form='unformatted',status='replace')
	
	! eta_cell(maxele), eta_cell_new(maxele)
	do i=1,maxele
		write(pw_restart_out) eta_cell(i), eta_cell_new(i)
	end do
	
	! eta_node(maxnod)
	do i=1,maxnod
		write(pw_restart_out) eta_node(i)
	end do
	
	! un_face(maxlayer,maxface), vn_face(maxlayer,maxface)
	! Note, I will write this as:
	! 		(maxface,maxlayer) format
	do j=1,maxface
		write(pw_restart_out) (un_face(k,j), k=1,maxlayer), (vn_face(k,j), k=1,maxlayer)
	end do
	
	! u_node(0:maxlayer,maxnod), v_node(0:maxlayer,maxnod)
	do i=1,maxnod
		write(pw_restart_out) (u_node(k,i), k=0,maxlayer), (v_node(k,i), k=0,maxlayer)
	end do

	! wn_cell(0:maxlayer,maxele)
	! Note, I will write this as:
	! 		(maxele,0:maxlayer) format
	do i=1,maxele
		write(pw_restart_out) (wn_cell(k,i), k=0,maxlayer)
	end do
	
	! w_node(0:maxlayer,maxnod)
	do i=1,maxnod
		write(pw_restart_out) (w_node(k,i), k=0,maxlayer)
	end do
	
	! Transport variables at cell, node, and face =============================!
	! at cell
	do i=1,maxele
		write(pw_restart_out) (salt_cell(k,i), k=1,maxlayer), (temp_cell(k,i), k=1,maxlayer), (rho_cell(k,i), k=1,maxlayer)
	end do
	
	! at node
	do i=1,maxnod
		write(pw_restart_out) (salt_node(k,i), k=1,maxlayer), (temp_node(k,i), k=1,maxlayer), (rho_node(k,i), k=1,maxlayer)
	end do
	
	! at face
	do j=1,maxface
		write(pw_restart_out) (salt_face(k,j), k=1,maxlayer), (temp_face(k,j), k=1,maxlayer), (rho_face(k,j), k=1,maxlayer)
	end do
		
	close(pw_restart_out)
end subroutine write_restart