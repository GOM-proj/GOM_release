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
   ! jw
   
   format_string = '(I7.7)'	! jw
	write(File_num_buff,format_string) it
	restart_file_name = trim(id_restart_out)//trim(File_num_buff)//'.out'
		
	open(pw_restart_out,file=trim(restart_file_name),form='unformatted',status='replace')
	
	! jw
	do i=1,maxele
		write(pw_restart_out) eta_cell(i), eta_cell_new(i)
	end do
	
	! jw
	do i=1,maxnod
		write(pw_restart_out) eta_node(i)
	end do
	
	! jw
	! jw
	! jw
	do j=1,maxface
		write(pw_restart_out) (un_face(k,j), k=1,maxlayer), (vn_face(k,j), k=1,maxlayer)
	end do
	
	! jw
	do i=1,maxnod
		write(pw_restart_out) (u_node(k,i), k=0,maxlayer), (v_node(k,i), k=0,maxlayer)
	end do

	! jw
	! jw
	! jw
	do i=1,maxele
		write(pw_restart_out) (wn_cell(k,i), k=0,maxlayer)
	end do
	
	! jw
	do i=1,maxnod
		write(pw_restart_out) (w_node(k,i), k=0,maxlayer)
	end do
	
	! jw
	! jw
	do i=1,maxele
		write(pw_restart_out) (salt_cell(k,i), k=1,maxlayer), (temp_cell(k,i), k=1,maxlayer), (rho_cell(k,i), k=1,maxlayer)
	end do
	
	! jw
	do i=1,maxnod
		write(pw_restart_out) (salt_node(k,i), k=1,maxlayer), (temp_node(k,i), k=1,maxlayer), (rho_node(k,i), k=1,maxlayer)
	end do
	
	! jw
	do j=1,maxface
		write(pw_restart_out) (salt_face(k,j), k=1,maxlayer), (temp_face(k,j), k=1,maxlayer), (rho_face(k,j), k=1,maxlayer)
	end do
		
	close(pw_restart_out)
end subroutine write_restart