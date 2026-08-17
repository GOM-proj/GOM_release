!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This will check if maximum/minimum values satisfiy the given criteria.
!! If it doesn't satisfy, GOM will be terminated.
!!
subroutine model_termination_check
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, j, k
	real(dp):: min_eta, max_eta, max_u, max_v, max_salt, max_temp
	integer :: min_eta_loc(1), max_eta_loc(1)
	integer :: max_u_loc(2), max_v_loc(2), max_salt_loc(2), max_temp_loc(2)
	! jw
		
	! jw
	! jw
	! jw
	! jw
	min_eta = minval(eta_cell)
	max_eta = maxval(eta_cell)
	min_eta_loc = minloc(eta_cell)
	max_eta_loc = maxloc(eta_cell)
	
	! jw
	max_u = maxval(un_face)
	max_v = maxval(vn_face)
	max_u_loc = maxloc(un_face)
	max_v_loc = maxloc(vn_face)
	
	! jw
	max_salt = maxval(salt_cell)
	max_temp = maxval(temp_cell)
	max_salt_loc = maxloc(salt_cell)
	max_temp_loc = maxloc(temp_cell)
	
	! jw
	if(min_eta < eta_min_terminate) then
		i = min_eta_loc(1)
		write(*,*) 				'it = ', it, ', Minimum water surface elevation: ', min_eta, ' < ', eta_min_terminate, ', at cell#: ',i
		write(*,*) 				'Stop'
		write(pw_run_log,*) 	'it = ', it, ', Minimum water surface elevation: ', min_eta, ' < ', eta_min_terminate, ', at cell#: ',i
		write(pw_run_log,*) 	'Stop'
		stop
	else if(max_eta > eta_max_terminate) then
		i = max_eta_loc(1)
		write(*,*) 				'it = ', it, ', Maximum water surface elevation: ', max_eta, ' > ', eta_max_terminate, ', at cell#: ',i
		write(*,*) 				'Stop'
		write(pw_run_log,*) 	'it = ', it, ', Maximum water surface elevation: ', max_eta, ' > ', eta_max_terminate, ', at cell#: ',i
		write(pw_run_log,*) 	'Stop'
		stop
	else if(max_u > uv_terminate) then
		k = max_u_loc(1)
		j = max_u_loc(2)
		write(*,*) 				'it = ', it, ', Maximum horizontal velocity (u): ', max_u, ' > ', uv_terminate, ', at (k,j): ', k, j
		write(*,*) 				'Stop'
		write(pw_run_log,*) 	'it = ', it, ', Maximum horizontal velocity (u): ', max_u, ' > ', uv_terminate, ', at (k,j): ', k, j
		write(pw_run_log,*) 	'Stop'
		stop
	else if(max_v > uv_terminate) then
		j = max_v_loc(1)
		k = max_v_loc(2)
		write(*,*) 				'it = ', it, ', Maximum horizontal velocity (v): ', max_v, ' > ', uv_terminate, ', at (k,j): ', k, j
		write(*,*) 				'Stop'
		write(pw_run_log,*) 	'it = ', it, ', Maximum horizontal velocity (v): ', max_v, ' > ', uv_terminate, ', at (k,j): ', k, j
		write(pw_run_log,*) 	'Stop'
		stop
	else if(max_salt > salt_terminate) then
		i = max_salt_loc(1)
		k = max_salt_loc(2)
		write(*,*) 				'it = ', it, ', Maximum salinity: ', max_salt, ' > ', salt_terminate, ', at (i,k): ', i, k
		write(*,*) 				'Stop'
		write(pw_run_log,*) 	'it = ', it, ', Maximum salinity: ', max_salt, ' > ', salt_terminate, ', at (i,k): ', i, k
		write(pw_run_log,*) 	'Stop'
		stop
	else if(max_temp > temp_terminate) then
		i = max_temp_loc(1)
		k = max_temp_loc(2)
		write(*,*) 				'it = ', it, ', Maximum temperature: ', max_temp, ' > ', temp_terminate, ', at (i,k): ', i, k
		write(*,*) 				'Stop'
		write(pw_run_log,*) 	'it = ', it, ', Maximum temperature: ', max_temp, ' > ', temp_terminate, ', at (i,k): ', i, k
		write(pw_run_log,*) 	'Stop'
		stop
	end if
end subroutine model_termination_check