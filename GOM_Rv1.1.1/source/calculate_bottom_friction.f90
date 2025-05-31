!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine calculates bottom drag coefficient:
!! 		Cdb: depending on depth
!! 		Gamma_B: depending on east/west velocities
!! 
subroutine calculate_bottom_friction
	use mod_global_variables
	use mod_file_definition
! jw
   implicit none
   
   integer :: j, k
   real(dp):: chezy_coef, hydraulic_radius, expo
   real(dp):: rtemp
   ! jw

	! jw
	Cdb 		= 0.0_dp	! jw
	Gamma_B 	= 0.0_dp ! jw
	
	if(bf_flag == 0) then	! jw
		! jw
	else if(bf_flag == 1) then	! jw
		expo = 1.0_dp/6.0_dp ! jw
		!$omp parallel do private(j,k,hydraulic_radius,chezy_coef,rtemp)
		do j=1,maxface
         if(h_face(j) > dry_depth) then
            ! jw
            hydraulic_radius = max(1.0,h_face(j)) ! jw
            chezy_coef = (hydraulic_radius**expo) / manning(j) ! jw
            rtemp = gravity / chezy_coef**2
            
            ! jw
            Cdb(j) = rtemp
            ! jw
            
            k = bottom_layer_at_face(j)
				
				! jw
           	Gamma_B(j) = Cdb(j)*dsqrt(un_face(k,j)**2+vn_face(k,j)**2) ! jw
         end if
		end do
		!$omp end parallel do
	else if(bf_flag == 2) then  ! jw
		!$omp parallel do private(j,k,rtemp)
		do j=1,maxface
			if(h_face(j) > dry_depth) then
				k = bottom_layer_at_face(j)
				! jw
				! jw
				rtemp = von_Karman**2/(log(max(0.5,dzhalf_face(k-1,j))/bf_height(j)))**2 ! jw
				
				! jw
				Cdb(j) = rtemp
				! jw
				Gamma_B(j) = Cdb(j)*dsqrt(un_face(k,j)**2+vn_face(k,j)**2) ! jw
			end if
		end do
 		!$omp end parallel do
	end if
	
	! jw
	if(dia_bottom_friction == 1) then
		write(pw_dia_bottom_friction,*) 'it = ', it, ', elapsed_time = ', elapsed_time
		write(pw_dia_bottom_friction,*) 'Cdb(j), Gamma_B(j)'
		do j=1,maxface
			write(pw_dia_bottom_friction,'(A3, I5, 2F10.4)') 'j=', j, Cdb(j), Gamma_B(j)
		end do
	end if
end subroutine calculate_bottom_friction
