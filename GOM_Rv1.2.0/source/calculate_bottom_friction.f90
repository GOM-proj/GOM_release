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
! 	use omp_lib
   implicit none
   
   integer :: j, k
   real(dp):: chezy_coef, hydraulic_radius, expo
   real(dp):: rtemp
   ! end of local variables ==================================================!

	! re-initialize following variables:	
	Cdb 		= 0.0_dp	! bottom drag coefficient
	Gamma_B 	= 0.0_dp ! gamma in equation
	
	if(bf_flag == 0) then	! no bottom friction
		! do nothing
	else if(bf_flag == 1) then	! Chezy-Manning equation
		expo = 1.0_dp/6.0_dp ! exponent part for Chezy equation
		!$omp parallel do private(j,k,hydraulic_radius,chezy_coef,rtemp)
		do j=1,maxface
         if(h_face(j) > dry_depth) then
            ! hydraulic_radius = h_face(j)	! In shallow estuaries, the hydraulic radius can be approximated by the total depth, H
            hydraulic_radius = max(1.0,h_face(j)) ! to prevent excessive bottom friction
            chezy_coef = (hydraulic_radius**expo) / manning(j) ! Cz = (1/n)*R^(1/6)
            rtemp = gravity / chezy_coef**2
            
            ! select either one here...
            Cdb(j) = rtemp
            ! Cdb(j) = min(rtemp,0.0025)
            
            k = bottom_layer_at_face(j)
				
				! bottom friction
           	Gamma_B(j) = Cdb(j)*dsqrt(un_face(k,j)**2+vn_face(k,j)**2) ! you can also use true u & v velocities, but both methods should be identical            
         end if
		end do
		!$omp end parallel do
	else if(bf_flag == 2) then  ! Quadratic fricition law (log law)
		!$omp parallel do private(j,k,rtemp)
		do j=1,maxface
			if(h_face(j) > dry_depth) then
				k = bottom_layer_at_face(j)
				! rtemp = von_Karman**2/(log(dzhalf_face(k-1,j)/bf_height(j)))**2 ! dzhalf_face locates at vertical level, and we have to use the height of the lower level at the bottom layer
				! this is to prevent excessive bottom friction
				rtemp = von_Karman**2/(log(max(0.5,dzhalf_face(k-1,j))/bf_height(j)))**2 ! dzhalf_face locates at vertical level, and we have to use the height of the lower level at the bottom layer
				
				! select either one here...
				Cdb(j) = rtemp
				! Cdb(j) = min(rtemp,0.0025)
				Gamma_B(j) = Cdb(j)*dsqrt(un_face(k,j)**2+vn_face(k,j)**2) ! you can also use true u & v velocities, but both methods should be identical
			end if
		end do
 		!$omp end parallel do
	end if
	
	! write diagnostic file ===================================================!
	if(dia_bottom_friction == 1) then
		write(pw_dia_bottom_friction,*) 'it = ', it, ', elapsed_time = ', elapsed_time
		write(pw_dia_bottom_friction,*) 'Cdb(j), Gamma_B(j)'
		do j=1,maxface
			write(pw_dia_bottom_friction,'(A3, I5, 2F10.4)') 'j=', j, Cdb(j), Gamma_B(j)
		end do
	end if
end subroutine calculate_bottom_friction
