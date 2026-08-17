!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! This subroutine includes turbulence closure models (k-epsilon)
!! 
subroutine turbulence_closure_ke
	use mod_global_variables
	use mod_file_definition
	
	implicit none
	integer :: i, j, k, kk, l
	real(dp):: rtemp01, rtemp02
	integer :: bm, tM, num_vertical_layers, num_vertical_levels
	real(dp):: rtemp1, rtemp2
	real(dp):: tau_bx, tau_by, tau_b, utau_b, tau_sx, tau_sy, tau_s, utau_s
   real(dp),dimension(maxlayer+1)  	:: a_lower_mat1, b_diagonal_mat1, c_upper_mat1, gam1 ! jw
   real(dp),dimension(maxlayer+1,1)	:: solution1, rhs1
   real(dp),dimension(maxlayer+1)  	:: a_lower_mat2, b_diagonal_mat2, c_upper_mat2, gam2 ! jw
   real(dp),dimension(maxlayer+1,1)	:: solution2, rhs2
   real(dp),dimension(maxlayer+1)	:: Ps, Pb, e_sink
   real(dp):: sum1

	! jw
! jw
! jw
! jw
	
	! jw
	rtemp01 = 1.0_dp/sigma_k
	rtemp02 = 1.0_dp/sigma_e

	!$omp parallel do private(j,k,kk,bm,tM,num_vertical_layers,num_vertical_levels) &
	!$omp& private(tau_bx,tau_by,tau_b,utau_b,tau_sx,tau_sy,tau_s,utau_s) &
	!$omp& private(rhs1,a_lower_mat1,b_diagonal_mat1,c_upper_mat1,gam1,solution1) &
	!$omp& private(rhs2,a_lower_mat2,b_diagonal_mat2,c_upper_mat2,gam2,solution2,e_sink) &
	!$omp& private(Ps,Pb,rtemp1,rtemp2)
	do j=1,maxface
		! jw
		a_lower_mat1 		= 0.0_dp
		b_diagonal_mat1	= 0.0_dp
		c_upper_mat1 		= 0.0_dp
		gam1 					= 0.0_dp
		solution1			= 0.0_dp
		rhs1					= 0.0_dp
		a_lower_mat2 		= 0.0_dp
		b_diagonal_mat2	= 0.0_dp
		c_upper_mat2 		= 0.0_dp
		gam2 					= 0.0_dp
		solution2			= 0.0_dp
		rhs2					= 0.0_dp

		Ps						= 0.0_dp ! jw
		Pb						= 0.0_dp ! jw
		e_sink				= 0.0_dp ! jw
		! jw

		bm = bottom_layer_at_face(j)
		tM = top_layer_at_face(j)
				
		num_vertical_layers = tM - bm + 1		
		num_vertical_levels = num_vertical_layers + 1 ! jw
		
		! jw
      if(top_layer_at_face(j) == 0) then
      	cycle ! jw
      end if
      
      
      if(num_vertical_layers == 1) then
      	! jw
      	cycle ! jw
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
		tau_bx = Gamma_B(j)*rho_o*un_face(bm,j)
		tau_by = Gamma_B(j)*rho_o*vn_face(bm,j)
		tau_b = sqrt(tau_bx**2 + tau_by**2)
		utau_b = sqrt(tau_b/rho_o) ! jw
		TKE_k_face(bm-1,j) = MAX(utau_b**2/sqrt(C_mu), k_min) ! jw
		
		! jw
		TKE_e_face(bm-1,j) = MAX(utau_b**3/(von_Karman*(0.5*dzhalf_face(bm-1,j)+bf_height(j))), eps_min) ! jw
				
		! jw
		! jw
		kk = num_vertical_levels

		rhs1(kk,1) = TKE_k_face(bm-1,j)
		a_lower_mat1(kk) = 0.0_dp
		b_diagonal_mat1(kk) = 1.0_dp
		c_upper_mat1(kk) = 0.0_dp

		rhs2(kk,1) = TKE_e_face(bm-1,j)
		a_lower_mat2(kk) = 0.0_dp
		b_diagonal_mat2(kk) = 1.0_dp
		c_upper_mat2(kk) = 0.0_dp
		
		
		! jw
		! jw
		tau_sx = wind_stress_normal(j) ! jw
		tau_sy = wind_stress_tangnt(j) ! jw
		tau_s = sqrt(tau_sx**2 + tau_sy**2)
		utau_s = sqrt(tau_s/rho_o)
		TKE_k_face(tM,j) = MAX(tau_s/(rho_o*sqrt(C_mu)), k_min)
		
		! jw
		! jw
		TKE_e_face(tM,j) = MAX(utau_s**3/(von_Karman*(0.5*dzhalf_face(tM,j)+0.0002_dp)), eps_min) ! jw
		
		! jw
		! jw
		rhs1(1,1) = TKE_k_face(tM,j) ! jw
		a_lower_mat1(1) = 0.0_dp
		b_diagonal_mat1(1) = 1.0_dp
		c_upper_mat1(1) = 0.0_dp

		rhs2(1,1) = TKE_e_face(tM,j) ! jw
		a_lower_mat2(1) = 0.0_dp
		b_diagonal_mat2(1) = 1.0_dp
		c_upper_mat2(1) = 0.0_dp
		
		! jw
		! jw
		! jw
		! jw
		! jw

      ! jw
      do k = bottom_layer_at_face(j), top_layer_at_face(j)-1	! jw
      	! jw
      	! jw
      	! jw
      	! jw
      	! jw
      	! jw
         kk = top_layer_at_face(j) - k + 1 	! jw
         
			! jw
         Ps(kk) = dt*Av_face(k,j)*(((un_face(k+1,j)-un_face(k,j))**2 + (vn_face(k+1,j)-vn_face(k,j))**2)/(dzhalf_face(k,j))**2)
         
         rhs1(kk,1) = rhs1(kk,1) + Ps(kk)
         rhs2(kk,1) = rhs2(kk,1) + (TKE_e_face(k,j)/max(TKE_k_face(k,j),k_min))*c1_eps*Ps(kk) ! jw
			
			! jw
			Pb(kk) = dt*Kv_face(k,j)*(gravity/rho_o)*((rho_face(k+1,j)-rho_face(k,j))/dzhalf_face(k,j))
			if(Pb(kk) >= 0.0_dp) then
				! jw
				rhs1(kk,1) = rhs1(kk,1) + Pb(kk)
				rhs2(kk,1) = rhs2(kk,1) + (TKE_e_face(k,j)/max(TKE_k_face(k,j),k_min))*c3_eps*Pb(kk)
			else ! jw
				! jw
				b_diagonal_mat1(kk) = b_diagonal_mat1(kk) - &
				&								Pb(kk)/max(TKE_k_face(k,j),k_min)
				b_diagonal_mat2(kk) = b_diagonal_mat2(kk) - &
 				&								(TKE_e_face(k,j)/max(TKE_k_face(k,j),k_min))*c3_eps*(Pb(kk)/max(TKE_k_face(k,j),k_min))
			end if
			
			! jw
			! jw
 			e_sink(kk) = dt*TKE_e_face(k,j)/max(TKE_k_face(k,j),k_min) ! jw
			b_diagonal_mat1(kk) = b_diagonal_mat1(kk) + e_sink(kk)
			b_diagonal_mat2(kk) = b_diagonal_mat2(kk) + c2_eps * e_sink(kk)
			
			! jw
			rtemp1 = 1.0_dp/dzhalf_face(k,j)
			! jw


			a_lower_mat1(kk)		= 	a_lower_mat1(kk) - dt*rtemp1*rtemp01* &
			&							  	(Av_face(k+1,j)+Av_face(k,j))/(2*dz_face(k+1,j))
			b_diagonal_mat1(kk)	= 	b_diagonal_mat1(kk) + &
			&							  	1.0_dp + dt*rtemp1*rtemp01* &
			&								((Av_face(k+1,j)+Av_face(k,j))/(2*dz_face(k+1,j)) + (Av_face(k,j)+Av_face(k-1,j))/(2*dz_face(k,j)))
			c_upper_mat1(kk)		= 	c_upper_mat1(kk) - dt*rtemp1*rtemp01* &
			&								(Av_face(k,j)+Av_face(k-1,j))/(2*dz_face(k,j))
			rhs1(kk,1) 				= 	rhs1(kk,1) + TKE_k_face(k,j)

			a_lower_mat2(kk)		= 	a_lower_mat2(kk) - dt*rtemp1*rtemp02* &
			&								(Av_face(k+1,j)+Av_face(k,j))/(2*dz_face(k+1,j))
			b_diagonal_mat2(kk)	= 	b_diagonal_mat2(kk) + &
			&							  	1.0_dp + dt*rtemp1*rtemp02* &
			&								((Av_face(k+1,j)+Av_face(k,j))/(2*dz_face(k+1,j)) + (Av_face(k,j)+Av_face(k-1,j))/(2*dz_face(k,j)))
			c_upper_mat2(kk)		= 	c_upper_mat2(kk) - dt*rtemp1*rtemp02* &
			&								(Av_face(k,j)+Av_face(k-1,j))/(2*dz_face(k,j))
			rhs2(kk,1) 				= 	rhs2(kk,1) + TKE_e_face(k,j)			
      end do
      
		! jw
		call tridiagonal_solver(maxlayer+1,num_vertical_levels,1, a_lower_mat1, b_diagonal_mat1, c_upper_mat1, rhs1, solution1, gam1)
      call tridiagonal_solver(maxlayer+1,num_vertical_levels,1, a_lower_mat2, b_diagonal_mat2, c_upper_mat2, rhs2, solution2, gam2)
			
! jw
! jw
! jw
! jw
! jw
		
		! jw
		do k=tM,bm-1,-1
	
			kk = top_layer_at_face(j) - k + 1
			
			TKE_k_face(k,j) = MAX(solution1(kk,1), k_min)
			TKE_e_face(k,j) = MAX(solution2(kk,1), eps_min)
			
			! jw
			Av_face(k,j) = C_mu*TKE_k_face(k,j)**2/TKE_e_face(k,j)
			Kv_face(k,j) = Av_face(k,j)/sigma_t
			
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
	end do ! jw
	!omp end parallel do
		
	! jw
	! jw
	!$omp parallel do private(i,j,k,l,sum1)
	do i=1,maxele
		do k=1,maxlayer
			sum1 = 0.0_dp
			do l = 1,tri_or_quad(i)
				j = facenum_at_cell(l,i)
				sum1 = sum1 + Kv_face(k,j)
			end do
			Kv_cell(k,i) = sum1/tri_or_quad(i)
		end do
	end do
	!omp end parallel do
		
! jw
! jw
! jw
end subroutine turbulence_closure_ke