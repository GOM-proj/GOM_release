!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! solve transport equation with TVD
!! 
subroutine solve_transport_equation_v19
	use mod_global_variables
	! jw
	implicit none

	integer :: i, j, k, l, ii, n1, n2
	integer :: i2, j2, ii2, i3, ie
	integer :: kk, num_vertical_layer
	integer :: t_layer, b_layer
	integer :: count1

	real(dp),dimension(0:maxlayer+1,maxele,maxtran2) :: con_cell ! jw
	real(dp),dimension(maxlayer,maxele,maxtran2) :: con_cell_new
	real(dp),dimension(maxlayer,maxnod,maxtran2) :: con_node_new
	real(dp),dimension(maxlayer,maxface,maxtran2) :: con_face_new

	! jw
	integer :: ith_ob, ith_Qb, ith_WR, ith_SS ! jw
	real(dp),dimension(num_ob_cell,maxtran2) :: con_at_ob
	real(dp),dimension(num_Qb_cell,maxtran2) :: con_at_Qb
	real(dp),dimension(num_WR_cell,maxtran2) :: con_at_WR
	real(dp),dimension(num_SS_cell,maxtran2) :: con_at_SS
	real(dp),dimension(maxlayer,num_ob_cell,maxtran2) :: con_at_obk ! jw
	real(dp),dimension(maxlayer,num_Qb_cell,maxtran2) :: con_at_Qbk ! jw
	real(dp),dimension(maxlayer,num_WR_cell,maxtran2) :: con_at_WRk
	real(dp),dimension(maxlayer,num_SS_cell,maxtran2) :: con_at_SSk
	
	real(dp),dimension(maxlayer,maxtran2):: con_cell_iik_oQb ! jw
	real(dp),dimension(maxtran2):: con_cell_iik ! jw

	! jw
	real(dp),dimension(maxlayer,maxface):: Q_jk_theta, D_jk_old, d_jk_theta ! jw
	real(dp),dimension(0:maxlayer,maxele):: Q_ik_theta, D_ik_old, d_ik_theta
	real(dp) :: Q_jk_theta2(4) ! jw

	! jw
	real(dp),dimension(maxlayer,maxface,maxtran2) :: r_jk
	real(dp),dimension(0:maxlayer,maxele,maxtran2):: r_ik
	real(dp),dimension(maxtran2):: sum1, sum1_all, sum2, sum2_all
	
	! jw
	real(dp),dimension(maxlayer,maxface) :: small_phi_jk
	real(dp),dimension(0:maxlayer,maxele):: small_phi_ik
	real(dp),dimension(maxlayer,maxface,maxtran2):: large_PHI_jk 
	real(dp),dimension(0:maxlayer,maxele,maxtran2):: large_PHI_ik
	real(dp),dimension(maxlayer,maxface,maxtran2):: PSI_jk
	real(dp),dimension(0:maxlayer,maxele,maxtran2):: PSI_ik
	
	! jw
	real(dp),dimension(maxtran2):: h_flux_term, v_flux_term
	real(dp),dimension(maxtran2):: h_adv, v_adv, h_diff
	real(dp),dimension(maxtran2):: sum_QC_i, sum_QC_ii, sum_hdiff, sum_hlimiter

	! jw
	real(dp),dimension(maxlayer,maxtran2):: rhs, solution
	real(dp),dimension(maxlayer):: a_lower_mat, b_diagonal_mat, c_upper_mat, gam

	! jw
	integer :: t2
	real(dp):: u1,u2,u3,v1,v2,v3
	
	! jw
	integer :: subcycle, Nt
	real(dp):: dt2
	
	! jw
	real(dp):: cp = 4186.0 ! jw
	real(dp):: dz_old, dz_new ! jw
	real(dp):: rtemp, rtemp1, rtemp2
	real(dp):: sum11, sum22
	! jw
		
	! jw
	con_cell_new = 0.0_dp
	con_node_new = 0.0_dp
	con_face_new = 0.0_dp
	
	! jw
	! jw
	! jw
	! jw
	do i3=1,maxtran2
		if(tran_id(i3) == 1) then
			do i=1,maxele
				do k=1,maxlayer
					con_cell(k,i,i3) = salt_cell(k,i)
				end do
	 			! jw
				con_cell(0,i,i3) = 0.0_dp
				con_cell(maxlayer+1,i,i3) = 0.0_dp
			end do
		end if
		if(tran_id(i3) == 2) then
			if(heat_option == 1) then
				! jw
				do i=1,maxele
					do k=1,maxlayer
						con_cell(k,i,i3) = temp_cell(k,i)
					end do
		 			! jw
					con_cell(0,i,i3) = 0.0_dp
					con_cell(maxlayer+1,i,i3) = 0.0_dp
					
					! jw
					t_layer = top_layer_at_element(i)
					b_layer = bottom_layer_at_element(i)
					
					
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
					
					
					do k=b_layer,t_layer
						con_cell(k,i,i3) = con_cell(k,i,i3) + dt*phi_sz(k,i)/(rho_o*cp*dz_cell(k,i))
					end do
					
					
					! jw
					! jw
					! jw
				end do
			else if(heat_option == 2) then ! jw
				! jw
				do i=1,maxele
					do k=1,maxlayer
						con_cell(k,i,i3) = temp_cell(k,i)
					end do
		 			! jw
					con_cell(0,i,i3) = 0.0_dp
					con_cell(maxlayer+1,i,i3) = 0.0_dp
					
					! jw
					t_layer = top_layer_at_element(i)
					b_layer = bottom_layer_at_element(i)
					
					! jw
					con_cell(t_layer,i,i3) = con_cell(t_layer,i,i3) ! jw

					! jw
					do k=1,maxlayer
						con_cell(k,i,i3) = con_cell(k,i,i3) + dt*phi_sw(i)/(rho_o*cp*dz_cell(t_layer,i))
					end do
	
					! jw
					! jw
					! jw
				end do				
			end if
			
			! jw
			! jw
			! jw
			! jw
			! jw
		end if
	end do
	

	Nt = trans_sub_iter ! jw
	
	! jw
	! jw
	! jw
	if(Nt == 0) then
		Nt = 1
	end if

	! jw
	
	
	dt2 = dt/Nt ! jw
	do subcycle=0,(Nt-1) ! jw
		! jw
		! jw
		! jw
		! jw
		! jw
		
		
		! jw
		! jw
		con_at_ob = 0.0_dp
		con_at_Qb = 0.0_dp
		con_at_WR = 0.0_dp
		con_at_SS = 0.0_dp
		con_at_obk = 0.0_dp
		con_at_Qbk = 0.0_dp
		con_at_WRk = 0.0_dp
		con_at_SSk = 0.0_dp

		! jw
  		r_ik = 1.0_dp
 		r_jk = 1.0_dp

		! jw
		! jw
	   a_lower_mat = 0.0_dp
	   b_diagonal_mat = 0.0_dp
	   c_upper_mat = 0.0_dp
	   gam = 0.0_dp   
		rhs = 0.0_dp
		solution = 0.0_dp

	   ! jw
	   ! jw
	   ! jw
	   ! jw
	   ! jw
		! jw
		! jw
		! jw
		u2 = (julian_day*86400.0 - dt) + (subcycle+1)*dt2 ! jw
		
		! jw
		do i=1,num_ob_cell
			! jw
			! jw
			! jw
			! jw
			! jw
			
			! jw
			! jw
			do i3=1,maxtran2
				! jw
				if(tran_id(i3) == 1 .and. salt_ser_id(i) > 0) then
					ii = salt_ser_id(i)
					do t2=2,salt_ser_data_num(ii)
						! jw
						! jw
						! jw
						u1 = salt_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! jw
						u3 = salt_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! jw
						
						if(u2 >= u1 .and. u2 <= u3) then
							v1 = salt_ser_salt(t2-1,ii)		! jw
							v3 = salt_ser_salt(t2  ,ii)		! jw
							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
							exit
						end if			
					end do
					! jw
					con_at_ob(i,i3) = spinup_function_baroclinic * v2
				end if
			
				! jw
				if(tran_id(i3) == 2 .and. temp_ser_id(i) > 0) then
					ii = temp_ser_id(i)
					do t2=2,temp_ser_data_num(ii)
						! jw
						! jw
						! jw
						u1 = temp_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! jw
						u3 = temp_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! jw
						
						if(u2 >= u1 .and. u2 <= u3) then
							v1 = temp_ser_temp(t2-1,ii)		! jw
							v3 = temp_ser_temp(t2  ,ii)		! jw
							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
							exit
						end if			
					end do
					! jw
					con_at_ob(i,i3) = spinup_function_baroclinic * v2 ! jw
				end if
			end do ! jw

			! jw
			! jw
			! jw
			! jw
			! jw
			! jw
			i2 = ob_cell_id(i) ! jw
			! jw
			t_layer = top_layer_at_element(i2)
			b_layer = bottom_layer_at_element(i2)

			do l=1,tri_or_quad(i2)
				j2 = facenum_at_cell(l,i2)
				if(boundary_type_of_face(j2) > 0) then ! jw
					do k=b_layer,t_layer
						! jw
						! jw
						! jw
						! jw
						dz_old = (dz_face(k,j2) + (dz_face_new(k,j2) - dz_face(k,j2))*(subcycle)/Nt)
						Q_jk_theta2(l) = face_length(j2)* dz_old &
						&				*((1.0-theta)*un_face(k,j2)*sign_in_outflow(l,i2) + theta*un_face_new(k,j2)*sign_in_outflow(l,i2))
						
						if(Q_jk_theta2(l) <= 0.0_dp) then ! jw
							do i3=1,maxtran2
								con_at_obk(k,i,i3) = con_at_ob(i,i3) ! jw
							end do
						else ! jw
							do i3=1,maxtran2
								con_at_obk(k,i,i3) = con_cell(k,i2,i3) ! jw
								! jw
							end do
						end if
					end do
				end if
			end do ! jw
		end do ! jw

		! jw
		! jw
		! jw
		! jw
		do i=1,num_Qb_cell
			! jw
			! jw
			! jw
			! jw
			
			! jw
			do i3=1,maxtran2
				! jw
				if(tran_id(i3) == 1 .and. Q_salt_ser_id(i) > 0) then
					ii = Q_salt_ser_id(i)
					do t2=2,salt_ser_data_num(ii)
						! jw
						! jw
						! jw
						u1 = salt_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! jw
						u3 = salt_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! jw
						
						if(u2 >= u1 .and. u2 <= u3) then
							v1 = salt_ser_salt(t2-1,ii) 		! jw
							v3 = salt_ser_salt(t2  ,ii) 		! jw
							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
							exit
						end if
					end do
					! jw
					con_at_Qb(i,i3) = spinup_function_baroclinic * v2
					
					! jw
				end if
				
				! jw
				if(tran_id(i3) == 2 .and. Q_temp_ser_id(i) > 0) then
					ii = Q_temp_ser_id(i)
					do t2=2,temp_ser_data_num(ii)
						! jw
						! jw
						! jw
						u1 = temp_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! jw
						u3 = temp_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! jw
						
						if(u2 >= u1 .and. u2 <= u3) then
							v1 = temp_ser_temp(t2-1,ii) 		! jw
							v3 = temp_ser_temp(t2  ,ii) 		! jw
							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
							exit
						end if
					end do
					! jw
					con_at_Qb(i,i3) = spinup_function_baroclinic * v2 ! jw
					! jw
					! jw
				end if
			end do ! jw
			
			! jw
			! jw
			! jw
			i2 = Q_boundary(i,1) ! jw
			
			t_layer = top_layer_at_element(i2)
			b_layer = bottom_layer_at_element(i2)
			
			do k=b_layer,t_layer
				do i3=1,maxtran2
					con_at_Qbk(k,i,i3) = con_at_Qb(i,i3)
					
					! jw
				end do
			end do
		end do ! jw

		! jw
		! jw
		! jw
		do i=1,num_WR_cell
			! jw
			! jw
			! jw
			! jw
			
			! jw
			do i3=1,maxtran2
				! jw
				if(tran_id(i3) == 1 .and. WR_salt_ser_id(i) > 0) then
					ii = WR_salt_ser_id(i)
					do t2=2,salt_ser_data_num(ii)
						! jw
						! jw
						! jw
						u1 = salt_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! jw
						u3 = salt_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! jw
						
						if(u2 >= u1 .and. u2 <= u3) then
							v1 = salt_ser_salt(t2-1,ii) 		! jw
							v3 = salt_ser_salt(t2  ,ii) 		! jw
							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
							exit
						end if
					end do
					! jw
					con_at_WR(i,i3) = spinup_function_baroclinic * v2
					! jw
				end if
				
				! jw
				if(tran_id(i3) == 2 .and. WR_temp_ser_id(i) > 0) then
					ii = WR_temp_ser_id(i)
					do t2=2,temp_ser_data_num(ii)
						! jw
						! jw
						! jw
						u1 = temp_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! jw
						u3 = temp_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! jw
						
						if(u2 >= u1 .and. u2 <= u3) then
							v1 = temp_ser_temp(t2-1,ii) 		! jw
							v3 = temp_ser_temp(t2  ,ii) 		! jw
							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
							exit
						end if
					end do
					! jw
					con_at_WR(i,i3) = spinup_function_baroclinic * v2 ! jw
					! jw
					! jw
				end if
			end do ! jw

			! jw
			! jw
			! jw
			i2 = WR_boundary(i,1)! jw
			j = WR_boundary(i,2) ! jw

			if(WR_layer(i) == 999) then
				! jw
				k = top_layer_at_face(j)
			else if(WR_layer(i) == 0) then
				! jw
				k = bottom_layer_at_face(j)
			else
				k = WR_layer(i)
			end if
			
			if(WRu_boundary(i) >= 0.0) then ! jw
				! jw
				do i3=1,maxtran2
					con_at_WRk(k,i,i3) = con_cell(k,i2,i3)
					! jw
				end do
			else ! jw
				! jw
				do i3=1,maxtran2
					con_at_WRk(k,i,i3) = con_at_WR(i,i3)
					! jw
				end do
			end if
		end do ! jw
		

		! jw
 		do i=1,num_SS_cell			
 			! jw
 			do i3=1,maxtran2
 				! jw
 				if(tran_id(i3) == 1 .and. SS_salt_ser_id(i) > 0) then
 					ii = SS_salt_ser_id(i)
 					do t2=2,salt_ser_data_num(ii)
 						! jw
 						! jw
 						! jw
 						u1 = salt_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! jw
 						u3 = salt_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! jw
 						
 						if(u2 >= u1 .and. u2 <= u3) then
 							v1 = salt_ser_salt(t2-1,ii) 		! jw
 							v3 = salt_ser_salt(t2  ,ii) 		! jw
 							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
 							exit
 						end if
 					end do
 					! jw
 					con_at_SS(i,i3) = spinup_function_baroclinic * v2
 				end if
 				
 				! jw
 				if(tran_id(i3) == 2 .and. SS_temp_ser_id(i) > 0) then
 					ii = SS_temp_ser_id(i)
 					do t2=2,temp_ser_data_num(ii)
 						! jw
 						! jw
 						! jw
 						u1 = temp_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! jw
 						u3 = temp_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! jw
 						
 						if(u2 >= u1 .and. u2 <= u3) then
 							v1 = temp_ser_temp(t2-1,ii) 		! jw
 							v3 = temp_ser_temp(t2  ,ii) 		! jw
 							v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
 							exit
 						end if
 					end do
 					! jw
 					con_at_SS(i,i3) = spinup_function_baroclinic * v2
 				end if
 			end do ! jw
 		end do ! jw
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
 		!$omp parallel
		! jw
		!$omp do private(j, k, dz_old, b_layer, t_layer)
			do j=1,maxface
				b_layer = bottom_layer_at_face(j)
				t_layer = top_layer_at_face(j)
				
				do k=b_layer,t_layer
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
					dz_old = dz_face(k,j) + (dz_face_new(k,j) - dz_face(k,j))*(subcycle)/Nt
					Q_jk_theta(k,j) = face_length(j)*dz_old &
					&						*((1.0-theta)*un_face(k,j) + theta*un_face_new(k,j))
					D_jk_old(k,j) = face_length(j) * dz_old * Kh(k,j) &
					&					/delta_j(j)
					d_jk_theta(k,j) = max(0.0, D_jk_old(k,j) - 0.5*abs(Q_jk_theta(k,j)))
					
					! jw
					if(Q_jk_theta(k,j) == 0.0_dp) then
						small_phi_jk(k,j) = 1.0_dp
					else
						small_phi_jk(k,j) = min(1.0, 2.0*D_jk_old(k,j)/abs(Q_jk_theta(k,j)))
					end if
				end do ! jw
			end do ! jw
		!$omp end do
		
		! jw
		!$omp do private(i, k, dz_old, b_layer, t_layer)
			do i=1,maxele
				b_layer = bottom_layer_at_element(i)
				t_layer = top_layer_at_element(i)
				! jw
				! jw
				do k=b_layer,t_layer-1 ! jw
					! jw
					! jw
					! jw
					Q_ik_theta(k,i) = area(i)*((1.0-theta)*wn_cell(  k,i) + theta*wn_cell_new(  k,i))
					
					! jw
					dz_old = dzhalf_cell(k,i) + (dzhalf_cell_new(k,i) - dzhalf_cell(k,i))*(subcycle)/Nt
					D_ik_old(k,i) = area(i)*Kv(k,i)/dz_old
					d_ik_theta(k,i) = max(0.0, D_ik_old(k,i) - 0.5*abs(Q_ik_theta(k,i)))
						
					! jw
					if(Q_ik_theta(k,i) == 0.0_dp) then
						small_phi_ik(k,i) = 1.0_dp
					else
						small_phi_ik(k,i) = min(1.0, 2.0*D_ik_old(k,i)/abs(Q_ik_theta(k,i)))
					end if
				end do ! jw
			end do ! jw
		!$omp end do
		!$omp end parallel
		
		! jw
		! jw
		if(h_flux_limiter > 0 .or. v_flux_limiter > 0) then
			!$omp parallel
	 		!$omp do private(i, j, k, l, ii, i3,											&
	 		!$omp & 			  b_layer, t_layer, num_vertical_layer,					&
	 		!$omp &			  ith_ob, ith_Qb,	ith_WR,										&
	 		!$omp &			  con_cell_iik, con_cell_iik_oQb,							&
	 		!$omp & 			  sum1, sum1_all, sum2, sum2_all,							&
	 		!$omp &			  Q_jk_theta2,														&
	 		!$omp &			  rtemp, rtemp1, rtemp2)
			do i=1,maxele
				t_layer = top_layer_at_element(i)
				b_layer = bottom_layer_at_element(i)
				
				if(t_layer == 0) then
					! jw
					cycle
				end if
				num_vertical_layer = t_layer - b_layer + 1 ! jw
				
				! jw
				! jw
				! jw
				! jw
				ith_ob = ob_element_flag(i) ! jw
				ith_Qb = Qb_element_flag(i) ! jw
				ith_WR = WR_element_flag(i) ! jw
				
				if(ith_ob > 0) then ! jw
		 			do k=b_layer,t_layer
						do i3=1,maxtran2
							con_cell_iik_oQb(k,i3) = con_at_obk(k,ith_ob,i3) ! jw
						end do
					end do
				end if
				if(ith_Qb > 0) then ! jw
					do k=b_layer,t_layer
						do i3=1,maxtran2
							con_cell_iik_oQb(k,i3) = con_at_Qbk(k,ith_Qb,i3) ! jw
							! jw
						end do
					end do
				end if
				if(ith_WR > 0) then ! jw
					! jw
					j = WR_boundary(i,2) ! jw
					if(WR_layer(ith_WR) == 999) then
						! jw
						k = top_layer_at_face(j)
					else if(WR_layer(ith_WR) == 0) then
						! jw
						k = bottom_layer_at_face(j)
					else
						! jw
						k = WR_layer(ith_WR)
					end if
					
					do i3=1,maxtran2
						con_cell_iik_oQb(k,i3) = con_at_WRk(k,ith_WR,i3)
						! jw
					end do
				end if

				! jw
				do k=b_layer,t_layer
					! jw
					! jw
					sum1 		= 0.0_dp ! jw
					sum2 		= 0.0_dp ! jw
					sum1_all = 0.0_dp ! jw
					sum2_all = 0.0_dp ! jw
					
					! jw
					! jw
					! jw
					if(num_vertical_layer == 1) then
						! jw
					else
						! jw
						if(k == b_layer) then
							do i3=1,maxtran2
								rtemp2 = abs(0.5*(Q_ik_theta(k  ,i) + abs(Q_ik_theta(k  ,i))))
								sum1_all(i3) = rtemp2 * (con_cell(k,i,i3) - con_cell(k+1,i,i3))
								sum2_all(i3) = rtemp2
							end do
						else if(k == t_layer) then
							do i3=1,maxtran2
								rtemp1 = abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))
								sum1_all(i3) = rtemp1 * (con_cell(k,i,i3) - con_cell(k-1,i,i3))
								sum2_all(i3) = rtemp1
							end do						
						else
							do i3=1,maxtran2
								rtemp1 = abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))
								rtemp2 = abs(0.5*(Q_ik_theta(k  ,i) + abs(Q_ik_theta(k  ,i))))
								sum1_all(i3) = rtemp1 * (con_cell(k,i,i3) - con_cell(k-1,i,i3)) &
								&				 + rtemp2 * (con_cell(k,i,i3) - con_cell(k+1,i,i3))
								sum2_all(i3) = rtemp1 + rtemp2
							end do
						end if
					end if
					
					! jw
					! jw
					do l=1,tri_or_quad(i)
						j = facenum_at_cell(l,i)
						ii = adj_cellnum_at_cell(l,i)
						if(boundary_type_of_face(j) /= -1 .or. isflowside3(j) > 0 .or. isflowside4(j) > 0) then ! jw
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
							Q_jk_theta2(l) = Q_jk_theta(k,j)*sign_in_outflow(l,i)
								
							if(Q_jk_theta2(l) < 0.0_dp) then ! jw
								! jw
								! jw
								! jw
								if(ii > 0) then	! jw
									do i3=1,maxtran2
										con_cell_iik(i3) = con_cell(k,ii,i3)
									end do
								else ! jw
									do i3=1,maxtran2
										con_cell_iik(i3) = con_cell_iik_oQb(k,i3) ! jw
									end do
								end if
								
								do i3=1,maxtran2 ! jw
									sum1(i3) = sum1(i3) + abs(Q_jk_theta2(l))*(con_cell(k,i,i3) - con_cell_iik(i3))
									sum2(i3) = sum2(i3) + abs(Q_jk_theta2(l))
								end do
							end if ! jw
						end if ! jw
					end do ! jw
					
					! jw
					do i3=1,maxtran2
						sum1_all(i3) = sum1_all(i3) + sum1(i3)
						sum2_all(i3) = sum2_all(i3) + sum2(i3)
					end do
					
					! jw
					do l=1,tri_or_quad(i)
						j = facenum_at_cell(l,i)
						ii = adj_cellnum_at_cell(l,i)
						
						if(boundary_type_of_face(j) /= -1 .or. isflowside3(j) > 0 .or. isflowside4(j) > 0) then ! jw
							! jw
							! jw
							! jw
							if(Q_jk_theta2(l) >= 0.0_dp) then ! jw
								if(ii > 0) then	! jw
									do i3=1,maxtran2
										con_cell_iik(i3) = con_cell(k,ii,i3)
									end do
								else ! jw
									do i3=1,maxtran2
										con_cell_iik(i3) = con_cell_iik_oQb(k,i3) ! jw
									end do
								end if
								
								! jw
								do i3=1,maxtran2
									rtemp = sum2_all(i3)*(con_cell_iik(i3) - con_cell(k,i,i3))
									if(rtemp == 0.0_dp) then ! jw
										r_jk(k,j,i3) = 1.0 ! jw
									else
										r_jk(k,j,i3) = sum1_all(i3)/rtemp
									end if
								end do
							end if
						end if
					end do ! jw
				end do ! jw
				! jw
				
				! jw
				! jw
				! jw
				if(Q_ik_theta(t_layer,i) >= 0.0_dp) then ! jw
					do k=b_layer,t_layer-1
						do i3=1,maxtran2
							rtemp1 = sum1(i3) + abs(Q_ik_theta(k-1,i)) * (con_cell(k,i,i3)-con_cell(k-1,i,i3))
							rtemp2 = sum2(i3) + abs(Q_ik_theta(k-1,i))
							rtemp  = (con_cell(k+1,i,i3)-con_cell(k,i,i3))*rtemp2
							if(rtemp == 0.0) then ! jw
								r_ik(k,i,i3) = 1.0_dp
							else
								r_ik(k,i,i3) = rtemp1/rtemp
							end if						
						end do
					end do
					do i3=1,maxtran2
						r_ik(t_layer  ,i,i3) = 1.0_dp
						r_ik(b_layer-1,i,i3) = 1.0_dp
					end do
				else ! jw
					do k=b_layer+1,t_layer
						do i3=1,maxtran2
							rtemp1 = sum1(i3) + abs(Q_ik_theta(k,i)) * (con_cell(k,i,i3)-con_cell(k+1,i,i3))
							rtemp2 = sum2(i3) + abs(Q_ik_theta(k,i))
							rtemp  = (con_cell(k-1,i,i3)-con_cell(k,i,i3))*rtemp2
							if(rtemp == 0.0_dp) then ! jw
								r_ik(k-1,i,i3) = 1.0_dp
							else
								r_ik(k-1,i,i3) = rtemp1/rtemp
							end if
						end do
					end do
					do i3=1,maxtran2
						r_ik(t_layer  ,i,i3) = 1.0_dp
						r_ik(b_layer-1,i,i3) = 1.0_dp
					end do
				end if
				! jw
			end do ! jw
			!$omp end do
			
			! jw
			!$omp do private(i,k,i3,b_layer,t_layer)
				do i=1,maxele
					b_layer = bottom_layer_at_element(i)
					t_layer = top_layer_at_element(i)
					
					do k=b_layer,t_layer-1 ! jw
						do i3=1,maxtran2
							if(h_flux_limiter == 0) then
								PSI_ik(k,i,i3) = 0.0_dp
							else if(h_flux_limiter == 1) then ! jw
								large_PHI_ik(k,i,i3) = max(small_phi_ik(k,i), &
								&									min(1.0,r_ik(k,i,i3)))
								PSI_ik(k,i,i3) = large_PHI_ik(k,i,i3) - small_phi_ik(k,i)
							else if(h_flux_limiter == 2) then ! jw
								large_PHI_ik(k,i,i3) = max(small_phi_ik(k,i), &
								&									(r_ik(k,i,i3) + abs(r_ik(k,i,i3)))/(1.0+abs(r_ik(k,i,i3))))
								PSI_ik(k,i,i3) = large_PHI_ik(k,i,i3) - small_phi_ik(k,i)
							else if(h_flux_limiter == 3) then ! jw
								large_PHI_ik(k,i,i3) = max(small_phi_ik(k,i), &
								&									min(1.0,2.0*r_ik(k,i,i3)), &
								&									min(2.0,r_ik(k,i,i3)))
								PSI_ik(k,i,i3) = large_PHI_ik(k,i,i3) - small_phi_ik(k,i)
							end if
						end do
					end do
				end do
			!$omp end do
			
			! jw
			!$omp do private(j,k,i3,b_layer,t_layer)
				do j=1,maxface
					b_layer = bottom_layer_at_face(j)
					t_layer = top_layer_at_face(j)
					
					do k=b_layer,t_layer
						do i3=1,maxtran2
							if(h_flux_limiter == 0) then
								PSI_jk(k,j,i3) = 0.0_dp
							else if(h_flux_limiter == 1) then ! jw
								large_PHI_jk(k,j,i3) = max(small_phi_jk(k,j), &
								&									min(1.0,r_jk(k,j,i3)))
								PSI_jk(k,j,i3) = large_PHI_jk(k,j,i3) - small_phi_jk(k,j)
							else if(h_flux_limiter == 2) then ! jw
								large_PHI_jk(k,j,i3) = max(small_phi_jk(k,j), &
								&									(r_jk(k,j,i3) + abs(r_jk(k,j,i3)))/(1.0+abs(r_jk(k,j,i3))))
								PSI_jk(k,j,i3) = large_PHI_jk(k,j,i3) - small_phi_jk(k,j)
							else if(h_flux_limiter == 3) then ! jw
								large_PHI_jk(k,j,i3) = max(small_phi_jk(k,j), &
								&									min(1.0,2.0*r_jk(k,j,i3)), &
								&									min(2.0,r_jk(k,j,i3)))
								PSI_jk(k,j,i3) = large_PHI_jk(k,j,i3) - small_phi_jk(k,j)
							end if
						end do
					end do
				end do
			!$omp end do
			!$omp end parallel
			! jw
		end if ! jw


		! jw
		! jw
		!$omp parallel
		!$omp do private(i,j,k,l,ii,kk,i3, 						&
		!$omp & 			b_layer,t_layer,num_vertical_layer, &
		!$omp &			dz_old, dz_new,							&
		!$omp &			Q_jk_theta2, 								&
		!$omp & 			ith_ob, ith_Qb, ith_WR, ith_SS,		&
		!$omp &			con_cell_iik,							 	&
		!$omp &			con_cell_iik_oQb,							&
		!$omp &			sum_QC_i,	sum_QC_ii,					&
		!$omp &			sum_hdiff, sum_hlimiter, 				&
		!$omp &			h_adv, h_diff, h_flux_term, v_adv, v_flux_term, &
		!$omp &			a_lower_mat,b_diagonal_mat,c_upper_mat,rhs,solution,gam)
		do i=1,maxele
			b_layer = bottom_layer_at_element(i)
			t_layer = top_layer_at_element(i)
			
			! jw
			! jw
			if(t_layer == 0) then
				cycle ! jw
			end if
			
			! jw
			! jw
			num_vertical_layer = t_layer - b_layer + 1 ! jw
			
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
			do k = b_layer+1, t_layer-1 ! jw
				! jw
				! jw
				! jw
				kk = t_layer + 1 - k ! jw
				! jw
				! jw
				
				! jw
				! jw
				a_lower_mat(kk) = -dt2*d_ik_theta(k  ,i)
				c_upper_mat(kk) = -dt2*d_ik_theta(k-1,i)
				! jw
				! jw
				dz_new = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle+1)/Nt
				b_diagonal_mat(kk) = - a_lower_mat(kk) &
				&							+ area(i)*dz_new &
				&							- c_upper_mat(kk)
			end do ! jw
			
			
			! jw
			if(t_layer == b_layer) then ! jw
				! jw
				! jw
				! jw
				! jw
				kk = 1 ! jw
				k = t_layer ! jw
				
				! jw
				! jw
				! jw
				! jw
				dz_new = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle+1)/Nt
				b_diagonal_mat(kk) = area(i)*dz_new
				! jw
			else ! jw
				! jw
				! jw
				! jw
				! jw
				kk = 1
				k = t_layer ! jw
				
				! jw
				! jw
				! jw
				! jw
				c_upper_mat(kk) = -dt2*d_ik_theta(k-1,i)
				dz_new = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle+1)/Nt
				b_diagonal_mat(kk) = area(i)*dz_new &
				&						  -c_upper_mat(kk)
				
				! jw
				! jw
				kk = num_vertical_layer ! jw
				k = b_layer
				
				! jw
	        	! jw
	        	! jw
	      	! jw
	         a_lower_mat(kk) = -dt2*d_ik_theta(k,i)
	         dz_new = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle+1)/Nt
	         b_diagonal_mat(kk) = -a_lower_mat(kk) &
	         &						   + area(i)*dz_new
	         ! jw
	      end if ! jw
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
			ith_ob = ob_element_flag(i) ! jw
			ith_Qb = Qb_element_flag(i) ! jw
			ith_WR = WR_element_flag(i) ! jw
			if(ith_ob > 0) then ! jw
	 			do k=b_layer,t_layer
					do i3=1,maxtran2
						con_cell_iik_oQb(k,i3) = con_at_obk(k,ith_ob,i3) ! jw
					end do
				end do
			end if
			if(ith_Qb > 0) then ! jw
				do k=b_layer,t_layer
					do i3=1,maxtran2
						con_cell_iik_oQb(k,i3) = con_at_Qbk(k,ith_Qb,i3) ! jw
						! jw
					end do
				end do
			end if
			if(ith_WR > 0) then ! jw
				! jw
				j = WR_boundary(i,2) ! jw
				if(WR_layer(ith_WR) == 999) then
					! jw
					k = top_layer_at_face(j)
				else if(WR_layer(ith_WR) == 0) then
					! jw
					k = bottom_layer_at_face(j)
				else
					! jw
					k = WR_layer(ith_WR)
				end if
					
				do i3=1,maxtran2
					con_cell_iik_oQb(k,i3) = con_at_WRk(k,ith_WR,i3)
					! jw
				end do				
			end if

			do k=b_layer,t_layer
				kk = t_layer + 1 - k ! jw
				
	         sum_QC_i = 0.0_dp ! jw
	         sum_QC_ii = 0.0_dp ! jw
	         sum_hdiff = 0.0_dp ! jw
	         sum_hlimiter = 0.0_dp ! jw
				do l=1,tri_or_quad(i) ! jw
					j = facenum_at_cell(l,i)
					ii = adj_cellnum_at_cell(l,i) ! jw

	        		if(ii > 0) then	! jw
	        			do i3=1,maxtran2
	        				con_cell_iik(i3) = con_cell(k,ii,i3)
	        			end do
	        		else ! jw
	        			do i3=1,maxtran2
	        				con_cell_iik(i3) = con_cell_iik_oQb(k,i3)
	        			end do
	        		end if

					! jw
					! jw
					! jw
					! jw
					! jw
					if(boundary_type_of_face(j) /= -1  .or. isflowside3(j) > 0 .or. isflowside4(j) > 0) then ! jw
						! jw
			      	! jw
			      	Q_jk_theta2(l) = Q_jk_theta(k,j)*sign_in_outflow(l,i)

	 					! jw
	 					! jw
	 					! jw
	 					! jw
	 					! jw
						if(Q_jk_theta2(l) > 0.0_dp) then
							! jw
							! jw
							do i3=1,maxtran2
								sum_QC_i(i3) = sum_QC_i(i3) +  abs(Q_jk_theta2(l))*con_cell(k,i,i3) ! jw
							end do
						else if(Q_jk_theta2(l) < 0.0_dp) then ! jw
							! jw
							! jw
							do i3=1,maxtran2
								sum_QC_ii(i3) = sum_QC_ii(i3) + abs(Q_jk_theta2(l))*con_cell_iik(i3) ! jw
							end do
						end if ! jw
						
						! jw
	 					do i3=1,maxtran2
	 						! jw
	 						sum_hdiff(i3) = sum_hdiff(i3) + d_jk_theta(k,j)*(con_cell_iik(i3) - con_cell(k,i,i3))
	 					end do

						! jw
						do i3=1,maxtran2
		 					sum_hlimiter(i3) = sum_hlimiter(i3) + &
		 					&						 PSI_jk(k,j,i3)*abs(Q_jk_theta2(l))*(con_cell(k,i,i3) - con_cell_iik(i3))
		 				end do
					end if ! jw
				end do ! jw
	         
	         ! jw
	         do i3=1,maxtran2
		         ! jw
		         ! jw
		         ! jw
		         h_adv(i3) = dt2*(sum_QC_ii(i3) - sum_QC_i(i3)) ! jw
		         
		         ! jw
		         h_diff(i3) = dt2*sum_hdiff(i3)
		         
		         ! jw
		         h_flux_term(i3) = 0.5*dt2*sum_hlimiter(i3)
	         
	         
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
	         	dz_old = dz_cell(k,i) + (dz_cell_new(k,i) - dz_cell(k,i))*(subcycle)/Nt
	         	rhs(kk,i3) = area(i)*dz_old*con_cell(k,i,i3) &
	         	&				+ h_adv(i3) + h_diff(i3)

		         if(h_flux_limiter > 0) then
		         	rhs(kk,i3) = rhs(kk,i3) + h_flux_term(i3)
		         end if
		      end do ! jw
		   end do ! jw
		   ! jw
			
			
			! jw
			! jw
			! jw
			! jw
			if(num_vertical_layer == 1) then
				! jw
				! jw
			else
				do k=b_layer,t_layer
					kk = t_layer + 1 - k ! jw
					
					if(k == b_layer) then ! jw
						do i3=1,maxtran2
							! jw
							! jw
							! jw
							! jw
							! jw
							v_adv(i3) = dt2*( &
							&	&
							&	-abs(0.5*(Q_ik_theta(k,i) + abs(Q_ik_theta(k,i))))*con_cell(k  ,i,i3) &
							&	+abs(0.5*(Q_ik_theta(k,i) - abs(Q_ik_theta(k,i))))*con_cell(k+1,i,i3) &
							&	)
														
							rhs(kk,i3) = rhs(kk,i3) + v_adv(i3)
							
							if(v_flux_limiter > 0) then
								v_flux_term(i3) = 0.5*dt2*( &
								& &
								& -PSI_ik(k,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k,i) + abs(Q_ik_theta(k,i))))*(con_cell(k+1,i,i3)-con_cell(k  ,i,i3)) &
								& +PSI_ik(k,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k,i) - abs(Q_ik_theta(k,i))))*(con_cell(k  ,i,i3)-con_cell(k+1,i,i3)) &
								&  )
								
								rhs(kk,i3) = rhs(kk,i3) + v_flux_term(i3)
							end if
						end do
					else if(k == t_layer) then ! jw
						do i3=1,maxtran2
							! jw
							! jw
							! jw
							! jw
							v_adv(i3) = dt2*( &
							&	 abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))*con_cell(k-1,i,i3) &
							&	&
							&	&
							&	-abs(0.5*(Q_ik_theta(k-1,i) - abs(Q_ik_theta(k-1,i))))*con_cell(k  ,i,i3) )
														
							rhs(kk,i3) = rhs(kk,i3) + v_adv(i3)
							
							if(v_flux_limiter > 0) then
								v_flux_term(i3) = 0.5*dt2*( &
								&  PSI_ik(k-1,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))*(con_cell(k  ,i,i3)-con_cell(k-1,i,i3)) &
								&  &
								&  &
								& -PSI_ik(k-1,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k-1,i) - abs(Q_ik_theta(k-1,i))))*(con_cell(k-1,i,i3)-con_cell(k  ,i,i3)) )
								
								rhs(kk,i3) = rhs(kk,i3) + v_flux_term(i3)
							end if
						end do
					else ! jw
						do i3=1,maxtran2
							! jw
							! jw
							! jw
							! jw
							v_adv(i3) = dt2*( &
							&	 abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))*con_cell(k-1,i,i3) &
							&	-abs(0.5*(Q_ik_theta(k  ,i) + abs(Q_ik_theta(k  ,i))))*con_cell(k  ,i,i3) &
							&	+abs(0.5*(Q_ik_theta(k  ,i) - abs(Q_ik_theta(k  ,i))))*con_cell(k+1,i,i3) &
							&	-abs(0.5*(Q_ik_theta(k-1,i) - abs(Q_ik_theta(k-1,i))))*con_cell(k  ,i,i3) )
														
							rhs(kk,i3) = rhs(kk,i3) + v_adv(i3)
							
							if(v_flux_limiter > 0) then
								v_flux_term(i3) = 0.5*dt2*( &
								&	PSI_ik(k-1,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k-1,i) + abs(Q_ik_theta(k-1,i))))*(con_cell(k  ,i,i3)-con_cell(k-1,i,i3)) &
								& -PSI_ik(k  ,i,i3)*&
								&		abs(0.5*(Q_ik_theta(  k,i) + abs(Q_ik_theta(k  ,i))))*(con_cell(k+1,i,i3)-con_cell(k  ,i,i3)) &
								& +PSI_ik(k  ,i,i3)*&
								&		abs(0.5*(Q_ik_theta(  k,i) - abs(Q_ik_theta(k  ,i))))*(con_cell(k  ,i,i3)-con_cell(k+1,i,i3)) &
								& -PSI_ik(k-1,i,i3)*&
								&		abs(0.5*(Q_ik_theta(k-1,i) - abs(Q_ik_theta(k-1,i))))*(con_cell(k-1,i,i3)-con_cell(k  ,i,i3)) )
								
								rhs(kk,i3) = rhs(kk,i3) + v_flux_term(i3)	
							end if							
						end do
					end if ! jw
				end do ! jw
			end if ! jw
			! jw
			
	      			
			! jw
	      ! jw
	      ! jw
			! jw
			! jw
	      ! jw
	      call tridiagonal_solver(maxlayer, num_vertical_layer, maxtran2, a_lower_mat, b_diagonal_mat, c_upper_mat, rhs, solution, gam)
	      
	      do k = b_layer,t_layer
	      	kk = t_layer - k + 1 ! jw
	      	
	      	! jw
	      	! jw
	      	do i3=1,maxtran2
	      		if(tran_id(i3) == 1 .or. tran_id(i3) == 2) then ! jw
			      	if(solution(kk,i3) >= 0.0_dp .and. solution(kk,i3) <= 50.0_dp) then
			      		con_cell_new(k,i,i3) = solution(kk,i3)
			      	else if(solution(kk,i3) < 0.0_dp) then
			      		con_cell_new(k,i,i3) = 0.0_dp
			      	else if(solution(kk,i3) > 50.0_dp) then
			      		con_cell_new(k,i,i3) = 50.0_dp
			      	end if
			      end if
		      end do
	      end do
	      
	      ! jw
	      ! jw
	      if(SS_element_flag(i) > 0) then
	      	ith_SS = SS_element_flag(i)
	      	k = SS_layer(ith_SS)
	      	if(Q_add_SS(ith_SS) >= 0) then
	      		! jw
		      	do i3=1,maxtran2
			      	if(k == 999) then
			      		! jw
			      		con_cell_new(t_layer,i,i3) = con_cell_new(t_layer,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_at_SS(ith_SS,i3)*dt2/(area(i)*dz_cell_new(t_layer,i))
			      	else if(k == 0) then
			      		! jw
			      		con_cell_new(b_layer,i,i3) = con_cell_new(b_layer,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_at_SS(ith_SS,i3)*dt2/(area(i)*dz_cell_new(b_layer,i))
			      	else
			      		! jw
			      		con_cell_new(k,i,i3) = con_cell_new(k,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_at_SS(ith_SS,i3)*dt2/(area(i)*dz_cell_new(k,i))	
			      	end if
			      end do
			   else
			   	! jw
			   	do i3=1,maxtran2
			      	if(k == 999) then
			      		! jw
			      		con_cell_new(t_layer,i,i3) = con_cell_new(t_layer,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_cell(t_layer,i,i3)*dt2/(area(i)*dz_cell_new(t_layer,i))
			      	else if(k == 0) then
			      		! jw
			      		con_cell_new(b_layer,i,i3) = con_cell_new(b_layer,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_cell(b_layer,i,i3)*dt2/(area(i)*dz_cell_new(b_layer,i))
			      	else
			      		! jw
			      		con_cell_new(k,i,i3) = con_cell_new(k,i,i3) + &
			      		&	Q_add_SS(ith_SS)*con_cell(k,i,i3)*dt2/(area(i)*dz_cell_new(k,i))	
			      	end if
			      end do
		      end if
	      end if
	   end do ! jw
		!$omp end do
		! jw
	   
	   ! jw
	   ! jw
 	   !$omp do private(i,k,i3)
		do i=1,maxele
			do k=bottom_layer_at_element(i),top_layer_at_element(i)
		   	do i3=1,maxtran2
		   		con_cell(k,i,i3) = con_cell_new(k,i,i3)
		   	end do
		   end do
	 		
	 		! jw
	 		do k=0,bottom_layer_at_element(i)-1
	 			do i3=1,maxtran2
					con_cell(k,i,i3) = 0.0_dp
				end do
			end do
			do k=top_layer_at_element(i)+1,maxlayer+1
	 			do i3=1,maxtran2
					con_cell(k,i,i3) = 0.0_dp
				end do				
			end do
		end do
 		!$omp end do
 	   !$omp end parallel
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
	end do ! jw
	
	
	! jw
	! jw
	!$omp parallel
	!$omp do private(i,k,l,ii,i3,sum11,sum22,count1)
	do i=1,maxnod
		! jw
		! jw
		! jw
		! jw
		! jw
		do k=top_layer_at_node(i),bottom_layer_at_node(i),-1
			do i3=1,maxtran2
				sum11 = 0.0_dp ! jw
				sum22 = 0.0_dp ! jw
				count1 = 0 ! jw
				do l = 1,adj_cells_at_node(i)
					ii = adj_cellnum_at_node(l,i)
					
					! jw
					if(k>=bottom_layer_at_element(ii) .and. k<=top_layer_at_element(ii)) then ! jw
						sum11 = sum11 + area(ii)*dz_cell_new(k,ii) ! jw
						sum22 = sum22 + area(ii)*dz_cell_new(k,ii) * con_cell_new(k,ii,i3) ! jw
					else
						count1 = count1 + 1
					end if
				end do

				if(sum11 == 0.0_dp) then
					con_node_new(k,i,i3) = 0.0_dp
				else
					con_node_new(k,i,i3) = sum22/sum11 ! jw
				end if
				
				! jw
				! jw
				if(count1 == adj_cells_at_node(i)) then
					if(k /= top_layer_at_node(i)) then
						con_node_new(k,i,i3) = con_node_new(k+1,i,i3)
					end if
				end if				
			end do	
		end do
		
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
	!$omp end do

	! jw
	!$omp do private(j,n1,n2,k,i3)
	do j=1,maxface
		n1 = nodenum_at_face(1,j)
		n2 = nodenum_at_face(2,j)
		do k=bottom_layer_at_face(j),top_layer_at_face(j)
			do i3=1,maxtran2
				con_face_new(k,j,i3) = (con_node_new(k,n1,i3) + con_node_new(k,n2,i3)) * 0.5
			end do
		end do
	end do
	!$omp end do
	!$omp end parallel

	! jw
	do i3=1,maxtran2
		if(tran_id(i3) == 1) then ! jw
			salt_cell_new = con_cell_new(:,:,1)
			salt_node = con_node_new(:,:,1)
			salt_face = con_face_new(:,:,1)
		end if
		if(tran_id(i3) == 2) then ! jw
			temp_cell_new = con_cell_new(:,:,2)
			temp_node = con_node_new(:,:,2)
			temp_face = con_face_new(:,:,2)
		end if
	end do
	! jw

! jw
! jw
! jw
! jw
	! jw
end subroutine solve_transport_equation_v19