!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! 
!! This will find open boundary salinity and temperature condition for transport.
!! ===========================================================================!
subroutine find_openboundary_salt_temp_v3
	use mod_global_variables
	implicit none

	integer :: i, k, l, ii
	integer :: ie, j2, ii2
	! jw
	integer :: t_layer, b_layer

	! jw
	real(dp):: Q_jk_theta
	real(dp),dimension(num_ob_cell) :: salt_at_obc, temp_at_obc
! jw
	real(dp),dimension(maxlayer) :: sum1, sum2
	integer :: icount
	
	! jw
	integer :: t2
	real(dp):: u1,u2,u3,v1,v2,v3
	! jw

	
   ! jw
   ! jw
   ! jw
   ! jw
   ! jw
	! jw
	u2 = julian_day*86400.0 ! jw
	
	salt_at_obc = ref_salt
	temp_at_obc = ref_temp
	
	! jw
	! jw
	! jw
	! jw
	do i=1,num_ob_cell
		! jw
		! jw
		! jw
		! jw
		! jw
		
		
		! jw
		! jw
		! jw
		! jw
		ie = ob_cell_id(i) ! jw
		! jw
		! jw
		! jw
				
		! jw
		! jw
		if(is_salt == 1 .and. salt_ser_id(i) > 0) then 
			! jw
			ii = salt_ser_id(i)
			do t2=2,salt_ser_data_num(ii)
				! jw
				! jw
				! jw
				u1 = salt_ser_time(t2-1,ii) - reference_diff_days * 86400.0 	! jw
				u3 = salt_ser_time(t2  ,ii) - reference_diff_days * 86400.0		! jw
				
				! jw
				if(u2 >= u1 .and. u2 <= u3) then
					v1 = salt_ser_salt(t2-1,ii)		! jw
					v3 = salt_ser_salt(t2  ,ii)		! jw
					v2 = (u2-u1)*(v3-v1)/(u3-u1)+v1	! jw
					exit
				end if			
			end do
		else
			v2 = ref_salt ! jw
		end if
		! jw
		salt_at_obc(i) = spinup_function_baroclinic * v2
	
		! jw
		if(is_temp == 1 .and. temp_ser_id(i) > 0) then
			! jw
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
		else
			v2 = ref_temp ! jw
		end if		
		! jw
		temp_at_obc(i) = spinup_function_baroclinic * v2 ! jw
			
		do k=1,maxlayer
			salt_at_obck(k,i) = salt_at_obc(i)
			temp_at_obck(k,i) = temp_at_obc(i)
		end do

		! jw
		! jw
		! jw
		t_layer = top_layer_at_element(ie)
		b_layer = bottom_layer_at_element(ie)
		
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
		do k=1,b_layer-1
			salt_at_obck(k,i) = salt_at_obck(b_layer,i)
			temp_at_obck(k,i) = temp_at_obck(b_layer,i)
		end do
		do k=t_layer,maxlayer
			salt_at_obck(k,i) = salt_at_obck(t_layer,i)
			temp_at_obck(k,i) = temp_at_obck(t_layer,i)
		end do		
		
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
	! jw
	! jw
	do i=1,num_Qb_cell
		! jw
		! jw
		! jw
		! jw
		
		! jw
		if(is_salt == 1 .and. Q_salt_ser_id(i) > 0) then
			! jw
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
			salt_at_Qbc(i) = spinup_function_baroclinic * v2
		end if

		! jw
		if(is_temp == 1 .and. Q_temp_ser_id(i) > 0) then
			! jw
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
			temp_at_Qbc(i) = spinup_function_baroclinic * v2 ! jw
		end if
	end do ! jw
	! jw

end subroutine find_openboundary_salt_temp_v3