!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! initialization for the cold start alone
!! jw, not yet done, check wind options
subroutine setup_cold_start
   use mod_global_variables 
   implicit none

   ! jw
   ! jw
   
   ! jw
   eta_cell = 0.0
   eta_cell_new = 0.0
   wn_cell = 0.0
   
   eta_node = 0.0
   wetdry_node = 0

	if(i_transport_model_flag == 1)then
	   velo_u_transport = 0.0
      velo_v_transport = 0.0
      velo_w_transport = 0.0
   end if
	
	un_face = 0.0
   vn_face = 0.0

end subroutine setup_cold_start
