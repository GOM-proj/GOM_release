!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! 
!! subroutine for initialize heat model variables
!! jw, not yet finished
!! jw, this does nothing now
! 
subroutine initialize_transport_model
   use mod_global_variables
   implicit none
   
   ! integer :: i
   ! End of local variables ==================================================!
   
!   if(windp_flag == 2 .and. heat_model_flag /= 0) then
 !
 ! to eleminate heat for debug jl 3/22/04
 !
 !!!        call surf_fluxes(wtime1,wind_u1,wind_v1,   &
 !                   &       air_p1,air_temperature1,shum1,                 &
 !!                  &       srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz)
!      do i = 1, maxnod
!         sflux(i) = -fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i))
!      enddo

!      if(ishow == 1)then
!         write(*,*)'heat budge model completes...'
!      endif
!      write(520,*) 'heat budge model completes...'
!   endif

end subroutine initialize_transport_model
