!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!*************************************************************************************************
!*******************             check flux conservation             *****************************
!*******************   conservation of fluxes and total volume etc.  *****************************
subroutine check_conservation
   use mod_global_variables
	use mod_file_definition
   
   implicit none

   integer :: i, j, k, l, node, i_face_no
   integer :: ntmp, npflag
   real(dp):: total_volume, total_volume12, total_mass,   &
   &  	     total_potent_energy, total_kinetic_energy,  &
   &          zup, u_average, v_average, w_average

   real(dp) :: ftmp, xtmp, ytmp, wild, fluxbnd, fluxchan, fluxchan1, fluxchan2, fluxvert
	
	real(dp):: adv_depth(maxnod), adv_depth2(maxele)
	! jw
	
	! jw
   open(10,file='total.dat')
   total_volume=0 !total volume
   total_mass=0 !total mass
   total_potent_energy=0 !total potential energy
   total_kinetic_energy=0 !total kinetic energy
   do i = 1, maxele
      do k = bottom_layer_at_element(i), top_layer_at_element(i) ! jw
         total_volume=total_volume+area(i)*dz_cell(k,i)  
         total_mass=total_mass+rho_cell(k,i)*area(i)*dz_cell(k,i)  
         if(k==top_layer_at_element(i)) then
            zup=MSL+eta_cell_new(i)
         else
            zup=z_level(k)
         end if
         total_potent_energy = total_potent_energy   &
                 &           + 0.5*rho_cell(k,i)*gravity*area(i)*dz_cell(k,i)   &
                 &           * (2*zup-dz_cell(k,i))

         u_average=0 ! jw
         v_average=0 ! jw
         do l=1,tri_or_quad(i)
            node=nodenum_at_cell(l,i)
            u_average=u_average+(u_node(k,node)+u_node(k-1,node))/2/tri_or_quad(i)
            v_average=v_average+(v_node(k,node)+v_node(k-1,node))/2/tri_or_quad(i)
         end do
         w_average=(wn_cell(k,i)+wn_cell(k-1,i))/2 ! jw
         total_kinetic_energy = total_kinetic_energy   &
                  &           + rho_cell(k,i)*area(i)*dz_cell(k,i)   &
                  &           * (u_average**2+v_average**2+w_average**2)
      end do ! jw
   end do ! jw

   write(10,*) elapsed_time/86400.0,total_volume,total_mass,total_potent_energy,total_kinetic_energy

	! jw
   open(13,file='fluxflag.gr3')
   read(13,*)
   read(13,*)ntmp,npflag
   if(npflag/=maxnod) then
      write(pw_run_log,*)'# of pts in fluxflag should = maxnod'
      stop
   end if

   do i=1,maxnod
      read(13,*)j,xtmp,ytmp,wild
      adv_depth(i)=wild
   end do
   close(13)

   do i=1,maxele
      adv_depth2(i)=0
      do l=1,tri_or_quad(i)
         node=nodenum_at_cell(l,i)
         if(adv_depth2(i)<adv_depth(node)) adv_depth2(i)=adv_depth(node)
      end do
   end do

   open(9,file='flux.dat')
   open(8,file='salt.dat')
   total_volume12=0 !total volume inside rgns 1 and 2
   fluxbnd=0
   fluxchan=0  !total flux across unnatural bnds
   fluxchan1=0  !flux at bnd 1 (positive outward)
   fluxchan2=0  !flux at bnd 2 (positive outward)
   fluxvert=0
   do i=1,maxele
      if(adv_depth2(i) == 0) then
      	cycle
      end if
      if(top_layer_at_element(i)/=0) then
         fluxvert=fluxvert+wn_cell(top_layer_at_element(i),i)*area(i)
      end if
      do k=bottom_layer_at_element(i),top_layer_at_element(i) ! jw
         total_volume12=total_volume12+area(i)*dz_cell(k,i)  
      end do

      do l=1,tri_or_quad(i)
         i_face_no=facenum_at_cell(l,i)
         if(adj_cellnum_at_cell(l,i)==0) then ! jw
            do k=bottom_layer_at_face(i_face_no),top_layer_at_face(i_face_no) ! jw
               fluxbnd=fluxbnd+un_face(k,i_face_no)*face_length(i_face_no)*dz_face(k,i_face_no)
            end do
         else if(adv_depth2(adj_cellnum_at_cell(l,i))==0) then ! jw
            do k=bottom_layer_at_face(i_face_no),top_layer_at_face(i_face_no) ! jw
               ftmp=sign_in_outflow(l,i)*un_face(k,i_face_no)*face_length(i_face_no)*dz_face(k,i_face_no)
               fluxchan=fluxchan+ftmp
               if(adv_depth2(i)==1) then
                  fluxchan1=fluxchan1+ftmp
               else ! jw
                  fluxchan2=fluxchan2+ftmp
               end if
            end do
         end if ! jw
      end do ! jw
   end do ! jw

   write(9,*) elapsed_time/86400.0,fluxchan1,-fluxchan2,fluxbnd,total_volume12,fluxchan,&
     &       fluxvert,fluxbnd+fluxchan+fluxvert
   if(ishow == 1) then
   	write(*,*)'done computing fluxes...'
   end if
   write(pw_run_log,*) 'done computing fluxes...'
   ! jw

end subroutine check_conservation
