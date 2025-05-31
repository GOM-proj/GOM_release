!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!*************************************************************************************************
!*****                             turbulence closure schemes                                *****
!*****  compute turbulence diffusivities vertical_eddy_viscosity, vertical_eddy_diffusivity, *****
!*****                              and in my-g, also qdiff.                                 *****
!*************************************************************************************************
subroutine turbulence
   use mod_global_variables
   use mod_file_definition
   
   implicit none
! jw

   integer :: i, j, k, l, n1, n2, klev, icount, ie, isd
   integer :: nqdim

   real(8) :: drhodz, bvf, dundz, dutdz, shear2, critical_richardson, vdiffc, tdiffc, vmax, vmin, st,   &
      &       udown, uup, vdown, vup, dudz, dvdz, rho_up, rho_down, etam, zctr, dists,    &
	  &       distb, q2fs, q2bot, xlfs, xlbot, eleft, eright, prod, buoy, diss,    &
	  &       pplus, pminus, psi, cpsi2p, cpsi3, tmp, q2l, upper, xl_max, en

   real(8),allocatable,dimension(:,:) :: richardson_number
   real(8),allocatable,dimension(:,:) :: q2bt, xlbt, rzbt, shearbt
   real(8),allocatable,dimension( : ) :: q2tmp, xltmp, xlmax


   real(8),allocatable,dimension( : ) :: a_lower_mat, b_diagonal_mat, c_upper_mat, gam
   real(8),allocatable,dimension(:,:) :: soln, rrhs

   real(8),allocatable,dimension(:,:) :: qdiff, qdiff2
   ! jw

!
! jw
!
      allocate(richardson_number(max_no_faces,maxlayer+1))

      allocate(q2bt(max_no_faces,maxlayer+1))
      allocate(xlbt(max_no_faces,maxlayer+1))
      allocate(rzbt(max_no_faces,maxlayer+1))
      allocate(shearbt(max_no_faces,maxlayer+1))

      allocate(q2tmp(maxlayer+1))
      allocate(xltmp(maxlayer+1))
      allocate(xlmax(maxlayer+1))

      allocate(qdiff(max_no_faces,maxlayer+1))
      allocate(qdiff2(max_no_faces,maxlayer+1))

!**************************************************************************************************
!*****             initialize vertical diffusivity
!*****             for ub & myg, compute vertical diffusivities at whole levels 
!**************************************************************************************************
      if(i_turbulence_flag == 3) then
         do i = 1, face_num
            do j = num_bottom_layer_at_face(i),num_top_layer_at_face(i) !wet sides
               call algebraic_stress(gravity,i,j)
               vertical_eddy_viscosity(i,j)   = dmin1(diffmax(i),dmax1(bgdiff,vd))
               vertical_eddy_diffusivity(i,j) = dmin1(diffmax(i),dmax1(bgdiff,td))
               qdiff(i,j)  = dmin1(diffmax(i),dmax1(bgdiff,qd))
               qdiff2(i,j) = dmin1(diffmax(i),dmax1(bgdiff,qd2))
            enddo !j=num_bottom_layer_at_face(i),num_top_layer_at_face(i)
         enddo !i=1,face_num
      endif !i_turbulence_flag=3
!**************************************************************************************************
!**************************************************************************************************

      do i = 1, face_num
         do k = num_bottom_layer_at_face(i), num_top_layer_at_face(i) !wet sides
            if(k == num_top_layer_at_face(i)) then 
               richardson_number(i,k) = 0.0
            else !k <= mj-1; implies m > m
               drhodz = (srho(i,k+1)-srho(i,k))/dzhalf(i,k)
               bvf    = -gravity*(drhodz/ref_water_density+gravity/1.5e3**2)
               dundz  = (velocity_normal(i,k+1)-velocity_normal(i,k))/dzhalf(i,k)
               dutdz  = (velocity_tangnt(i,k+1)-velocity_tangnt(i,k))/dzhalf(i,k)
               shear2 =  dundz**2+dutdz**2 !same form in 2 frames
               shear2=dmax1(shear2,1.0d-10)
               richardson_number(i,k)=dmax1(bvf/shear2,0.0d0)
            endif
         enddo !k=num_bottom_layer_at_face(i),num_top_layer_at_face(i)	
      enddo !i=1,face_num

!*************************************************************************************************
!... scheme 1: step function
!*************************************************************************************************
      if(i_turbulence_flag == 1) then
         critical_richardson = 0.25  !critical richardson #
         vdiffc = 1.0d-3 !sub-critical mixing coefficient
         tdiffc = 1.0*vdiffc
         do i = 1, face_num
            do k = num_bottom_layer_at_face(i), num_top_layer_at_face(i) !wet sides
               if(richardson_number(i,k) < critical_richardson) then
                  vertical_eddy_viscosity(i,k)   = vdiffc
                  vertical_eddy_diffusivity(i,k) = tdiffc
               else
                  vertical_eddy_viscosity(i,k)   = 1.0d-5
                  vertical_eddy_diffusivity(i,k) = 1.0d-5
               endif
            enddo !k
         enddo !i=1,face_num

         if(ishow == 1) write(*,*) 'done turbulence closure (step)...'
         write(pw_run_log,*) 'done turbulence closure (step)...'
      endif !i_turbulence_flag=1

!*************************************************************************************************
!... scheme 2: pacanowski and philander (1981)
!*************************************************************************************************
      if(i_turbulence_flag == 2) then
         do i = 1, face_num
            if(h_at_face_center(i) <= h1_pp) then
               vmax = vdmax_pp1
               vmin = vdmin_pp1
            elseif(h_at_face_center(i)<h2_pp) then
              vmax =  vdmax_pp1   &
			    &  + (vdmax_pp2-vdmax_pp1)*(h_at_face_center(i)-h1_pp)/(h2_pp-h1_pp)
              vmin =  vdmin_pp1   &
			    &  + (vdmin_pp2-vdmin_pp1)*(h_at_face_center(i)-h1_pp)/(h2_pp-h1_pp)
            else !h_at_face_center >= h2
               vmax = vdmax_pp2
               vmin = vdmin_pp2
            endif

            do k = num_bottom_layer_at_face(i),num_top_layer_at_face(i) !wet sides
!	    vmax >= vmin
               vertical_eddy_viscosity(i,k)   =  vmax/(1+5*richardson_number(i,k))**2+vmin
               vertical_eddy_diffusivity(i,k) =  vertical_eddy_viscosity(i,k)   &
			               &                  / (1.0 + 5.0*richardson_number(i,k))   &
						   &                  +  tdmin_pp
            enddo !k	
         enddo !i=1,face_num

         if(ishow == 1) write(*,*) 'done turbulence closure (pp)...'
         write(pw_run_log,*) 'done turbulence closure (pp)...'
      endif !i_turbulence_flag=2

!*************************************************************************************************
!... scheme 3: mellor-yamada-galperin & umlauf-burchard scheme
!*************************************************************************************************
      if(i_turbulence_flag == 3) then
!------------------------------------------------------------
         allocate(a_lower_mat(maxlayer+1))
         allocate(b_diagonal_mat(maxlayer+1))
         allocate(c_upper_mat(maxlayer+1))
         allocate(gam(maxlayer+1))

         allocate(soln(maxlayer+1,5), rrhs(maxlayer+1,5))
		 
         st = 0.0 !secnds(0.0) !timing the process
!...
!... solve for q2 and q2l (or psi in ub)
!...
! jw
         do j = 1, face_num
            if(num_top_layer_at_face(j) == 0) cycle
               n1 = nodenum_at_face(j,1)
               n2 = nodenum_at_face(j,2)
            do k = num_bottom_layer_at_face(j),num_top_layer_at_face(j)
               q2bt(j,k)=(q2(j,k-1)+q2(j,k))/2.0
               xlbt(j,k)=(xl(j,k-1)+xl(j,k))/2.0
               if(xlbt(j,k)<xlmin2(j)) then
                  write(pw_run_log,*)'xlbt < xlmin2',it,j,k,xlbt(j,k),xlmin2(j)
                  do l=num_bottom_layer_at_face(j)-1,num_top_layer_at_face(j)
                     write(pw_run_log,*)l,xlbt(j,l)
                  enddo 
                  stop
               endif
               udown=(u_at_node(n1,k-1)+u_at_node(n2,k-1))/2.0
               uup=(u_at_node(n1,k)+u_at_node(n2,k))/2.0
               vdown=(v_at_node(n1,k-1)+v_at_node(n2,k-1))/2.0
               vup=(v_at_node(n1,k)+v_at_node(n2,k))/2.0
               dudz=(uup-udown)/dz(j,k)
               dvdz=(vup-vdown)/dz(j,k)
               shearbt(j,k)=dudz**2+dvdz**2
               if(k==num_top_layer_at_face(j).or.k==num_bottom_layer_at_face(j)) then
                  drhodz=0
               else ! jw
                  rho_up=srho(j,k)+dz(j,k)/2/dzhalf(j,k)*(srho(j,k+1)-srho(j,k))
                  rho_down=srho(j,k-1)+dz(j,k-1)/2/dzhalf(j,k-1)*(srho(j,k)-srho(j,k-1))
                  drhodz=(rho_up-rho_down)/dz(j,k)
               endif
               rzbt(j,k)=gravity/ref_water_density*drhodz
            enddo !k=num_bottom_layer_at_face(j),num_top_layer_at_face(j)
         enddo !j=1,face_num


         do j = 1, face_num
            if(num_top_layer_at_face(j)==0) cycle

!	compute upper bound for xl (not used)
            if(adj_cellnum_at_face(j,2)/=0) then
               etam = dmax1(eta_old(adj_cellnum_at_face(j,1)),   &
                  &         eta_old(adj_cellnum_at_face(j,2)))
            else
               etam = eta_old(adj_cellnum_at_face(j,1))
            endif

            do k=num_bottom_layer_at_face(j),num_top_layer_at_face(j)
               if(k==num_top_layer_at_face(j)) then
                  zctr=z_mean_sea_level+etam-dz(j,k)/2
               else !m > m
                  zctr=z_level(k)-dz(j,k)/2
               endif
               dists=z_mean_sea_level+etam-zctr
               distb=zctr-(z_mean_sea_level-h_at_face_center(j))
!	  xlmax(k)=0.4*dmin1(dists,distb)
               xlmax(k)=0.4*dists*distb/(h_at_face_center(j)+etam)
               if(xlmax(k)<=0) then
                  write(pw_run_log,*)'dist<0 in my-g',j,k,h_at_face_center(j)+etam,dists,distb
                  stop
               endif
            enddo !k


!	b.c. (computed using values from previous time except wind)
            q2fs=16.6**(2.0/3)*dsqrt(wind_stress_normal(j)**2+wind_stress_tangnt(j)**2)/2
            q2fs=dmax1(q2fs,q2min)
            q2bot = 16.6**(2.0/3)*bottom_friction_coefficient(j)     &
             &    * dsqrt(velocity_normal(j,num_bottom_layer_at_face(j))**2.0d0   & 
             &          + velocity_tangnt(j,num_bottom_layer_at_face(j))**2) / 2.0d0
            q2bot=dmax1(q2bot,q2min)
            xlfs  = dmax1(xlmin1(j),dz(j,num_top_layer_at_face(j))/2*0.4)
            xlbot = dmax1(xlmin2(j),dz(j,num_bottom_layer_at_face(j))/2*0.4)

            nqdim=num_top_layer_at_face(j)-num_bottom_layer_at_face(j)+1 !>=1


!...    need to solve eq. system since nqdim >=1 (or m >= m)
!	matrix q
            do k=1,nqdim !row #
               klev=num_top_layer_at_face(j)-k+1 !level #; m <= klev <= m
               eleft  = 0.0 !for b_diagonal_mat
               eright = 0.0
               if(k>1) then
                  eleft = -qdiff(j,klev)*dt/dzhalf(j,klev)
                  a_lower_mat(k) = eleft
               endif
               if(k<nqdim) then
                  eright=-qdiff(j,klev-1)*dt/dzhalf(j,klev-1)
                  c_upper_mat(k)=eright
               endif

               prod = (vertical_eddy_viscosity(j,klev)+vertical_eddy_viscosity(j,klev-1))/2   &
                &   *  shearbt(j,klev)
               buoy = (vertical_eddy_diffusivity(j,klev)+vertical_eddy_diffusivity(j,klev-1))/2   &
                &   *  rzbt(j,klev)
               diss = cmiu0**3*dsqrt(q2bt(j,klev))/xlbt(j,klev) !diss/k
               if(prod+buoy>0) then
                  pplus=prod+buoy
                  pminus=0
               else
                  pplus=0
                  pminus=-prod-buoy
               endif
               b_diagonal_mat(k)=dz(j,klev)-eleft-eright+(diss+pminus/q2bt(j,klev))*dt*dz(j,klev)

!	  r.h.s. for q2
               rrhs(k,1)=dz(j,klev)*(q2bt(j,klev)+dt*pplus)
            enddo !k=1,nqdim

!	soln for q2 at new level
            call tridiagonal_solver(maxlayer+1, nqdim, 1,   &
                       &                a_lower_mat,b_diagonal_mat,c_upper_mat,rrhs,soln,gam)
            do k=1,nqdim
               klev=num_top_layer_at_face(j)-k+1 !level #; m <= klev <= m
               if(klev==num_top_layer_at_face(j)) then !fs value previals for m=m
                  q2tmp(klev)=q2fs
               elseif(klev==num_bottom_layer_at_face(j)) then
                  q2tmp(klev)=q2bot
               else
                  q2tmp(klev)=dmax1(soln(k,1),q2min)
               endif
            enddo !k


!	matrix ql
            do k=1,nqdim !row #
               klev  = num_top_layer_at_face(j)-k+1 !level #; m <= klev <= m
               eleft = 0
               eright=0
               if(k>1) then !klev <= m-1
                  eleft=-qdiff2(j,klev)*dt/dzhalf(j,klev)
                  a_lower_mat(k)=eleft
               endif
               if(k<nqdim) then !klev >= m+1
                  eright=-qdiff2(j,klev-1)*dt/dzhalf(j,klev-1)
                  c_upper_mat(k)=eright
               endif

!	  diagonal and rhs
               b_diagonal_mat(k)=dz(j,klev)-eleft-eright
               if(i_turbulence_model_name=='my'.or.i_turbulence_model_name=='kl') then
                  zctr=z_level(klev)-dz(j,klev)/2
                  dists=dmax1(z_mean_sea_level+etam-zctr,xlmin2(j))
                  distb=dmax1(zctr-(z_mean_sea_level-h_at_face_center(j)),xlmin2(j))
                  psi=1+1.33*(xlbt(j,klev)/0.4/distb)**2+0.25*(xlbt(j,klev)/0.4/dists)**2
                  cpsi2p=psi*cpsi2 !f_wall*cpsi2
               else !other gls
                  cpsi2p=cpsi2
               endif
               diss=cpsi2p*cmiu0**3*dsqrt(q2bt(j,klev))/xlbt(j,klev) !diss/k
               b_diagonal_mat(k)=b_diagonal_mat(k)+dt*dz(j,klev)*diss
  
               prod =  cpsi1   &
                &   * (vertical_eddy_viscosity(j,klev)+vertical_eddy_viscosity(j,klev-1))/2   &
                &   *  shearbt(j,klev)
               if(i_turbulence_model_name=='my') then
                  cpsi3=0.9
               else !gls models
                  if(rzbt(j,klev)>0) then !unstable
                     cpsi3=1
                  else !stable
                     select case(i_turbulence_model_name)
                        case('kl')
                           cpsi3=2.53
                        case('ke')
                           cpsi3=-0.52
                        case('kw')
                           cpsi3=-0.58
                        case('ub')
                           cpsi3=0.1
                        case default
                        write(pw_run_log,*)'unknown closure model:',i_turbulence_model_name
                        stop
                     end select
                  endif
               endif

               buoy = cpsi3   &
                 &  *(vertical_eddy_diffusivity(j,klev)+vertical_eddy_diffusivity(j,klev-1))/2   &
                 &  * rzbt(j,klev)
               if(prod+buoy>0) then
                  pplus=prod+buoy
                  pminus=0
               else
                  pplus=0
                  pminus=-prod-buoy
               endif
               b_diagonal_mat(k)=b_diagonal_mat(k)+pminus/q2bt(j,klev)*dt*dz(j,klev)
               if(klev == num_bottom_layer_at_face(j) .or.   &
               &  klev == num_top_layer_at_face(j))then
                  b_diagonal_mat(k)=b_diagonal_mat(k)+dt*0.4*rnub*qdiff2(j,klev)/xlbt(j,klev)
               endif
               tmp=cmiu0**rpub*q2bt(j,klev)**rmub*xlbt(j,klev)**rnub !psi^n_{j,k}
               rrhs(k,1)=dz(j,klev)*tmp*(1+dt*pplus/q2bt(j,klev))
            enddo !k=1,nqdim

!	soln for q2l and xl at new level
            call tridiagonal_solver(maxlayer+1,nqdim,1,a_lower_mat,b_diagonal_mat,c_upper_mat,rrhs,soln,gam)
            do k=1,nqdim
               klev=num_top_layer_at_face(j)-k+1
               q2l=dmax1(soln(k,1),psimin)
               if(klev==num_top_layer_at_face(j)) then
                  xltmp(klev)=xlfs
               elseif(klev==num_bottom_layer_at_face(j)) then
                  xltmp(klev)=xlbot
               else
                  xltmp(klev)=(q2l*cmiu0**(-rpub)*q2tmp(klev)**(-rmub))**(1/rnub)
               endif
!	  galperin's clipping 
               if(rzbt(j,klev)<0) then
                  upper=dsqrt(-0.56*q2tmp(klev)/rzbt(j,klev))
                  xltmp(klev)=dmin1(xltmp(klev),upper)
               endif
!	  max. length based on dissipation; xlmin2 prevails
               xl_max=(cmiu0*dsqrt(q2tmp(klev)))**3/eps_min
               xltmp(klev)=dmax1(xlmin2(j),dmin1(xl_max,xltmp(klev)))
            enddo ! jw


!	convert to whole levels
            do k=num_bottom_layer_at_face(j)-1,num_top_layer_at_face(j)
               if(k==num_bottom_layer_at_face(j)-1) then
                  q2(j,k)=q2bot 
                  xl(j,k)=xlbot
               elseif(k==num_top_layer_at_face(j)) then
                  q2(j,k)=q2fs
                  xl(j,k)=xlfs
               else !m > m and m<= k <= m-1
                  q2(j,k)=q2tmp(k)+dz(j,k)/2/dzhalf(j,k)*(q2tmp(k+1)-q2tmp(k))
                  xl(j,k)=xltmp(k)+dz(j,k)/2/dzhalf(j,k)*(xltmp(k+1)-xltmp(k))
                  if(q2(j,k)<0.or.xl(j,k)<xlmin2(j)) then
                     write(pw_run_log,*)'negative q2/xl',q2(j,k),xl(j,k)
                     stop
                  endif
               endif
            enddo !k=num_bottom_layer_at_face(j)-1,num_top_layer_at_face(j)

         enddo !j=1,face_num

!...
!... compute vertical diffusivities at new time
!...
         do i=1,face_num
            do j=num_bottom_layer_at_face(i),num_top_layer_at_face(i) !wet sides
               call algebraic_stress(gravity,i,j)
               vertical_eddy_viscosity(i,j)   = dmin1(diffmax(i),dmax1(bgdiff,vd))
               vertical_eddy_diffusivity(i,j) = dmin1(diffmax(i),dmax1(bgdiff,td))
               qdiff(i,j)=dmin1(diffmax(i),dmax1(bgdiff,qd))
               qdiff2(i,j)=dmin1(diffmax(i),dmax1(bgdiff,qd2))
            enddo !j=num_bottom_layer_at_face(i),num_top_layer_at_face(i)
         enddo !i=1,face_num

!...  extend q2 and xl etc. to account for level changes next time
         do i=1,face_num
            if(num_top_layer_at_face(i)==0) cycle
            do j=0,num_bottom_layer_at_face(i)-2
               q2(i,j)=q2(i,num_bottom_layer_at_face(i)-1)
               xl(i,j)=xl(i,num_bottom_layer_at_face(i)-1)
            enddo 
            do j=1,num_bottom_layer_at_face(i)-1
               vertical_eddy_diffusivity(i,j) = vertical_eddy_diffusivity(i,num_bottom_layer_at_face(i))
               vertical_eddy_viscosity(i,j)   = vertical_eddy_viscosity(i,num_bottom_layer_at_face(i))
               qdiff(i,j)=qdiff(i,num_bottom_layer_at_face(i))
               qdiff2(i,j)=qdiff2(i,num_bottom_layer_at_face(i))
            enddo !j=1,num_bottom_layer_at_face(i)-1

            do j=num_top_layer_at_face(i)+1,maxlayer
               q2(i,j)=q2(i,num_top_layer_at_face(i))
               xl(i,j)=xl(i,num_top_layer_at_face(i))
               vertical_eddy_diffusivity(i,j) = vertical_eddy_diffusivity(i,num_top_layer_at_face(i))
               vertical_eddy_viscosity(i,j)   = vertical_eddy_viscosity(i,num_top_layer_at_face(i))
               qdiff(i,j)=qdiff(i,num_top_layer_at_face(i))
               qdiff2(i,j)=qdiff2(i,num_top_layer_at_face(i))
            enddo !j=num_top_layer_at_face(i)+1,maxlayer
         enddo !i


         deallocate(a_lower_mat,b_diagonal_mat,c_upper_mat,gam)
         deallocate(soln, rrhs)

         en=0 !secnds(st)
         if(ishow == 1) write(*,*)'myg-ub took',en,'seconds'
         write(pw_run_log,*) 'myg-ub took',en,'seconds'
      endif !i_turbulence_flag=3
!*************************************************************************************************
!*****                               end of my-g scheme
!*************************************************************************************************

!*************************************************************************************************
!*****      compute diffusivity at point and whole level: tdiffp for transport equation      *****
!*************************************************************************************************
      if(baroclinic_flag == 0 .or. barotropic_flag == 1)then
         do i = 1, maxnod
            do k = i_bottom_index_of_node(i), i_surface_index_of_node(i) !wet nodes only
               tdiffp(i,k) = 0.0d0
               icount = 0
               do j = 1, adj_cells_at_node(i)
                  ie = adj_cellnum_at_node(i,j) 
                  do l = 1, tri_or_quad(ie)
                     isd = facenum_at_cell(ie,l)
                     if(num_top_layer_at_face(isd) /= 0) then
                        icount = icount + 1
                        klev = min(max(k,num_bottom_layer_at_face(isd)),num_top_layer_at_face(isd))
                        tdiffp(i,k) = tdiffp(i,k)+vertical_eddy_diffusivity(isd,klev)
                     endif
                  enddo  ! jw
               enddo  ! jw

               if(icount == 0) then
                  tdiffp(i,k) = 1.0e-2
               else
                  tdiffp(i,k) = tdiffp(i,k)/icount
               endif
            enddo !k=i_bottom_index_of_node(i),i_surface_index_of_node(i)
!...	extend
!
!******  this is original.
!******  this loop changed, because in 2d, the i_surface_index_of_node(i) = 0 
! jw
! jw
! jw
!******  changed like this
            if(i_surface_index_of_node(i) /= 0) then
               do k = i_surface_index_of_node(i)+1, maxlayer
                  tdiffp(i,k) = tdiffp(i,i_surface_index_of_node(i))
               enddo !k
            endif
         enddo !i=1,maxnod
      endif


      deallocate(richardson_number)
      deallocate(q2bt,xlbt,rzbt,shearbt)
      deallocate(q2tmp,xltmp,xlmax)

      deallocate(qdiff, qdiff2)

end subroutine turbulence


