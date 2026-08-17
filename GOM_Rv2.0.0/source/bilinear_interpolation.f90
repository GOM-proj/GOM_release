!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! bi-linear interpolation for ELM
!! inverse bilinear mapping for quadrangles
!! convexity of the quad must have been checked, and (x,y) must not be outside the quad
!! ===========================================================================!
subroutine bilinear_interpolation(elem,x1,x2,x3,x4,y1,y2,y3,y4,x,y,xi,eta,shapef)
	use mod_global_variables
	use mod_file_definition	
	implicit none
	
	real(dp),parameter   :: small3=1.e-5
	real(dp),intent(in)  :: elem,x1,x2,x3,x4,y1,y2,y3,y4,x,y ! jw
	real(dp),intent(out) :: xi, eta, shapef(4)
	real(dp),dimension(2):: axi, aet, bxy, root_xi, root_et
	integer :: i, j, icaseno, icount
	real(dp):: x0, y0, dxi, deta, dd, beta, gamma, delta
	! jw

	! jw
	x0 = (x1+x2+x3+x4)/4.0
	y0 = (y1+y2+y3+y4)/4.0
	axi(1) = x2-x1+x3-x4
	axi(2) = y2-y1+y3-y4
	aet(1) = x3+x4-x1-x2
	aet(2) = y3+y4-y1-y2
	bxy(1) = x1-x2+x3-x4
	bxy(2) = y1-y2+y3-y4
	! jw
	! jw
	dxi  = 2.0*((x3-x4)*(y1-y2)-(y3-y4)*(x1-x2))
	deta = 2.0*((x4-x1)*(y3-y2)-(y4-y1)*(x3-x2))
	! jw
	if(dabs(bxy(1)) < small3 .and. dabs(bxy(2)) < small3  &
	&                        .or.                         &
	&  dabs(dxi)    < small3 .and. dabs(deta)<small3) then
		icaseno = 1      
		! jw
		dd = axi(1)*aet(2)-axi(2)*aet(1)
		if(dd == 0.0) then
			write(pw_run_log,*)'case 1 error:',dd
			stop
		end if
		
		xi  = 4.0*(aet(2)*(x-x0)-aet(1)*(y-y0))/dd
		eta = 4.0*(axi(1)*(y-y0)-axi(2)*(x-x0))/dd
   else if(dabs(dxi) < small3 .and. dabs(deta) >= small3) then   
		icaseno = 2      
		! jw
		eta = 4.0*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/deta
		dd  = (axi(1)+eta*bxy(1))**2+(axi(2)+eta*bxy(2))**2
		if(dd == 0.0) then
			write(pw_run_log,*)'case 2 error:',dd
			stop
		end if
		xi = (  (4.0*(x-x0)-eta*aet(1))*(axi(1)+eta*bxy(1))   &
		&     + (4.0*(y-y0)-eta*aet(2))*(axi(2)+eta*bxy(2)))/dd
   else if(dabs(dxi) >= small3 .and. dabs(deta) < small3) then   
		icaseno=3      
		! jw
		xi =  4.0*(bxy(2)*(x-x0)-bxy(1)*(y-y0))/dxi
		dd = (aet(1)+xi*bxy(1))**2+(aet(2)+xi*bxy(2))**2
		if(dd == 0.0) then
			write(pw_run_log,*)'case 3 error:',dd
			stop
		end if
		eta = ( (4.0*(x-x0)-xi*axi(1))*(aet(1)+xi*bxy(1))   &
		&   	+ (4.0*(y-y0)-xi*axi(2))*(aet(2)+xi*bxy(2)))/dd
	else ! jw
		icaseno = 4      
		! jw
		beta  = aet(2)*axi(1)-aet(1)*axi(2)-4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))
		gamma = 4.0*(aet(1)*(y-y0)-aet(2)*(x-x0))
		delta = beta*beta-4*gamma*dxi
		if(delta == 0.0) then
			xi  = -beta/2/dxi
			eta = (4.0*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-xi*dxi)/deta
		else if(delta>0) then
			! jw
			root_xi(1)=(-beta+dsqrt(delta))/2/dxi
			root_xi(2)=(-beta-dsqrt(delta))/2/dxi
			icount=0
			do i = 1, 2
				root_et(i)=(4*(bxy(2)*(x-x0)-bxy(1)*(y-y0))-root_xi(i)*dxi)/deta
				if(dabs(root_xi(i))<=1.1 .and. dabs(root_et(i))<=1.1) then
					xi=root_xi(i)
					eta=root_et(i)
					icount=icount+1
				end if
			end do
		
			if(icount==2 .and. dabs(root_xi(1)-root_xi(2))<small_06) then
				! jw
				! jw
				! jw
			else if(icount/=1) then
				write(pw_run_log,*)'abnormal instances',(root_xi(j),root_et(j),j=1,2),icount,elem
				write(pw_run_log,'(10e17.10)')x,y,x1,x2,x3,x4,y1,y2,y3,y4
				write(pw_run_log,*)dxi,deta,bxy(1),bxy(2)
				stop
			end if
      else
			write(pw_run_log,*)'no roots',delta,elem
			stop
		end if
	end if

	if(dabs(xi)>1.1 .or. dabs(eta)>1.1) then
		write(pw_run_log,*)'out of bound in ibilinear:',xi,eta,elem,icaseno
		write(pw_run_log,'(2e17.10)')x,y
		stop
	end if
	
	xi=dmin1(1.d0,dmax1(xi,-1.d0))
	eta=dmin1(1.d0,dmax1(eta,-1.d0))
	shapef(1)=(1-xi)*(1-eta)/4
	shapef(2)=(1+xi)*(1-eta)/4
	shapef(3)=(1+xi)*(1+eta)/4
	shapef(4)=(1-xi)*(1+eta)/4
end subroutine bilinear_interpolation
