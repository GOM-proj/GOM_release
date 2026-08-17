!! ===========================================================================! 
!! GOM is developed by Jungwoo Lee & Jun Lee
!! ===========================================================================! 
!! This subroutine is based on "utm_geo.f90" & 
!! Lat_Long_to_UTM_transformation.pdf
subroutine coordinate_conversion(x1,y1,UTM_PROJECTION_ZONE,option, x2,y2)
	implicit none

	! jw
	integer, parameter :: sp = kind(0.0) ! jw
	integer, parameter :: dp = kind(0.d0) ! jw
	! jw
	
	real(dp), intent(in) :: x1,y1
	integer , intent(in) :: UTM_PROJECTION_ZONE, option
	real(dp), intent(out):: x2,y2
	
	! jw
	real(dp), parameter :: semi_major_axis = 6378137.0_dp
	real(dp), parameter :: semi_minor_axis = 6356752.314245_dp

	! jw
	! jw
	real(dp), parameter :: scfa = 0.9996_dp ! jw
	real(dp), parameter :: north = 0.0_dp, east = 500000.0_dp
	
	real(dp), parameter :: pi = 3.141592653589793_dp
	real(dp), parameter :: deg2rad = pi/180.0_dp, rad2deg = 180.0_dp/pi
	logical :: lsouth	
	

	integer :: zone
	real(dp):: rlon,rlat
	real(dp):: e2,e4,e6,ep2,xx,yy,dlat,dlon,cm,cmr,delam
	real(dp):: f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
	! jw
	! jw
	
	! jw
	e2 = 1-(semi_minor_axis/semi_major_axis)**2 ! jw
	e4 = e2*e2
	e6 = e2*e4
	ep2=e2/(1.d0-e2)
	

   !
   !---- Set Zone parameters
   !
   lsouth = .false.   
   if(UTM_PROJECTION_ZONE < 0) then
   	lsouth = .true.
   end if
   zone = abs(UTM_PROJECTION_ZONE)
   cm = zone*6.0d0 - 183.d0
   cmr = cm*deg2rad
   
   
   if(option == 1) then ! jw
		dlon = x1 ! jw
		dlat = y1 ! jw
   else if(option == 2) then ! jw
   	xx = x1
     	yy = y1
     	if(lsouth) then
     		yy = yy - 1.d7
   	end if
   end if
   
   ! jw
   if(option == 1) then ! jw
      rlon = deg2rad*dlon ! jw
      rlat = deg2rad*dlat ! jw
      
      delam = dlon - cm
      if(delam < -180.d0) then
      	delam = delam + 360.d0
      end if
      if(delam > 180.d0) then
      	delam = delam - 360.d0
      end if
      delam = delam*deg2rad
      
      f1 = (1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256d0)*rlat
      f2 = 3.d0*e2/8.d0 + 3.d0*e4/32.d0 + 45.d0*e6/1024.d0
      f2 = f2*sin(2.d0*rlat)
      f3 = 15.d0*e4/256.d0*45.d0*e6/1024.d0
      f3 = f3*sin(4.d0*rlat)
      f4 = 35.d0*e6/3072.d0
      f4 = f4*sin(6.d0*rlat)
      rm = semi_major_axis*(f1 - f2 + f3 - f4)
      if(dlat == 90.d0 .or. dlat == -90.d0) then
      	xx = 0.d0
      	yy = scfa*rm
      else
         rn = semi_major_axis/sqrt(1.d0 - e2*sin(rlat)**2)
         t = tan(rlat)**2
         c = ep2*cos(rlat)**2
         a = cos(rlat)*delam
       
         f1 = (1.d0 - t + c)*a**3/6.d0
         f2 = 5.d0 - 18.d0*t + t**2 + 72.d0*c - 58.d0*ep2
         f2 = f2*a**5/120.d0
         xx = scfa*rn*(a + f1 + f2)
         f1 = a**2/2.d0
         f2 = 5.d0 - t + 9.d0*c + 4.d0*c**2
         f2 = f2*a**4/24.d0
         f3 = 61.d0 - 58.d0*t + t**2 + 600.d0*c - 330.d0*ep2
         f3 = f3*a**6/720.d0
         yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3))
      endif
      xx = xx + east
      yy = yy + north   
   
   ! jw
   else
      xx = xx - east
      yy = yy - north
      e1 = sqrt(1.d0 - e2)
      e1 = (1.d0 - e1)/(1.d0 + e1)
      rm = yy/scfa
      u = 1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256.d0
      u = rm/(semi_major_axis*u)
    
      f1 = 3.d0*e1/2.d0 - 27.d0*e1**3.d0/32.d0
      f1 = f1*sin(2.d0*u)
      f2 = 21.d0*e1**2/16.d0 - 55.d0*e1**4/32.d0
      f2 = f2*sin(4.d0*u)
      f3 = 151.d0*e1**3.d0/96.d0
      f3 = f3*sin(6.d0*u)
      rlat1 = u + f1 + f2 + f3
      dlat1 = rlat1*rad2deg
      if(dlat1 >= 90.d0 .or. dlat1 <= -90.d0) then
      	dlat1 = dmin1(dlat1,90.d0)
         dlat1 = dmax1(dlat1,-90.d0)
         dlon = cm
      else
         c1 = ep2*cos(rlat1)**2
         t1 = tan(rlat1)**2
         f1 = 1.d0 - e2*sin(rlat1)**2
         rn1 = semi_major_axis/sqrt(f1)
         r1 = semi_major_axis*(1.d0 - e2)/sqrt(f1**3)
         d = xx/(rn1*scfa)
     
         f1 = rn1*tan(rlat1)/r1
         f2 = d**2/2.d0
         f3 = 5.d0*3.d0*t1 + 10.d0*c1 - 4.d0*c1**2 - 9.d0*ep2
         f3 = f3*d**2*d**2/24.d0
         f4 = 61.d0 + 90.d0*t1 + 298.d0*c1 + 45.d0*t1**2 - 252.d0*ep2 - 3.d0*c1**2
         f4 = f4*(d**2)**3.d0/720.d0
         rlat = rlat1 - f1*(f2 - f3 + f4)
         dlat = rlat*rad2deg
     
         f1 = 1.d0 + 2.d0*t1 + c1
         f1 = f1*d**2*d/6.d0
         f2 = 5.d0 - 2.d0*c1 + 28.d0*t1 - 3.d0*c1**2 + 8.d0*ep2 + 24.d0*t1**2
         f2 = f2*(d**2)**2*d/120.d0
         rlon = cmr + (d - f1 + f2)/cos(rlat1)
         dlon = rlon*rad2deg
        	if (dlon < -180.d0) dlon = dlon + 360.d0
        	if (dlon > 180.d0) dlon = dlon - 360.d0
     	endif
   endif

  	if(option == 1) then ! jw
   	x2 = xx
   	if(lsouth) then
   		yy = yy + 1.d7
   	end if
    	y2 = yy
	else if(option == 2) then ! jw
   	x2 = dlon
    	y2 = dlat
   endif	
end subroutine coordinate_conversion