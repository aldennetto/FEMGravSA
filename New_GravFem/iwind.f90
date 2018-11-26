
!****************************************************************************
    integer function iwind(x,y,m,x0,y0)
!****************************************************************************
!   winding number algorithm for closed polygonal paths,
!	Kenneth Lelend, Mathematics of computation, vol. 29 (130),
!	554-558, 1975.  Determines whether the point x0,y0
!	falls outside, inside, or on the edge of the polygon 
!	defined by the array of points x,y 
!
!	for x,y counterclockwise:
!	iwind = 0  => outside
!	iwind = 1  => inside
!	iwind = 2  => on a polygonal edge
!	iwind = 3  => on a polygonal point
!
!	for x,y clockwise:
!	iwind =  0  => outside
!	iwind = -1  => inside
!	iwind =  2  => on a polygonal edge
!	iwind =  3  => on a polygonal point
!****************************************************************************
    implicit none
	integer :: m,mx,m1,s,k,t,up
	double precision :: zero,one,eps,det,x1,x2,y1,y2,dx,dy
	double precision :: x0,y0,x(m),y(m)
	zero=0.0
	eps=1.0d-5
	one=1.0
	s=0
	m1=m-1
	x2=x(1)-x0
	y2=y(1)-y0
	k=0
   10	k=k+1
		if(k>m1) goto 60
		up=0
		t=0
		x1=x2
		y1=y2
		x2=x(k+1)-x0
		y2=y(k+1)-y0
		if(x1==zero.and.y1==zero) goto 80
		if(y2==zero) goto 15
		if(y2>zero) goto 20
!	y2 < 0
		if(y1<zero) goto 30
!	y1 >/= 0, y2 < 0
		up=1
		goto 40
   15	if(x2==zero) goto 10
   20	If(y1>=zero) goto 30
!	y2 >/= 0, y1 < 0
		up=-1
		go to 40
!	end step 1; begin step 2
   30	if(y1/=zero) goto 10
		if(y2/=zero) goto 10
		if(x2==zero) goto 10
		if(sign(one,x1)/=sign(one,x2)) goto 70
		go to 10
!	end step 2; begin step 3
   40	continue
		dx=x2-x1
		dy=y2-y1
        det=(x2*y1-x1*y2)/(dx*dx+dy*dy)
		if(abs(det)<eps) goto 70
		if(float(up)*det>eps) goto 50
		t=-up
   50	continue
		s=s+t
	goto 10
   60	continue
	    iwind=-s
	return
   70   continue
    	iwind=2
	return
   80   continue
     	iwind=3
	return
	end function iwind




