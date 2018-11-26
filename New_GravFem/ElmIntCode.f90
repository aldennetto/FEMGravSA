!--------------------------------------------------------------------- 
  subroutine SpecGrav3d(ielv,dgr,iobs,n)
!---------------------------------------------------------------------
! Computes the volume element stiffness matrix for the 3D DC electrical
! resistivty problem for a unit local conductivity. To complete the
! stiffness calculation, each term in the the stiffness matrix must be
! scaled by the conductivity of the element.
!---------------------------------------------------------------------
  use global_grav
  implicit none
  integer :: ielv,iobs
  integer :: nn,ni,i,j,k,ip,n
  double precision :: node(8,3),xi(8,3),wt(8)
  double precision :: xip,yip,zip,dist3
  double precision :: dpsix(3,8),dsdx(3,3),detj
  double precision :: psi(8),dpsi(8,3)
  double precision :: dgr,dg

  !Initialize element stiffness matrix
  dgr=0.0d0 ! gravitational contribution due to the current element

  !Get the gauss integration points for the current element
  nn=8      !nndelv(ielv) = list of number of nodes in each element
  call gaus3d(nn,ni,xi,wt)  !calculates the location of the gauss integration points in local coordinates

        !xi = location
        !wt = weight
    ! nodelv = global list of nodes associated with each element for all elements

  !Recover the node coordinates and conductivities
  if(n==1) then
    do i=1,nn
      node(i,1)=element(ielv)%x(i)
      node(i,2)=element(ielv)%y(i)
      node(i,3)=element(ielv)%z(i)
    end do
  else
    do i=1,nn
      node(i,1)=elemento(ielv)%x(i)
      node(i,2)=elemento(ielv)%y(i)
      node(i,3)=elemento(ielv)%z(i)
    end do
  endif
  !Began integration of the sub-stiffness matrices
  do ip=1,ni
    !compute the value of the trilinear shape function at this integration point
    call shape3d(nn,xi(ip,1),xi(ip,2),xi(ip,3),psi,dpsi)

    !Calculate dxds and detj
    call trnfrm3d(nn,node,dpsi,dsdx,detj)
    if(detj.le.0.0) then
      write(*,*) 'Bad Jacobian in volume element ',ielv,element(ielv)%id,element(ielv)%type
      write(*,*) 'Program terminated.'
      write(*,*) 'Press enter to exit the program...'
      read(*,*)
!      close(outunt)
      stop
    end if


    xip=0.0d0
    yip=0.0d0
    zip=0.0d0
    do i=1,8
      xip=xip+node(i,1)*psi(i)
      yip=yip+node(i,2)*psi(i)
      zip=zip+node(i,3)*psi(i)
    end do


    !Accumulate integration point value of integrals
    dist3=(sqrt((observation(iobs)%x-xip)**2+(observation(iobs)%y-yip)**2&
              +(observation(iobs)%z-zip)**2))**3

    dg=((observation(iobs)%x-xip)*observation(iobs)%ux)+((observation(iobs)%y-yip)&
        *observation(iobs)%uy)+((observation(iobs)%z-zip)*observation(iobs)%uz)

    dgr=dgr+(ugc*detj*wt(ip)*dg/dist3)

!    if(dgr<0) stop

  end do


  return
  end subroutine SpecGrav3d

!---------------------------------------------------------------------

  subroutine gaus3d(nn,ni,xi,wt)
!---------------------------------------------------------------------
! set quadrature points for 3D elements
!---------------------------------------------------------------------
  implicit none
  integer i,nn,ni
  double precision xi(8,3),wt(8),rt1o3,rt3,a,b

  if(nn==8) then
    !set eight-point quadrature for 8-node trilinear brick elements
    !Set the gauss point weights
    ni=8
    wt=1.0d0

    !set gauss point coordinates
    rt1o3=dsqrt(1.0d0/3.0d0)
    xi(1,1)=rt1o3
    xi(1,2)=rt1o3
	xi(1,3)=rt1o3

	xi(2,1)=-rt1o3
	xi(2,2)=rt1o3
	xi(2,3)=rt1o3

	xi(3,1)=-rt1o3
	xi(3,2)=-rt1o3
	xi(3,3)=rt1o3

	xi(4,1)=rt1o3
	xi(4,2)=-rt1o3
	xi(4,3)=rt1o3

	xi(5,1)=rt1o3
	xi(5,2)=rt1o3
	xi(5,3)=-rt1o3

	xi(6,1)=-rt1o3
	xi(6,2)=rt1o3
	xi(6,3)=-rt1o3

	xi(7,1)=-rt1o3
	xi(7,2)=-rt1o3
	xi(7,3)=-rt1o3

	xi(8,1)=rt1o3
	xi(8,2)=-rt1o3
	xi(8,3)=-rt1o3

  else if(nn==6) then
    !set 6-point quardrature for 6-node trilinear wedge elements
    ni=6
    rt1o3=1.0d0/3.0d0*dsqrt(3.0d0)
    wt=1.0d0/6.0d0 !given as 1/3 in Dhatt & Touzot
    xi(1,1)=1.0d0/6.0d0
    xi(1,2)=2.0d0/3.0d0
    xi(1,3)=-rt1o3

    xi(2,1)=1.0d0/6.0d0
    xi(2,2)=1.0d0/6.0d0
    xi(2,3)=-rt1o3

    xi(3,1)=2.0d0/3.0d0
    xi(3,2)=1.0d0/6.0d0
    xi(3,3)=-rt1o3

    xi(4,1)=1.0d0/6.0d0
    xi(4,2)=2.0d0/3.0d0
    xi(4,3)=rt1o3

    xi(5,1)=1.0d0/6.0d0
    xi(5,2)=1.0d0/6.0d0
    xi(5,3)=rt1o3

    xi(6,1)=2.0d0/3.0d0
    xi(6,2)=1.0d0/6.0d0
    xi(6,3)=rt1o3

  else if(nn==4) then
    !set one-point quadrature for 4-node trilinear tetrahedral elements
    ni=4
    wt=0.25d0
    xi(1,1)=0.5854102d0
    xi(1,2)=0.1381966d0
    xi(1,3)=0.1381966d0

    xi(2,1)=0.1381966d0
    xi(2,2)=0.5854102d0
    xi(2,3)=0.1381966d0

    xi(3,1)=0.1381966d0
    xi(3,2)=0.1381966d0
    xi(3,3)=0.5854102d0

    xi(4,1)=0.1381966d0
    xi(4,2)=0.1381966d0
    xi(4,3)=0.1381966d0
  else
!    write(outunt,*) 'Element integration rule not recognized in routine quas3d.'
!    write(outunt,*) 'Program terminated.'
!    close(outunt)
    stop
  end if

  return
  end subroutine gaus3d

!---------------------------------------------------------------------
  subroutine shape3d(nn,eps,eta,zet,psi,dpsi)
!---------------------------------------------------------------------
! Compute shape functions and their derivatives
! for 3D body elements at the point eps,eta,zet
!---------------------------------------------------------------------
  implicit none
  integer nn
  double precision eps,eta,zet,lamda,psi(8),dpsi(8,3)
  double precision a1,a2,b1,b2,c1,c2,d
  if(nn==8) then
    !trilinear shape functions for a hexahedron element base on
    !Dhatt and Touzot, 1984, section 2.6.1, page 114
    a1=1.0d0+eps
    a2=1.0d0-eps
    b1=1.0d0+eta
    b2=1.0d0-eta
    c1=1.0d0+zet
    c2=1.0d0-zet
    d=0.125d0

    !Shape functions
    psi(1)=a2*b2*c2*d
    psi(2)=a1*b2*c2*d
    psi(3)=a1*b1*c2*d
    psi(4)=a2*b1*c2*d
    psi(5)=a2*b2*c1*d
    psi(6)=a1*b2*c1*d
    psi(7)=a1*b1*c1*d
    psi(8)=a2*b1*c1*d

    !Shape function derivatives with respect to xi(1)
    dpsi(1,1)=-b2*c2*d
    dpsi(2,1)= b2*c2*d
    dpsi(3,1)= b1*c2*d
    dpsi(4,1)=-b1*c2*d
    dpsi(5,1)=-b2*c1*d
    dpsi(6,1)= b2*c1*d
    dpsi(7,1)= b1*c1*d
    dpsi(8,1)=-b1*c1*d

    !Shape function derivatives with respect to xi(2)
    dpsi(1,2)=-a2*c2*d
    dpsi(2,2)=-a1*c2*d
    dpsi(3,2)= a1*c2*d
    dpsi(4,2)= a2*c2*d
    dpsi(5,2)=-a2*c1*d
    dpsi(6,2)=-a1*c1*d
    dpsi(7,2)= a1*c1*d
    dpsi(8,2)= a2*c1*d

    !Shape function derivatives with respect to xi(3)
    dpsi(1,3)=-a2*b2*d
    dpsi(2,3)=-a1*b2*d
    dpsi(3,3)=-a1*b1*d
    dpsi(4,3)=-a2*b1*d
    dpsi(5,3)= a2*b2*d
    dpsi(6,3)= a1*b2*d
    dpsi(7,3)= a1*b1*d
    dpsi(8,3)= a2*b1*d

  else if(nn==6) then
    !trilinear shape functions for a six-node wdge element base on
    !Dhatt and Touzot, 1984, section 2.7.1, page 120
    lamda=1.0d0-eps-eta
    a1=0.5d0*(1.0d0-zet)
    b1=0.5d0*(1.0d0+zet)
    psi(1)=lamda*a1
    psi(2)=eps*a1
    psi(3)=eta*a1
    psi(4)=lamda*b1
    psi(5)=eps*b1
    psi(6)=eta*b1

    !Shape function derivatives with respect to xi(1)
    dpsi(1,1)=-a1
    dpsi(2,1)=a1
    dpsi(3,1)=0.0d0
    dpsi(4,1)=-b1
    dpsi(5,1)=b1
    dpsi(6,1)=0.0d0

    !Shape function derivatives with respect to xi(2)
    dpsi(1,2)=-a1
    dpsi(2,2)=0.0d0
    dpsi(3,2)=a1
    dpsi(4,2)=-b1
    dpsi(5,2)=0.0d0
    dpsi(6,2)=b1

    !Shape function derivatives with respect to xi(3)
    dpsi(1,3)=-0.5d0*lamda
    dpsi(2,3)=-0.5d0*eps
    dpsi(3,3)=-0.5d0*eta
    dpsi(4,3)=0.5d0*lamda
    dpsi(5,3)=0.5d0*eps
    dpsi(6,3)=0.5d0*eta

  else if(nn==4) then
    !trilinear shape functions for a tetrahedral element base on
    !Dhatt and Touzot, 1984, section 2.5.2, page 111
    psi(1)=1.0d0-eps-eta-zet
    psi(2)=eps
    psi(3)=eta
    psi(4)=zet

    !Shape function derivatives with respect to xi(1)
    dpsi(1,1)=-1.0d0
    dpsi(2,1)=1.0d0
    dpsi(3,1)=0.0d0
    dpsi(4,1)=0.0d0

    !Shape function derivatives with respect to xi(2)
    dpsi(1,2)=-1.0d0
    dpsi(2,2)=0.0d0
    dpsi(3,2)=1.0d0
    dpsi(4,2)=0.0d0

    !Shape function derivatives with respect to xi(3)
    dpsi(1,3)=-1.0d0
    dpsi(2,3)=0.0d0
    dpsi(3,3)=0.0d0
    dpsi(4,3)=1.0d0
  else
!    write(outunt,*) 'Error:  element type not recognized in routine shape3d.'
!    write(outunt,*) 'Program terminated.'
!    close(outunt)
    stop
  end if

  return
  end subroutine shape3d

!---------------------------------------------------------------------
  subroutine trnfrm3d(n,x,dpsi,dsdx,detj)
!---------------------------------------------------------------------
! Compute the transformation matrix and its determinent in a 3D
! element at the local point for which dpsi has been computed.
!---------------------------------------------------------------------
  implicit none
  integer i,j,k,n
  double precision x(8,3),dpsi(8,3)
  double precision dsdx(3,3),detj,dxds(3,3)
  double precision  detji,sum

  !Calculate the jacobian matrix of derivatives of the global
  !nodal coordinates x, relative to local coordinates
  do i=1,3
    do j=1,3
      sum=0.0d0
	  do k=1,n
        sum=sum+dpsi(k,j)*x(k,i)
      end do
      dxds(i,j)=sum
    end do
  end do

  !Calculate the determinant of the jacobian matrix
  detj=dxds(1,1)*(dxds(2,2)*dxds(3,3)-dxds(2,3)*dxds(3,2)) &
    & -dxds(1,2)*(dxds(2,1)*dxds(3,3)-dxds(2,3)*dxds(3,1)) &
    & +dxds(1,3)*(dxds(2,1)*dxds(3,2)-dxds(2,2)*dxds(3,1))

  !Compute the transformation matrix dsdx
  detji=1.0d0/detj
  dsdx(1,1)=(dxds(2,2)*dxds(3,3)-dxds(3,2)*dxds(2,3))*detji
  dsdx(1,2)=(dxds(1,3)*dxds(3,2)-dxds(1,2)*dxds(3,3))*detji
  dsdx(1,3)=(dxds(1,2)*dxds(2,3)-dxds(1,3)*dxds(2,2))*detji
  dsdx(2,1)=(dxds(3,1)*dxds(2,3)-dxds(2,1)*dxds(3,3))*detji
  dsdx(2,2)=(dxds(1,1)*dxds(3,3)-dxds(1,3)*dxds(3,1))*detji
  dsdx(2,3)=(dxds(2,1)*dxds(1,3)-dxds(2,3)*dxds(1,1))*detji
  dsdx(3,1)=(dxds(2,1)*dxds(3,2)-dxds(3,1)*dxds(2,2))*detji
  dsdx(3,2)=(dxds(1,2)*dxds(3,1)-dxds(3,2)*dxds(1,1))*detji
  dsdx(3,3)=(dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1))*detji

  return
  end subroutine trnfrm3d
!---------------------------------------------------------------------
! Initialize grid

  subroutine initialize_parameter_grids
  use global_grav
  implicit none

  integer           :: n,err,a,b,i
  double precision  :: rdum1,rdum2,latitude,longitude,depth,den
  character(len=17) :: type
  character         :: str1,line*200

  open(unit=10,file='GravModelv8.txt',status='old')
  read(10,*) nlat,nlong
  n_paragrids=nlat*nlong*7
  allocate(paragrid(n_paragrids))

  do i=1,n_paragrids
    paragrid(i)%moho=.false.
  enddo
  n=1

  do
    read(10,'(a)',iostat=err) line

    if(err.eq.0) then
      if(line(1:1).eq.'#') then
        read(line,*) str1,longitude,latitude
        paragrid(n)%lat = latitude*pi/180.0d0
        paragrid(n)%long = (360.0d0+longitude)*pi/180.0d0
        a=1
        b=1
      else
        if(a==1) then
           paragrid(n)%h=-100.0d0
           paragrid(n)%layer='As'
           n=n+1
           a=a+1
        endif
        read(line,*)depth,den,rdum1,rdum2,type
        depth=(depth/1000.0d0)*-1
        if ((depth.ge.refmoho) .and. (b==1)) then
           paragrid(n)%lat=paragrid(n-1)%lat
           paragrid(n)%long=paragrid(n-1)%long
           paragrid(n)%h=refmoho
           paragrid(n)%moho=.true.
           paragrid(n)%density=den
           a=a+1
           n=n+1
           b=99
           if(a==3) then
             paragrid(n-1)%layer='Ma'
           elseif(a==4) then
             paragrid(n-1)%layer='Cr'
           else
             stop
             print *, 'Incorrect Moho'
           endif
        endif
        select case (type)
            case ('Lid')
                paragrid(n)%lat=paragrid(n-1)%lat
                paragrid(n)%long=paragrid(n-1)%long
                paragrid(n)%h=depth
                paragrid(n)%layer='Ma'
                paragrid(n)%density=den
                a=a+1
                n=n+1
            case ('Crust')
                paragrid(n)%lat=paragrid(n-1)%lat
                paragrid(n)%long=paragrid(n-1)%long
                paragrid(n)%h=depth
                paragrid(n)%density=den
                paragrid(n)%layer='Cr'
                a=a+1
                n=n+1
            case ('Seds3')
                paragrid(n)%lat=paragrid(n-1)%lat
                paragrid(n)%long=paragrid(n-1)%long
                paragrid(n)%h=depth
                paragrid(n)%density=den
                paragrid(n)%layer='S3'
                a=a+1
                n=n+1
            case ('Seds2')
                paragrid(n)%lat=paragrid(n-1)%lat
                paragrid(n)%long=paragrid(n-1)%long
                paragrid(n)%h=depth
                paragrid(n)%density=den
                paragrid(n)%layer='S2'
                a=a+1
                n=n+1
            case ('Seds1')
                paragrid(n)%lat=paragrid(n-1)%lat
                paragrid(n)%long=paragrid(n-1)%long
                paragrid(n)%h=depth
                paragrid(n)%density=den
                paragrid(n)%layer='S1'
                n=n+1
        end select

      endif

    else
      exit
    endif
  end do

  close(10)

  end subroutine initialize_parameter_grids

!---------------------------------------------------------------------
  subroutine initialize_cartesian_grids

    use global_grav
    implicit none

    integer         :: n,i,j,k
    double precision   :: a,b,nphi

    a=6378.1370d0
    b=6356.7523d0
    n=1
    allocate(x(nlat,nlong,7))
    allocate(y(nlat,nlong,7))
    allocate(z(nlat,nlong,7))
    allocate(density(nlat,nlong,7))
    allocate(layer_type(nlat,nlong,7))
    allocate(MohoP(nlat,nlong,7))

    do i=1,nlat
       do j=1,nlong
          do k=1,7
            nphi=(a**2)/sqrt(((a**2)*(cos(paragrid(n)%lat))**2)+((b**2)*&
                    (sin(paragrid(n)%lat))**2))
            x(i,j,k)=(nphi+paragrid(n)%h)*cos(paragrid(n)%lat)*&
                                  cos(paragrid(n)%long)
            y(i,j,k)=(nphi+paragrid(n)%h)*cos(paragrid(n)%lat)*&
                                  sin(paragrid(n)%long)
            z(i,j,k)=((b**2/a**2)*nphi+paragrid(n)%h)*sin(paragrid(n)%lat)
            density(i,j,k)=paragrid(n)%density
            layer_type(i,j,k)=paragrid(n)%layer

            if(paragrid(n)%moho==.true.) then
              MohoP(i,j,k)=.true.
            else
              MohoP(i,j,k)=.false.
            endif

            n=n+1
          end do
       end do
    end do

  deallocate(paragrid)


  end subroutine initialize_cartesian_grids

!========================================================================
  subroutine initialize_elements
!========================================================================

    use global_grav
    implicit none
    integer          :: i,j,k,n,e1,e,a,m,c1,c2,saltele
    double precision :: t1,t2,t3,t4,t5,l1,l2,l3
    character(len=2) :: type1
    integer,dimension(:),allocatable :: salt

    open(unit=10,file='SaltElements.dat',status='old')
    read(10,*) saltele
    allocate(salt(saltele))

    do i=1,saltele
      read(10,*) salt(i)
    enddo
    close(10)

    c1=0
    c2=0

    e1=1
    e=1
    a=0

    do i=1,nlat-1
       do j=1,nlong-1
          do k=1,6
              if((layer_type(i+1,j,k+1)==layer_type(i,j,k+1)).and.(layer_type(i,j+1,k+1)==&
                  layer_type(i+1,j+1,k+1)).and.(layer_type(i,j,k+1)==layer_type(i,j+1,k+1))) then

                    type1=layer_type(i+1,j,k+1)

              else
                    type1='CE'        ! Combination element
              endif
              if((type1=='Cr').or.(type1=='CE')) then
                c1=c1+1
              else
                c2=c2+1
              endif
          enddo
       enddo
    enddo

    n_elements=c2
    c_elmnts=c1
    allocate(element(c_elmnts))
    allocate(elemento(n_elements))

    ncol=(nlong-1)*(nlat-1)
    allocate(cr(ncol))

    do i=1,ncol
      cr(i)%double=.false.
    enddo

!    allocate(object(nlat*nlong))
    do i=1,nlat-1
       do j=1,nlong-1
          a=a+1
          do k=1,6

! Split elements between crustal elements and sediment/mantle elements
              if((layer_type(i+1,j,k+1)==layer_type(i,j,k+1)).and.(layer_type(i,j+1,k+1)==&
                  layer_type(i+1,j+1,k+1)).and.(layer_type(i,j,k+1)==layer_type(i,j+1,k+1))) then

                    type1=layer_type(i+1,j,k+1)

              else
                    type1='CE'        ! Combination element

              endif
              if((type1=='Cr').or.(type1=='CE')) then

                  element(e)%x(2)=x(i,j,k)
                  element(e)%y(2)=y(i,j,k)
                  element(e)%z(2)=z(i,j,k)
                  element(e)%x(6)=x(i,j,k+1)
                  element(e)%y(6)=y(i,j,k+1)
                  element(e)%z(6)=z(i,j,k+1)
                  element(e)%x(3)=x(i,j+1,k)
                  element(e)%y(3)=y(i,j+1,k)
                  element(e)%z(3)=z(i,j+1,k)
                  element(e)%x(7)=x(i,j+1,k+1)
                  element(e)%y(7)=y(i,j+1,k+1)
                  element(e)%z(7)=z(i,j+1,k+1)
                  element(e)%x(1)=x(i+1,j,k)
                  element(e)%y(1)=y(i+1,j,k)
                  element(e)%z(1)=z(i+1,j,k)
                  element(e)%x(5)=x(i+1,j,k+1)
                  element(e)%y(5)=y(i+1,j,k+1)
                  element(e)%z(5)=z(i+1,j,k+1)
                  element(e)%x(4)=x(i+1,j+1,k)
                  element(e)%y(4)=y(i+1,j+1,k)
                  element(e)%z(4)=z(i+1,j+1,k)
                  element(e)%x(8)=x(i+1,j+1,k+1)
                  element(e)%y(8)=y(i+1,j+1,k+1)
                  element(e)%z(8)=z(i+1,j+1,k+1)

                  element(e)%rho => cr(a)%rhoc

                  cr(a)%rhoc = refrho

                  cr(a)%elmnts(1) = e

                  ! Density of element = average of density at nodes 5, 6, 7, and 8
!                  cr(a)%rhoc=(density(i+1,j,k+1)+density(i,j,k+1)+&
!                                density(i,j+1,k+1)+density(i+1,j+1,k+1))/4.0d0


                  m=0
                  element(e)%crust=.true.
                  if(MohoP(i+1,j,k+1)==.true.) m=m+1
                  if(MohoP(i,j,k+1)==.true.) m=m+1
                  if(MohoP(i,j+1,k+1)==.true.) m=m+1
                  if(MohoP(i+1,j+1,k+1)==.true.) m=m+1
                  if(m>2) element(e)%crust=.false.          ! element below reference moho (30 km)

                  element(e)%id=a
                  element(e)%type=type1

                  if(element(e)%id==element(e-1)%id) then
                    cr(a)%double = .true.
                    cr(a)%elmnts(1) = e-1
                    cr(a)%elmnts(2) = e
                  endif
                  e=e+1

              else

                  elemento(e1)%x(2)=x(i,j,k)
                  elemento(e1)%y(2)=y(i,j,k)
                  elemento(e1)%z(2)=z(i,j,k)
                  elemento(e1)%x(6)=x(i,j,k+1)
                  elemento(e1)%y(6)=y(i,j,k+1)
                  elemento(e1)%z(6)=z(i,j,k+1)
                  elemento(e1)%x(3)=x(i,j+1,k)
                  elemento(e1)%y(3)=y(i,j+1,k)
                  elemento(e1)%z(3)=z(i,j+1,k)
                  elemento(e1)%x(7)=x(i,j+1,k+1)
                  elemento(e1)%y(7)=y(i,j+1,k+1)
                  elemento(e1)%z(7)=z(i,j+1,k+1)
                  elemento(e1)%x(1)=x(i+1,j,k)
                  elemento(e1)%y(1)=y(i+1,j,k)
                  elemento(e1)%z(1)=z(i+1,j,k)
                  elemento(e1)%x(5)=x(i+1,j,k+1)
                  elemento(e1)%y(5)=y(i+1,j,k+1)
                  elemento(e1)%z(5)=z(i+1,j,k+1)
                  elemento(e1)%x(4)=x(i+1,j+1,k)
                  elemento(e1)%y(4)=y(i+1,j+1,k)
                  elemento(e1)%z(4)=z(i+1,j+1,k)
                  elemento(e1)%x(8)=x(i+1,j+1,k+1)
                  elemento(e1)%y(8)=y(i+1,j+1,k+1)
                  elemento(e1)%z(8)=z(i+1,j+1,k+1)

                  t1=sqrt((elemento(e1)%x(5)-elemento(e1)%x(1))**2+(elemento(e1)%y(5)-&
                          elemento(e1)%y(1))**2+(elemento(e1)%z(5)-elemento(e1)%z(1))**2)
                  t2=sqrt((elemento(e1)%x(6)-elemento(e1)%x(2))**2+(elemento(e1)%y(6)-&
                          elemento(e1)%y(2))**2+(elemento(e1)%z(6)-elemento(e1)%z(2))**2)
                  t3=sqrt((elemento(e1)%x(7)-elemento(e1)%x(3))**2+(elemento(e1)%y(7)-&
                          elemento(e1)%y(3))**2+(elemento(e1)%z(7)-elemento(e1)%z(3))**2)
                  t4=sqrt((elemento(e1)%x(8)-elemento(e1)%x(4))**2+(elemento(e1)%y(8)-&
                          elemento(e1)%y(4))**2+(elemento(e1)%z(8)-elemento(e1)%z(4))**2)
                  t5=min(t1,t2,t3,t4)

                  elemento(e1)%thickness=(t1+t2+t3+t4)/4.0d0

                  l1=sqrt((elemento(e1)%x(2)-elemento(e1)%x(1))**2+(elemento(e1)%y(2)-&
                          elemento(e1)%y(1))**2+(elemento(e1)%z(2)-elemento(e1)%z(1))**2)
                  l2=sqrt((elemento(e1)%x(3)-elemento(e1)%x(2))**2+(elemento(e1)%y(3)-&
                          elemento(e1)%y(2))**2+(elemento(e1)%z(3)-elemento(e1)%z(2))**2)
                  l3=(l1+l2)/2.0d0

                  if((l3/elemento(e1)%thickness)>10.0) then
                    elemento(e1)%sheet=.true.

                  else
                      elemento(e1)%sheet=.false.
                      elemento(e1)%thickness=0.0d0
                  endif

                  ! Density of elemento = average of density at nodes 5, 6, 7, and 8
!                  elemento(e1)%rho=(density(i+1,j,k+1)+density(i,j,k+1)+&
!                                density(i,j+1,k+1)+density(i+1,j+1,k+1))/4.0d0


                  m=0
                  elemento(e1)%crust=.false.
                  if(MohoP(i+1,j,k)==.true.) m=m+1
                  if(MohoP(i,j,k)==.true.) m=m+1
                  if(MohoP(i,j+1,k)==.true.) m=m+1
                  if(MohoP(i+1,j+1,k)==.true.) m=m+1
                  if(m>2) elemento(e1)%crust=.true.          ! Mantle element above reference moho (30 km)


                  elemento(e1)%id=a
                  elemento(e1)%type=type1

                  if(elemento(e1)%type=='Ma') elemento(e1)%rho=3374.58194d0
                  if(elemento(e1)%type=='S1') elemento(e1)%rho=2110.0d0
                  if(any(salt==a)) then
                    if(elemento(e1)%type=='S2') elemento(e1)%rho=2250.0d0
                    if(elemento(e1)%type=='S3') elemento(e1)%rho=2400.0d0
                  else
                    if(elemento(e1)%type=='S2') elemento(e1)%rho=2370.0d0
                    if(elemento(e1)%type=='S3') elemento(e1)%rho=2540.0d0
                  endif


                  e1=e1+1
              endif

          end do
       end do
    end do

    print *, c_elmnts, n_elements

    deallocate(x,y,z,density,layer_type,MohoP,salt)


  end subroutine initialize_elements

!========================================================================
  subroutine initialize_observations(bouguerfile)
!========================================================================

    use global_grav
    implicit none
    integer :: n
    double precision :: h,lat,long
    double precision :: a,b,nphi,distp
    character(25)    :: bouguerfile

    a=6378.1370d0
    b=6356.7523d0

    open(10,file=bouguerfile)
    read(10,*) n_observations


    allocate(observation(n_observations))

    do n=1,n_observations

       observation(n)%id=n

       read(10,*) observation(n)%lat,observation(n)%long,h,observation(n)%bouguer_anom
       h=h/1000.0d0
       long=(observation(n)%long)*pi/180.0d0
       lat=observation(n)%lat*pi/180.0d0
       nphi=a**2/sqrt(((a**2)*(cos(lat))**2)+((b**2)*&
                (sin(lat))**2))
       observation(n)%x=(nphi+h)*cos(lat)*cos(long)
       observation(n)%y=(nphi+h)*cos(lat)*sin(long)
       observation(n)%z=((b**2/a**2)*nphi+h)*sin(lat)

       distp=sqrt((observation(n)%x)**2+&
                    (observation(n)%y)**2+(observation(n)%z)**2)

       observation(n)%ux=observation(n)%x/distp
       observation(n)%uy=observation(n)%y/distp
       observation(n)%uz=observation(n)%z/distp


    enddo
    close(10)
  end subroutine initialize_observations

!========================================================================
  subroutine initialize_gravspec3d
!========================================================================

    use omp_lib
    use global_grav
    implicit none

    integer                :: i,j
    double precision       :: sed_eff,dgr,error

    ! First remove sediment effect

    allocate(specG(n_observations,c_elmnts))

    !$omp parallel &
    !$omp default(shared) &
    !$omp private(dgr,j,sed_eff,i)
    !$omp do schedule(static,50)
    do i=1,n_observations
      sed_eff=0.0d0
      do j=1,n_elements
        if((elemento(j)%type=='S1').or.(elemento(j)%type=='S2')&
            .or.(elemento(j)%type=='S3')) then
          if(elemento(j)%sheet==.true.) then
            call SpecGrav2d(j,dgr,i)
            sed_eff=sed_eff+(dgr*(elemento(j)%rho-refrho))
          else
            call SpecGrav3d(j,dgr,i,2)
            sed_eff=sed_eff+(dgr*(elemento(j)%rho-refrho))
          endif
        else
!        if(elemento(j)%type=='Ma') then
          if(elemento(j)%crust==.true.) then
            call SpecGrav3d(j,dgr,i,2)
            sed_eff=sed_eff+(dgr*(elemento(j)%rho-refrho))
          endif
        endif
      enddo
!      print *, observation(i)%bouguer_anom,sed_eff
      observation(i)%bouguer_anom=observation(i)%bouguer_anom-sed_eff
    enddo
    !$omp end do
    !$omp end parallel


    !$omp parallel &
    !$omp default(shared) &
    !$omp private(dgr,j,i)
    !$omp do schedule(static,50)
    do i=1,n_observations
      do j=1,c_elmnts
        call SpecGrav3d(j,dgr,i,1)
        specG(i,j)=dgr
      enddo
    enddo
    !$omp end do
    !$omp end parallel


    deallocate(elemento)

    return

  end subroutine initialize_gravspec3d

!========================================================================
  subroutine smoothing_factor(el,mu,eta)
!========================================================================
!---------------------------------------------------------------------
! Calculates the smoothness applied to the objective function
!---------------------------------------------------------------------
    use global_grav
    implicit none

    integer                :: el,a,b
    double precision       :: eta,d1,d2,d3,d4,d5,d6,d7,d8,mu

    a=(nlat-2)*(nlong-1)+1
    eta=0.0d0

      if (el==1) then
        d1=cr(el)%rhoc-cr(2)%rhoc
        d2=cr(el)%rhoc-cr(nlong)%rhoc
        d3=cr(el)%rhoc-cr(nlong+1)%rhoc
        eta=mu*(d1+d2+d3)/3.0d0
      elseif (el==(nlong-1)) then
        d1=cr(el)%rhoc-cr(nlong-2)%rhoc
        d2=cr(el)%rhoc-cr(2*(nlong-1))%rhoc
        d3=cr(el)%rhoc-cr(2*nlong-3)%rhoc
        eta=mu*(d1+d2+d3)/3.0d0
      elseif (el==a) then
        d1=cr(el)%rhoc-cr(el-(nlong-1))%rhoc
        d2=cr(el)%rhoc-cr(el+1)%rhoc
        d3=cr(el)%rhoc-cr(el-nlong+2)%rhoc
        eta=mu*(d1+d2+d3)/3.0d0
      elseif (el==ncol) then
        d1=cr(el)%rhoc-cr(el-(nlong-1))%rhoc
        d2=cr(el)%rhoc-cr(el-1)%rhoc
        d3=cr(el)%rhoc-cr(el-nlong)%rhoc
        eta=mu*(d1+d2+d3)/3.0d0
      elseif (el>1 .and. el<nlong-1) then
        d1=cr(el)%rhoc-cr(el+1)%rhoc
        d2=cr(el)%rhoc-cr(el-1)%rhoc
        d3=cr(el)%rhoc-cr(el+nlong-1)%rhoc
        d4=cr(el)%rhoc-cr(el+nlong-2)%rhoc
        d5=cr(el)%rhoc-cr(el+nlong)%rhoc
        eta=mu*(d1+d2+d3+d4+d5)/5.0d0
      elseif (modulo(el,(nlong-1))==1) then
        d1=cr(el)%rhoc-cr(el+1)%rhoc
        d2=cr(el)%rhoc-cr(el+nlong-1)%rhoc
        d3=cr(el)%rhoc-cr(el-nlong+1)%rhoc
        d4=cr(el)%rhoc-cr(el+nlong)%rhoc
        d5=cr(el)%rhoc-cr(el-nlong+2)%rhoc
        eta=mu*(d1+d2+d3+d4+d5)/5.0d0
      elseif (modulo(el,(nlong-1))==0) then
        d1=cr(el)%rhoc-cr(el-1)%rhoc
        d2=cr(el)%rhoc-cr(el+nlong-1)%rhoc
        d3=cr(el)%rhoc-cr(el-nlong+1)%rhoc
        d4=cr(el)%rhoc-cr(el-nlong)%rhoc
        d5=cr(el)%rhoc-cr(el+nlong-2)%rhoc
        eta=mu*(d1+d2+d3+d4+d5)/5.0d0
      elseif (el>a .and. el<ncol) then
        d1=cr(el)%rhoc-cr(el-1)%rhoc
        d2=cr(el)%rhoc-cr(el+1)%rhoc
        d3=cr(el)%rhoc-cr(el-nlong+1)%rhoc
        d4=cr(el)%rhoc-cr(el-nlong)%rhoc
        d5=cr(el)%rhoc-cr(el-nlong+2)%rhoc
        eta=mu*(d1+d2+d3+d4+d5)/5.0d0
      else
        d1=cr(el)%rhoc-cr(el-1)%rhoc
        d2=cr(el)%rhoc-cr(el+1)%rhoc
        d3=cr(el)%rhoc-cr(el-nlong+1)%rhoc
        d4=cr(el)%rhoc-cr(el-nlong)%rhoc
        d5=cr(el)%rhoc-cr(el-nlong+2)%rhoc
        d6=cr(el)%rhoc-cr(el+nlong-1)%rhoc
        d7=cr(el)%rhoc-cr(el+nlong)%rhoc
        d8=cr(el)%rhoc-cr(el+nlong-2)%rhoc
        eta=mu*(d1+d2+d3+d4+d5+d6+d7+d8)/8.0d0
      endif


      eta=dabs(eta)

  end subroutine smoothing_factor
!========================================================================
  subroutine data_misfit
!========================================================================
! Calculates the smoothness applied to the objective function
!---------------------------------------------------------------------
    use omp_lib
    use global_grav
    implicit none

    integer                :: i,j,a
    double precision       :: dgr
    double precision       :: error

    error=0.0d0

    !$omp parallel &
    !$omp default(shared) &
    !$omp private(dgr,j,i)

    !$omp do schedule(static,50)

    do i=1,n_observations
      dgr=0.0d0
      do j=1,c_elmnts
        if(element(j)%crust==.false.) then
          dgr=dgr+(specG(i,j)*(element(j)%rho-3374.58194d0))
        else
          dgr=dgr+(specG(i,j)*(element(j)%rho-refrho))
        endif
      enddo
      observation(i)%dgr=dgr

    enddo

    !$omp end do
    !$omp end parallel


    !$omp parallel do schedule(static,50) reduction(+:error)
    do i=1,n_observations
        error=error+(observation(i)%bouguer_anom-observation(i)%dgr)**2
    enddo
    !$omp end parallel do

    rmserr=sqrt(error/n_observations)     ! data misfit


    print *, 'The RMS is:'
    print *, rmserr
     

  end subroutine data_misfit

!========================================================================
  subroutine SpecGrav2d(iels,dgr,iobs)
!========================================================================
! Computes the surface element stiffness matrix for the 3D DC electrical
! resistivty problem for a unit local conductivity. To complete the
! stiffness calculation, each term in the the stiffness matrix must be
! scaled by the conductivity of the element.
!---------------------------------------------------------------------
  use global_grav
  implicit none
  integer :: iels,nn,ni,i,j,k,ip,iobs
  double precision :: x3d(3,4),node(2,4),xi(4,2),wt(4)
  double precision :: dpsix(2,4),dsdx(2,2),detj
  double precision :: psi(4),dpsi(4,2)
  double precision :: d12,d13,d23,d14,d24
  double precision :: dx,dy,dz
  double precision :: dgr,dg,xip,yip,zip,dist3


  !Initialize element stiffness matrix
  dgr=0.0d0 ! gravitational contribution due to the current element

  !Get the gauss integration points for the current element
  nn=4

  call gaus2d(nn,ni,xi,wt)


  !Recover the node coordinates and conductivities
  do i=1,nn
    x3d(1,i)=(elemento(iels)%x(i)+elemento(iels)%x(i+4))/2.0d0
    x3d(2,i)=(elemento(iels)%y(i)+elemento(iels)%y(i+4))/2.0d0
    x3d(3,i)=(elemento(iels)%z(i)+elemento(iels)%z(i+4))/2.0d0
  end do


  zip=(x3d(3,1)+x3d(3,2)+x3d(3,3)+x3d(3,4))/4.0d0
  !Began integration of the sub-stiffness matrices
  do ip=1,ni
    !compute the value of the trilinear shape function at this integration point
    call shape2d(nn,xi(ip,1),xi(ip,2),psi,dpsi)

    !Calculate dxds and detj
    call trnfrm2d(nn,x3d,dpsi,dsdx,detj)
    if(detj.le.0.0) then
      print *, x3d(1,1),x3d(2,1)
      print *, x3d(1,2),x3d(2,2)
      print *, x3d(1,3),x3d(2,3)
      print *, x3d(1,4),x3d(2,4)
      write(*,*) 'Bad Jacobian in surface element ',iels
      write(*,*) 'Program terminated.'
      write(*,*) 'Press enter to exit the program...'
      read(*,*)
!      close(outunt)
      stop
    end if

    xip=0.0d0
    yip=0.0d0
    do i=1,4
      xip=xip+x3d(1,i)*psi(i)
      yip=yip+x3d(2,i)*psi(i)
    end do
    !Accumulate integration point value of integrals
    dist3=(sqrt((observation(iobs)%x-xip)**2+(observation(iobs)%y-yip)**2&
              +(observation(iobs)%z-zip)**2))**3

    dg=((observation(iobs)%x-xip)*observation(iobs)%ux)+((observation(iobs)%y-yip)&
        *observation(iobs)%uy)+((observation(iobs)%z-zip)*observation(iobs)%uz)

    dgr=dgr+(ugc*detj*wt(ip)*dg*elemento(iels)%thickness/dist3)

  end do

  return
  end subroutine SpecGrav2d

!========================================================================
  subroutine shape2d(n,eps,eta,psi,dpsi)
!========================================================================
! Compute finite element shape functions for 2D elements
  implicit none
  integer :: n
  double precision :: eps,eta,psi(4),dpsi(4,2)
  double precision :: x,x2,e,e2,z1,z2,z3

  if(n==3) then
    !bilinear shape functions for a 3-node triangle element
    psi(1)=1.0d0-eps-eta
    psi(2)=eps
    psi(3)=eta
    dpsi(1,1)=0.0d0
    dpsi(1,2)=0.0d0
    dpsi(2,1)=1.0d0
    dpsi(2,2)=0.0d0
    dpsi(3,1)=0.0d0
    dpsi(3,2)=1.0d0
  else if(n==4) then
    !bilinear shape functions for the 4-node quadrilateral element
    x=eps
    e=eta
    psi(1)=0.25d0*(1.0d0-x)*(1.0d0-e)
    psi(2)=0.25d0*(1.0d0+x)*(1.0d0-e)
    psi(3)=0.25d0*(1.0d0+x)*(1.0d0+e)
    psi(4)=0.25d0*(1.0d0-x)*(1.0d0+e)
    dpsi(1,1)=0.25d0*(-1.0d0+e)
    dpsi(1,2)=0.25d0*(-1.0d0+x)
    dpsi(2,1)=0.25d0*( 1.0d0-e)
    dpsi(2,2)=0.25d0*(-1.0d0-x)
    dpsi(3,1)=0.25d0*( 1.0d0+e)
    dpsi(3,2)=0.25d0*( 1.0d0+x)
    dpsi(4,1)=0.25d0*(-1.0d0-e)
    dpsi(4,2)=0.25d0*( 1.0d0-x)
  else
    !report no element shapes found
    write(7,*) 'Error in call to shape,n=',n
    write(7,*) 'Run termianted.'
    stop
  end if
  return
  end subroutine shape2d

!========================================================================
  subroutine trnfrm2d(n,x,dpsi,dsdx,detj)
!========================================================================
! compute the transformation matrix and its determinent
! at the local point xi(1),xi(2) in a 2D element
!------------------------------------------------------------------------
  implicit none
  integer :: i,j,k,n
  double precision :: x(3,4),dpsi(4,2)
  double precision :: dsdx(2,2),detj,dxds(2,2)
  double precision :: detji

  !calculate the jacobian matrix of derivatives of the global
  !nodal coordinates x, relative to local coordinates xi
  do i=1,2
    do j=1,2
      dxds(i,j)=0.0d0
      do k=1,n
        dxds(i,j)=dxds(i,j)+dpsi(k,j)*x(i,k)
      end do
    end do
  end do

  !calculate the determinant of the jacobian matrix
  detj=dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)

  !test for a bad jacobian
  if(detj.le.0.0d0) return

  !compute the transformation matrix dsdx
  detji=1.0d0/detj
  dsdx(1,1)=dxds(2,2)*detji
  dsdx(2,2)=dxds(1,1)*detji
  dsdx(1,2)=-dxds(1,2)*detji
  dsdx(2,1)=-dxds(2,1)*detji

  return
  end subroutine trnfrm2d

!========================================================================
    subroutine gaus2d(nn,ni,xi,wt)
!========================================================================
!   routine to set quadrature points for 2-d elements
   implicit none
   integer :: i,j,k,ni,nn
   double precision :: xi(4,2),wt(4)
   double precision :: thrd,sq,xi1,xi2

   !one-point quadrature for 3-node triangle elements
   if(nn==3) then
     ni=1
     xi(1,1)=1.0d0/3.0d0
     xi(1,2)=1.0d0/3.0d0
     wt(1)=0.5d0
   !four-point quadrature for 4-node quad elements
   else if(nn==4) then
     ni=4
     sq=dsqrt(1.0d0/3.0d0)
     xi(1,1)=sq
     xi(1,1)=sq
     xi(2,1)=-sq
     xi(2,2)=sq
     xi(3,1)=-sq
     xi(3,2)=-sq
     xi(4,1)=sq
     xi(4,2)=-sq
     wt=1.0d0
   else
      write(7,*) 'Error: Unrecognized 2D integration point rule.'
      write(7,*) 'Run terminated.'
      write(*,*) 'Error:  Unrecognized 2D integration point rule.'
      write(*,*) 'Press enter to stop the program...'
      read(*,*)
      stop
    end if
    return
    end subroutine gaus2d

!========================================================================
   subroutine writeelements(crustfile)
!========================================================================
!   routine to set quadrature points for 2-d elements

    use global_grav
    implicit none
    integer :: i,j,k,l,m,e1
    double precision :: den,lat,long
    double precision :: a,b,nphi,distp
    double precision :: x1,y1,z1,xx,yy,zz
    character(25)    :: crustfile

    open(unit=11,file=crustfile)

    do i=1,ncol

      e1=cr(i)%elmnts(1)

      x1=0.0d0
      y1=0.0d0
      z1=0.0d0
      do j=1,4
        x1=x1+((element(e1)%x(j)+element(e1)%x(j+4))/2.0d0)
        y1=y1+((element(e1)%y(j)+element(e1)%y(j+4))/2.0d0)
        z1=z1+((element(e1)%z(j)+element(e1)%z(j+4))/2.0d0)
      enddo
      xx=x1/4.0d0
      yy=y1/4.0d0
      zz=z1/4.0d0

      call olson(lat,long,xx,yy,zz)
      write(11,*) lat,long,cr(i)%rhoc

    enddo

    close(11)
    close(12)
    close(13)
    close(14)

    return
   end subroutine writeelements


!========================================================================
   subroutine writeobservations(obsfile)
!========================================================================
!   routine to set quadrature points for 2-d elements
    use global_grav
    implicit none
    integer          :: i
    double precision :: lat,long
    character(25)    :: obsfile


!    !$omp parallel default(shared) private(lat,long)
!    !$omp do schedule(static,1200)
!    do i=1,n_observations
!        call olson(lat,long,observation(i)%x,&
!                observation(i)%y,observation(i)%z)
!        observation(i)%lat=lat
!        observation(i)%long=long
!    enddo
!    !$omp end do
!    !$omp end parallel

    open(unit=18,file=obsfile,status='new')
    do i=1,n_observations
      write(18,*) observation(i)%lat,observation(i)%long,observation(i)%bouguer_anom,observation(i)%dgr
    enddo

    close(18)

    return
   end subroutine writeobservations

!*****************************************************************************************
!  Heikkinen routine for cartesian to geodetic transformation
!
!# References
!  1. M. Heikkinen, "Geschlossene formeln zur berechnung raumlicher
!     geodatischer koordinaten aus rechtwinkligen Koordinaten".
!     Z. Ermess., 107 (1982), 207-211 (in German).
!  2. E. D. Kaplan, "Understanding GPS: Principles and Applications",
!     Artech House, 1996.

!    pure subroutine heikkinen(rvec, a, b, h, lon, lat)
!
!    implicit none
!
!    double precision,dimension(3),intent(in) :: rvec  !! position vector [km]
!    double precision,intent(in)  :: a                 !! geoid semimajor axis [km]
!    double precision,intent(in)  :: b                 !! geoid semiminor axis [km]
!    double precision,intent(out) :: h                 !! geodetic altitude [km]
!    double precision,intent(out) :: lon               !! longitude [rad]
!    double precision,intent(out) :: lat               !! geodetic latitude [rad]
!
!    double precision :: f,e_2,ep,r,e2,ff,g,c,s,pp,q,r0,u,v,z0,x,y,z,z2,r2,tmp,a2,b2
!
!    x   = rvec(1)
!    y   = rvec(2)
!    z   = rvec(3)
!    a2  = a*a
!    b2  = b*b
!    f   = (a-b)/a
!    e_2 = (2.0d0*f-f*f)
!    ep  = sqrt(a2/b2 - 1.0d0)
!    z2  = z*z
!    r   = sqrt(x**2 + y**2)
!    r2  = r*r
!    e2  = a2 - b2
!    ff  = 54.0d0 * b2 * z2
!    g   = r2 + (1.0d0 - e_2)*z2 - e_2*e2
!    c   = e_2**2 * ff * r2 / g**3
!    s   = (1.0d0 + c + sqrt(c**2 + 2.0d0*c))**(1.0d0/3.0d0)
!    pp  = ff / ( 3.0d0*(s + 1.0d0/s + 1.0d0)**2 * g**2 )
!    q   = sqrt( 1.0d0 + 2.0d0*e_2**2 * pp )
!    r0  = -pp*e_2*r/(1.0d0+q) + &
!            sqrt( max(0.0d0, 1.0d0/2.0d0 * a2 * (1.0d0 + 1.0d0/q) - &
!                ( pp*(1.0d0-e_2)*z2 )/(q*(1.0d0+q)) - &
!                1.0d0/2.0d0 * pp * r2) )
!    u   = sqrt( (r - e_2*r0)**2 + z2 )
!    v   = sqrt( (r - e_2*r0)**2 + (1.0d0 - e_2)*z2 )
!    z0  = b**2 * z / (a*v)
!
!    h   = u*(1.0d0 - b2/(a*v) )
!    lat = atan2( (z + ep**2*z0), r )
!    lon = atan2( y, x )
!
!    end subroutine heikkinen
!*****************************************************************************************

!*****************************************************************************************
!  Olson routine for cartesian to geodetic transformation.
!
!# References
!  1. Olson, D. K., Converting Earth-Centered, Earth-Fixed Coordinates to
!     Geodetic Coordinates, IEEE Transactions on Aerospace and Electronic
!     Systems, 32 (1996) 473-476.

    subroutine olson(lat,long,x,y,z)

    implicit none

    double precision  :: a                !!geoid semimajor axis [km]
    double precision  :: b                !!geoid semiminor axis [km]
    double precision :: h                !!geodetic altitude [km]
    double precision,intent(out) :: long             !!longitude [rad]
    double precision,intent(out) :: lat              !!geodetic latitude [rad]

    double precision :: f,x,y,z,e2,a1,a2,a3,a4,a5,a6,w,zp,&
                w2,r2,r,s2,c2,u,v,s,ss,c,g,rg,rf,m,p,z2
    double precision,parameter   :: pi=3.1415926535


    a=6378.1370d0
    b=6356.7523d0

    f  = (a-b)/a
    e2 = f * (2.0d0 - f)
    a1 = a * e2
    a2 = a1 * a1
    a3 = a1 * e2 / 2.0d0
    a4 = 2.5d0 * a2
    a5 = a1 + a3
    a6 = 1.0d0 - e2
    zp = abs(z)
    w2 = x*x + y*y
    w  = sqrt(w2)
    z2 = z * z
    r2 = z2 + w2
    r  = sqrt(r2)

    if (r < 100.0d0) then

        lat = 0.0d0
        long = 0.0d0
!        h = -1.0e7d0

    else

        s2 = z2 / r2
        c2 = w2 / r2
        u  = a2 / r
        v  = a3 - a4 / r

        if (c2 > 0.3d0) then
            s = (zp / r) * (1.0d0 + c2 * (a1 + u + s2 * v) / r)
            lat = asin(s)
            ss = s * s
            c = sqrt(1.0d0 - ss)
        else
            c = (w / r) * (1.0d0 - s2 * (a5 - u - c2 * v) / r)
            lat = acos(c)
            ss = 1.0d0 - c * c
            s = sqrt(ss)
        end if

        g   = 1.0d0 - e2 * ss
        rg  = a / sqrt(g)
        rf  = a6 * rg
        u   = w - rg * c
        v   = zp - rf * s
        f   = c * u + s * v
        m   = c * v - s * u
        p   = m / (rf / g + f)
        lat = lat + p
        if (z < 0.0d0) lat = -lat
        h = f + m * p / 2.0d0
        long = atan2( y, x )

    end if

    long=(long*180.0d0/pi)
    lat=(lat*180.0d0/pi)

    return

end subroutine olson

!========================================================================

