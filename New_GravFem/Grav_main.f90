! CODE FOR  SIMULATED ANNEALING OF GRAVITY DATA USING SIMULATED ANNEALING
! VERSION:
!      1.0
! AUTHOR:
!      Alden Netto
! PURPOSE:
!
!
! HISTORY:
!       2018/03/10
!


program Grav_main

  use global_grav
  implicit none

  integer                       :: i,j,n,iseed,b,m,ii,e1,e2,maxobs
  integer                       :: maxruns,maxiter,tempiter
  character(25)                 :: obsfile,crustfile,bouguerfile
  double precision              :: u,mu,eta,d,alpha
  double precision              :: dgr,ch,objc1,objc2
  double precision              :: rho1,rho2,q,T,T0,P,rhot
  double precision              :: obj0,rmstol,error,err1,maxerr
!  logical                      :: ichange

  open(10,file='GravFEM.in')
  do i=1,3
    read(10,*)
  enddo
  read(10,*) mu         ! Smoothing Factor
  read(10,*) rmstol
  read(10,*) maxiter
  read(10,*) tempiter
  read(10,*) alpha
  read(10,*) maxruns
  read(10,*) iseed
  read(10,*) refrho
  read(10,*) refmoho
  read(10,*) crustfile
  read(10,*) obsfile
  read(10,*) bouguerfile

  close(10)

  call initialize_parameter_grids

  call initialize_cartesian_grids

  call initialize_elements

  call initialize_observations(bouguerfile)

  call initialize_gravspec3d



!record the progress in reducing the objective function versus
  !iteration number
  write(*,*) '              Simulated Annealing Results'
!  write(*,*) 'Iteration No.  Temperature    RMS Error       Objective      Best Objective'

            ! xxxxxx xxxxxxxxxxx  xxxxxx  xxxxxxxxxxx



  !Compute the initial value of the objective function
  call smoothing_factor(1,mu,eta)
  call data_misfit

!  err1=0.0d0
!  maxerr=0.0d0
!  do ii=1,n_observations
!      err1=abs(observation(ii)%bouguer_anom-observation(ii)%dgr)
!      if(maxerr<err1) maxerr=err1
!  enddo
  obj0=rmserr+eta
  T0=obj0*3.0d0
  T=T0
  !report the initial status of the model
  write(*,*) 'Initial Temperature    RMS Error       Objective       Max Error'
  write(*,*) T0,rmserr,obj0,maxerr

  !Test for the unlikely meeting the RMS error tolerance

  !loop over temperature steps

  m=0

  do while (T>0)

    if((dabs(obj0)<=rmstol).or.(m>maxruns)) then
    !report convergence
      write(*,*) 'Convergence reached after',m,' iterations'
      call writeelements(crustfile)
      call writeobservations(obsfile)
      print *, 'Finished writing elements'
      print *, 'Smoothing Factor: ', mu
      print *, 'The Reference crustal density: ', refrho
      print *, 'Reference Moho depth: ', refmoho
      print *, 'Input gravity measurements: ', bouguerfile
      stop
    end if

    do i=1,ncol
      rho1=cr(i)%rhoc
      rho2=cr(i)%rhoc
iterloop: do n=1,maxiter
        u=ran(iseed)
        !Compute a random model perturbation size
        ch=T*sign(1.0d0,(u-0.5d0))*((1.0d0+1.0d0/T)**dabs(2.0d0*u-1.0d0)-1.0d0)

        rhot=cr(i)%rhoc
        cr(i)%rhoc=cr(i)%rhoc+(ch*380.0d0)
        if(cr(i)%rhoc<2670) cr(i)%rhoc=2670.0d0+(u*50.0d0)
        if(cr(i)%rhoc>3050) cr(i)%rhoc=3050.0d0-(u*50.0d0)


        if(cr(i)%double == .true.) then
          e1=cr(i)%elmnts(1)
          e2=cr(i)%elmnts(2)
          !$omp parallel default(shared)
          !$omp do schedule(static,400)
          do ii=1,n_observations
            observation(ii)%dgr=observation(ii)%dgr+((specG(ii,e1)+specG(ii,e2))*(cr(i)%rhoc-rhot))
          enddo
          !$omp end do
          !$omp end parallel
        else
          e1=cr(i)%elmnts(1)
          !$omp parallel default(shared)
          !$omp do schedule(static,400)
          do ii=1,n_observations
            observation(ii)%dgr=observation(ii)%dgr+(specG(ii,e1)*(cr(i)%rhoc-rhot))
          enddo
          !$omp end do
          !$omp end parallel
        endif

        error=0.0d0
        !$omp parallel do schedule(static,50) reduction(+:error)
        do ii=1,n_observations
            error=error+(observation(ii)%bouguer_anom-observation(ii)%dgr)**2
        enddo
        !$omp end parallel do

        rmserr=sqrt(error/n_observations)

        call smoothing_factor(i,mu,eta)

        objc2=rmserr+eta

        if (n==1) objc1=objc2

        if(objc2<obj0) then
          obj0=objc2
          exit iterloop
        elseif (objc2<objc1) then
          objc1=objc2
          rho2=cr(i)%rhoc
        endif


        if (n==maxiter) then
          u=ran(iseed)
          P=exp(-1*((objc1-obj0)/T))  !(1.0d0-(1.0d0-q)*(objc1-obj0)/T)**(1.0d0-q)
          rhot=cr(i)%rhoc
          if (u<P) then
            cr(i)%rhoc=rho2
            obj0=objc1
          elseif (u>=P) then
            cr(i)%rhoc=rho1
          endif
          if(cr(i)%double == .true.) then
            e1=cr(i)%elmnts(1)
            e2=cr(i)%elmnts(2)
            !$omp parallel default(shared)
            !$omp do schedule(static,50)
            do ii=1,n_observations
              observation(ii)%dgr=observation(ii)%dgr+((specG(ii,e1)+specG(ii,e2))*(cr(i)%rhoc-rhot))
            enddo
            !$omp end do
            !$omp end parallel
          else
            e1=cr(i)%elmnts(1)
            !$omp parallel default(shared)
            !$omp do schedule(static,50)
            do ii=1,n_observations
              observation(ii)%dgr=observation(ii)%dgr+(specG(ii,e1)*(cr(i)%rhoc-rhot))
            enddo
            !$omp end do
            !$omp end parallel
          endif
        endif

      enddo iterloop ! maxiter
    enddo

    err1=0.0d0
    maxerr=0.0d0
    do ii=1,n_observations
      err1=dabs(observation(ii)%bouguer_anom-observation(ii)%dgr)
      if(maxerr<err1) then
        maxerr=err1
        maxobs=ii
      endif
    enddo
    write(*,*) T,rmserr,obj0,maxerr,m
    m=m+1
    T=T0*exp(-1*alpha*m)
    write(*,*)

  enddo     ! do while



end program Grav_main
