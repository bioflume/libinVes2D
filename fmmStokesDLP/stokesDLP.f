c     forming the pressures should be done two at a time in order
c     to half the total number of FMM calls
c     In fact, all of the calls should be able to be done two
c     at a time

      subroutine stokesDLP(s1,s2,xs,ys,ns,
     $   dir1,dir2,u1,u2)
c     Evaluates the stresslet layer potential for Stokes' equation
c     (s1,s2) is the source strength of length ns
c     xs and ys are the source and target locations of length ns
c     (u1,u2) are the two components of the velocity field
      implicit real*8 (a-h,o-z)

      real *8 s1(ns),s2(ns)
c     strength of Stokes SLP
      real *8 xs(ns),ys(ns)
c     location of the source/target points
      real *8 dir1(ns),dir2(ns)
c     x and y components of the normal vector
      real *8 u1(ns),u2(ns)
c     x and y components of the velocity field
      complex *16 eye
      
      real *8, allocatable :: velTemp1(:),velTemp2(:)
      real *8, allocatable :: velTemp3(:),velTemp4(:)
      real *8, allocatable :: source(:,:) 
      real *8, allocatable :: density1(:)
      real *8, allocatable :: density2(:)
      complex *16, allocatable :: cdensity(:)

      allocate(velTemp1(ns),stat=ierr)
      allocate(velTemp2(ns),stat=ierr)
      allocate(velTemp3(ns),stat=ierr)
      allocate(velTemp4(ns),stat=ierr)
      allocate(source(2,ns),stat=ierr)
      allocate(density1(ns),stat=ierr)
      allocate(density2(ns),stat=ierr)
      allocate(cdensity(ns),stat=ierr)
c     allocate memory for temporary variables

      eye = (0.d0,1.d0)
      iprec = 4

      zeroCheck = 0.d0
      do i = 1,ns
        zeroCheck = zeroCheck + abs(s1(i)) + abs(s2(i))
        u1(i) = 0.d0
        u2(i) = 0.d0
      enddo
c     start with velocity field being 0

      if (zeroCheck .le. 1.0d-14) return 

      do i = 1,ns
        source(1,i) = xs(i)
        source(2,i) = ys(i)
      enddo

      do i = 1,ns
        density1(i) = s1(i)*dir1(i)
        density2(i) = s1(i)*dir2(i)
      enddo
      call derivGVer2(1,1,ns,source,density1,density2,
     $    iprec,velTemp1,velTemp2)
      call updateVel(ns,u1,velTemp1,2.d0)
      call updateVel(ns,u1,velTemp2,1.d0)

      call derivGVer2(1,2,ns,source,density1,density2,
     $    iprec,velTemp1,velTemp2)
      call updateVel(ns,u2,velTemp1,2.d0)
      call updateVel(ns,u2,velTemp2,1.d0)


      do i = 1,ns
        density1(i) = s2(i)*dir1(i)
        density2(i) = s2(i)*dir2(i)
      enddo
      call derivGVer2(2,1,ns,source,density1,density2,
     $    iprec,velTemp1,velTemp2)
      call updateVel(ns,u1,velTemp1,1.d0)
      call updateVel(ns,u1,velTemp2,2.d0)

      call derivGVer2(2,2,ns,source,density1,density2,
     $    iprec,velTemp1,velTemp2)
      call updateVel(ns,u2,velTemp1,1.d0)
      call updateVel(ns,u2,velTemp2,2.d0)


      do i = 1,ns
        density1(i) = s1(i)*dir2(i)
      enddo

      call derivG(2,1,1,ns,source,density1,iprec,velTemp1)
      call updateVel(ns,u1,velTemp1,1.d0)
      call derivG(2,2,1,ns,source,density1,iprec,velTemp1)
      call updateVel(ns,u2,velTemp1,1.d0)

      do i = 1,ns
        density1(i) = s2(i)*dir1(i)
      enddo

      call derivG(1,1,2,ns,source,density1,iprec,velTemp1)
      call updateVel(ns,u1,velTemp1,1.d0)
      call derivG(1,2,2,ns,source,density1,iprec,velTemp1)
      call updateVel(ns,u2,velTemp1,1.d0)

      do i = 1,ns
        cdensity(i) = s1(i)*dir1(i) + eye*s2(i)*dir2(i)
      enddo

      call press(ns,source,cdensity,iprec,
     $    velTemp1,velTemp2,velTemp3,velTemp4)
      call updateVel(ns,u1,velTemp1,-1.d0)
      call updateVel(ns,u1,velTemp3,-1.d0)
      call updateVel(ns,u2,velTemp2,-1.d0)
      call updateVel(ns,u2,velTemp4,-1.d0)

      do i=1,ns
        u1(i) = -7.957747154594767d-2*u1(i)
        u2(i) = -7.957747154594767d-2*u2(i)
      enddo


      deallocate(velTemp1)
      deallocate(velTemp2)
      deallocate(velTemp3)
      deallocate(velTemp4)
      deallocate(source)
      deallocate(density1)
      deallocate(density2)

 1020 format (ES9.2)

      end

c********************************************************
      subroutine updateVel(ns,vel,velTemp,const)
      implicit real *8 (a-h,o-z)

      dimension vel(ns),velTemp(ns)
      do i = 1,ns
        vel(i) = vel(i) + const*velTemp(i)
      enddo

      return
      end

c********************************************************
      subroutine press(ns,source,charge,iprec,
     $      pressure1,pressure2,pressure3,pressure4)

      implicit real *8 (a-h,o-z)

      complex *16 charge(ns)
      dimension pressure1(ns),pressure2(ns)
      dimension pressure3(ns),pressure4(ns)
      dimension source(2,ns)

      complex *16, allocatable :: grad(:,:)

      allocate(grad(2,ns),stat=ierr)

      ifpot = 0 ! don't need the potential 
      ifgrad = 1 ! need gradient
      ifhess = 0 ! don't need the hessian
      ifcharge = 1 ! have charge component
      ifdipole = 0 ! no dipole component

      call lfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)

      do n=1,ns
        pressure1(n) = 2.d0*real(grad(1,n))
        pressure2(n) = 2.d0*real(grad(2,n))
        pressure3(n) = 2.d0*imag(grad(1,n))
        pressure4(n) = 2.d0*imag(grad(2,n))
      enddo


      return
      end

c********************************************************
      subroutine derivGVer2(i,j,ns,source,density1,density2,
     $    iprec,DG1,DG2)

      implicit real *8 (a-h,o-z)

      dimension density1(ns),density2(ns)
      dimension DG1(ns),DG2(ns)
      dimension source(2,ns)

      complex *16 eye
      complex *16, allocatable :: charge(:),dipstr(:)
      real *8, allocatable :: dipvec(:,:)
      complex *16, allocatable :: pot1(:),pot2(:)
      complex *16, allocatable :: grad1(:,:),grad2(:,:)

      allocate(charge(ns),stat=ierr)
      allocate(dipstr(ns),stat=ierr)
      allocate(dipvec(2,ns),stat=ierr)
      allocate(pot1(ns),stat=ierr)
      allocate(pot2(ns),stat=ierr)
      allocate(grad1(2,ns),stat=ierr)
      allocate(grad2(2,ns),stat=ierr)

      eye = (0.d0,1.d0)

      ifhess = 0 ! don't need the hessian

c****************************
c     (i,j) = (1,1)
      if (i .eq. 1 .and. j .eq. 1) then
        ifpot = 0
        ifgrad = 1
        ifcharge = 1 ! need charge component
        ifdipole = 1 ! need dipole component
        do n=1,ns
          charge(n) = -density1(n) - eye*density2(n)
          dipstr(n) = density1(n) + eye*density2(n)
          dipvec(1,n) = source(1,n)
          dipvec(2,n) = 0.d0
        enddo

        call lfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess)

        ifcharge = 0 ! no longer need charge component
        ifpot = 1 ! need the potential

        do n=1,ns
          dipvec(1,n) = -1.d0
        enddo

        call lfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad2,ifhess,hess)
c     compute the second componet in the Stokes SLP

        do n=1,ns
          DG1(n) = real(grad1(1,n)) + source(1,n)*real(grad2(1,n)) + 
     $           real(pot2(n))
          DG2(n) = imag(grad1(2,n)) + source(1,n)*imag(grad2(2,n))
        enddo

c****************************
c     (i,j) = (1,2) or (i,j) = (2,1)
      elseif ((i .eq. 1 .and. j .eq. 2) .or. 
     $         i .eq. 2 .and. j .eq. 1)  then
        ifpot = 0
        ifgrad = 1
        ifcharge = 0 ! don't want charge component
        ifdipole = 1 ! need dipole component

        do n=1,ns
          dipstr(n) = density1(n) + eye*density2(n)
          dipvec(1,n) = 0.d0
          dipvec(2,n) = source(1,n)
        enddo

        call lfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess)

        ifpot = 1 ! need the potential

        do n=1,ns
          dipvec(2,n) = -1.d0
        enddo

        call lfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad2,ifhess,hess)
c     compute the second componet in the Stokes SLP

        do n=1,ns
          DG1(n) = real(grad1(1,n)) + source(1,n)*real(grad2(1,n)) +
     $          pot2(n)
          DG2(n) = imag(grad1(2,n)) + source(1,n)*imag(grad2(2,n))
        enddo


c****************************
c     (i,j) = (2,2)
      elseif (i .eq. 2 .and. j .eq. 2) then
        ifpot = 0
        ifgrad = 1
        ifcharge = 1 ! need charge component
        ifdipole = 1 ! need dipole component
        do n=1,ns
          charge(n) = -density1(n) - eye*density2(n)
          dipstr(n) = density1(n) + eye*density2(n)
          dipvec(1,n) = 0.d0
          dipvec(2,n) = source(2,n)
        enddo

        call lfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess)

        ifcharge = 0 ! no longer need charge component
        ifpot = 1 ! need the potential

        do n=1,ns
          dipvec(2,n) = -1.d0
        enddo

        call lfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad2,ifhess,hess)
c     compute the second componet in the Stokes SLP

        do n=1,ns
          DG1(n) = real(grad1(1,n)) + source(2,n)*real(grad2(1,n))
          DG2(n) = imag(grad1(2,n)) + source(2,n)*imag(grad2(2,n)) + 
     $           imag(pot2(n))
        enddo

      endif


      return
      end

c********************************************************
      subroutine derivG(i,j,k,ns,source,density,iprec,DG)

      implicit real *8 (a-h,o-z)

      dimension density(ns),DG(ns)
      dimension source(2,ns)

      real *8, allocatable :: charge(:),dipstr(:)
      real *8, allocatable :: dipvec(:,:)
      real *8, allocatable :: pot1(:),pot2(:)
      real *8, allocatable :: grad1(:,:),grad2(:,:)

      allocate(charge(ns),stat=ierr)
      allocate(dipstr(ns),stat=ierr)
      allocate(dipvec(2,ns),stat=ierr)
      allocate(pot1(ns),stat=ierr)
      allocate(pot2(ns),stat=ierr)
      allocate(grad1(2,ns),stat=ierr)
      allocate(grad2(2,ns),stat=ierr)

      ifhess = 0 ! don't need the hessian

c****************************
c     (i,j,k) = (1,1,1)
      if (i .eq. 1 .and. j .eq. 1 .and. k .eq. 1) then
        ifpot = 0
        ifgrad = 1
        ifcharge = 1 ! need charge component
        ifdipole = 1 ! need dipole component
        do n=1,ns
          charge(n) = -density(n)
          dipstr(n) = source(1,n)
          dipvec(1,n) = density(n)
          dipvec(2,n) = 0.d0
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess)

        ifcharge = 0 ! no longer need charge component
        ifpot = 1 ! need the potential

        do n=1,ns
          dipstr(n) = -1.d0
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad2,ifhess,hess)
c     compute the second componet in the Stokes SLP

        do n=1,ns
          DG(n) = grad1(1,n) + source(1,n)*grad2(1,n) + pot2(n)
        enddo

c****************************
c     (i,j,k) = (1,1,2)
      elseif (i .eq. 1 .and. j .eq. 1 .and. k .eq. 2) then
        ifpot = 0
        ifgrad = 1
        ifcharge = 1 ! need charge component
        ifdipole = 1 ! need dipole component
        do n=1,ns
          charge(n) = -density(n)
          dipstr(n) = source(1,n)
          dipvec(1,n) = density(n)
          dipvec(2,n) = 0.d0
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess)

        ifcharge = 0 ! no longer need charge component

        do n=1,ns
          dipstr(n) = -1.d0
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad2,ifhess,hess)
c     compute the second componet in the Stokes SLP

        do n=1,ns
          DG(n) = grad1(2,n) + source(1,n)*grad2(2,n)
        enddo

c****************************
c     (i,j,k) = (1,2,1) or (i,j,k) = (2,1,1)
      elseif ((i .eq. 1 .and. j .eq. 2 .and. k .eq. 1) .or.
     $        (i .eq. 2 .and. j .eq. 1 .and. k .eq. 1)) then
        ifpot = 0
        ifgrad = 1
        ifcharge = 0 ! don't want charge component
        ifdipole = 1 ! need dipole component
        do n=1,ns
          dipstr(n) = source(1,n)
          dipvec(1,n) = 0.d0
          dipvec(2,n) = density(n)
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess)

        ifpot = 1 ! need the potential

        do n=1,ns
          dipstr(n) = -1.d0
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad2,ifhess,hess)
c     compute the second componet in the Stokes SLP

        do n=1,ns
          DG(n) = grad1(1,n) + source(1,n)*grad2(1,n) + pot2(n)
        enddo

c****************************
c     (i,j,k) = (1,2,2) or (i,j,k) = (2,1,2)
      elseif ((i .eq. 1 .and. j .eq. 2 .and. k .eq. 2) .or.
     $        (i .eq. 2 .and. j .eq. 1 .and. k .eq. 2)) then
        ifpot = 0
        ifgrad = 1
        ifcharge = 0 ! don't want charge component
        ifdipole = 1 ! need dipole component
        do n=1,ns
          dipstr(n) = source(1,n)
          dipvec(1,n) = 0.d0
          dipvec(2,n) = density(n)
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess)

        do n=1,ns
          dipstr(n) = -1.d0
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad2,ifhess,hess)
c     compute the second componet in the Stokes SLP

        do n=1,ns
          DG(n) = grad1(2,n) + source(1,n)*grad2(2,n)
        enddo


c****************************
c     (i,j,k) = (2,2,1)
      elseif (i .eq. 2 .and. j .eq. 2 .and. k .eq. 1) then
        ifpot = 0
        ifgrad = 1
        ifcharge = 1 ! need charge component
        ifdipole = 1 ! need dipole component
        do n=1,ns
          charge(n) = -density(n)
          dipstr(n) = source(2,n)
          dipvec(1,n) = 0.d0
          dipvec(2,n) = density(n)
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess)

        ifcharge = 0 ! no longer need charge component

        do n=1,ns
          dipstr(n) = -1.d0
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad2,ifhess,hess)
c     compute the second componet in the Stokes SLP

        do n=1,ns
          DG(n) = grad1(1,n) + source(2,n)*grad2(1,n)
        enddo


c****************************
c     (i,j,k) = (2,2,2)
      elseif (i .eq. 2 .and. j .eq. 2 .and. k .eq. 2) then
        ifpot = 0
        ifgrad = 1
        ifcharge = 1 ! need charge component
        ifdipole = 1 ! need dipole component
        do n=1,ns
          charge(n) = -density(n)
          dipstr(n) = source(2,n)
          dipvec(1,n) = 0.d0
          dipvec(2,n) = density(n) 
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess)

        ifcharge = 0 ! no longer need charge component
        ifpot = 1 ! need the potential

        do n=1,ns
          dipstr(n) = -1.d0
        enddo

        call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot2,ifgrad,grad2,ifhess,hess)
c     compute the second componet in the Stokes SLP

        do n=1,ns
          DG(n) = grad1(2,n) + source(2,n)*grad2(2,n) + pot2(n)
        enddo


      endif



      end



