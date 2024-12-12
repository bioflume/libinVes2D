      subroutine stokesSLP(s1,s2,xs,ys,ns,u1,u2)
c     Evaluates the single-layer potential for Stokes' equation
c     (s1,s2) is the source strength of length ns
c     xs and ys are the source and target locations of length ns
c     (u1,u2) are the two components of the velocity field
      implicit real*8 (a-h,o-z)

      integer error
      real *8 s1(ns),s2(ns)
c     strength of Stokes SLP
      real *8 xs(ns),ys(ns)
c     location of the source/target points
      real *8 u1(ns),u2(ns)
c     x and y components of the velocity field
      
      real *8, allocatable :: source(:,:)
c     location of sources/targets
      real *8, allocatable :: charge(:)
c     charge strength of single-layer potential term
      real *8, allocatable :: dipstr(:),dipvec(:,:)
c     charge strength and direction of the 
c     double-layer potential term
      real *8, allocatable :: pot1(:),pot2(:)
c     room for two components that have to be summed
c     to form the velocity field

      allocate(source(2,ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(charge(ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(dipstr(ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(pot1(ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(pot2(ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(dipvec(2,ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
c     allocate memory for temporary variables

      twopi = 8.d0*datan(1.d0)

      iprec = 4 ! enough for 12 digits

      ifpot = 1 ! need the potential
      ifgrad = 0 ! don't need the gradient ... yet
      ifhess = 0 ! don't need the Hessian
      do i=1,ns
        source(1,i) = xs(i)
        source(2,i) = ys(i)
      enddo
c     set charge locations


c     START OF FORMING FIRST COMPONENET OF VELCOTIY
      ifcharge = 1 ! need a charge component
      ifdipole = 1 ! need a dipole componenet

      do i=1,ns
        charge(i) = -s1(i)
        dipstr(i) = xs(i)
        dipvec(1,i) = s1(i)
        dipvec(2,i) = s2(i)
      enddo

      call rfmm2dpart(ierr,iprec,ns,source,
     $   ifcharge,charge,ifdipole,dipstr,dipvec,
     $   ifpot,pot1,ifgrad,grad,ifhess,hess)
c     compute the first componet in the Stokes SLP
c      do i=1,5
c        print*,i,source(1,i),source(2,i),pot1(i)
c      enddo

      ifcharge = 0 ! don't need a charge component
      ifdipole = 1 ! need a dipole componenet

      do i=1,ns
        dipstr(i) = -1.d0
      enddo

      call rfmm2dpart(ierr,iprec,ns,source,
     $   ifcharge,charge,ifdipole,dipstr,dipvec,
     $   ifpot,pot2,ifgrad,grad,ifhess,hess)
c     compute the second componet in the Stokes SLP

      do i=1,ns
        u1(i) = pot1(i) + xs(i)*pot2(i)
        u1(i) = 5.d-1*u1(i)/twopi
      enddo
c     END OF FORMING FIRST COMPONENET OF VELCOTIY


c     START OF FORMING SECOND COMPONENET OF VELCOTIY
      ifcharge = 1 ! need a charge component
      ifdipole = 1 ! need a dipole componenet

      do i=1,ns
        charge(i) = -s2(i)
        dipstr(i) = ys(i)
        dipvec(1,i) = s1(i)
        dipvec(2,i) = s2(i)
      enddo

      call rfmm2dpart(ierr,iprec,ns,source,
     $   ifcharge,charge,ifdipole,dipstr,dipvec,
     $   ifpot,pot1,ifgrad,grad,ifhess,hess)
c     compute the first componet in the Stokes SLP

      ifcharge = 0 ! don't need a charge component
      ifdipole = 1 ! need a dipole componenet

      do i=1,ns
        dipstr(i) = -1.d0
      enddo

      call rfmm2dpart(ierr,iprec,ns,source,
     $   ifcharge,charge,ifdipole,dipstr,dipvec,
     $   ifpot,pot2,ifgrad,grad,ifhess,hess)
c     compute the second componet in the Stokes SLP

      do i=1,ns
        u2(i) = pot1(i) + ys(i)*pot2(i)
        u2(i) = 5.d-1*u2(i)/twopi
      enddo
c     END OF FORMING SECOND COMPONENET OF VELCOTIY




      end

