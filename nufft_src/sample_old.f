c     c Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
c     c Contact: greengard@cims.nyu.edu
c     c 
c     c This software is being released under a FreeBSD license
c     c (see license.txt in this directory). 
c     
      program testfft
      implicit none
c     
c     --- local variables
c     
      integer i,ier,iflag,j,k1,mx,ms,nj
      parameter (mx=10 000)
      real*8 xj(mx), sk(mx)
      real*8 err,eps,pi
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(mx),cj0(mx),cj1(mx)
      complex*16 fk0(mx),fk1(mx) 
      real*8, allocatable :: fw(:)

c     
c     --------------------------------------------------
c     create some test data
c     --------------------------------------------------
c      ms = 90
c      nj = 128
c      do k1 = -nj/2, (nj-1)/2
c         j = k1+nj/2+1
c         xj(j) = pi * dcos(-pi*j/nj)
c         cj(j) = dcmplx( dsin(pi*j/nj), dcos(pi*j/nj))
c      enddo

c
c              (ms-1)/2
c     cj(j) =    SUM      fk(k1) exp(-i k1 xj(j))  for j = 1,...,nj
c              k1= -ms/2                            
c
c xj: angular freqency (freq = 2 pi x)
c     ms- 

      ms=8
      nj=8 
      i=1
      do i = 1,nj
         xj(i) = 2*pi*real(i)/ms 
      enddo

      do j = 1,ms/2          
         fk0(j)= dcmplx(real(j+ms/2),0.0)
      enddo
      do j = ms/2+1,ms          
         fk0(j)= dcmplx(real(j-ms/2),0.0)
      enddo

      allocate ( fw(0:150) )
      do j = 1,ms
c         fw(2*j)=real(fk0(j))
c         fw(2*j+1)=imag(fk0(j))
         fw(2*j)=real(j)
         fw(2*j+1)=0.0
      enddo
      

      call zffti(nj,fw(150))
      call zfftf(nj,fw(0),fw(150))

c      stop

c      fk0(3)=1.0
c      fk0(4)=2.0
c      fk0(1)=3.0
c      fk0(2)=4.0

c     
c     --------------------------------------------------
c     start tests
c     --------------------------------------------------
c     
      iflag = -1
      print*,' Start 1D testing: ', ' nj =',nj, ' ms =',ms
      i = 6
c     do i = 1,4
      if (i.eq.1) eps=1d-4
      if (i.eq.2) eps=1d-8
      if (i.eq.3) eps=1d-12
      if (i.eq.4) eps=1d-16
c     extended/quad precision tests
      if (i.eq.5) eps=1d-20
      if (i.eq.6) eps=1d-24
      if (i.eq.7) eps=1d-28
      if (i.eq.8) eps=1d-32
      print*,' '
      print*,' Requested precision eps =',eps
      print*,' '
c     
c     -----------------------
c     call 1D Type1 method
c     -----------------------
c     
c      call dirft1d1(nj,xj,cj,iflag, ms,fk0)
c      print *, fk0

c      call nufft1d1f90(nj,xj,cj,iflag,eps, ms,fk1,ier)
c      call errcomp(fk0,fk1,ms,err)
c      print *,' ier = ',ier
c      print *,' type 1 error = ',err
c     
c     -----------------------
c     call 1D Type2 method
c     -----------------------
c     
c   in fk0, nj
c   out cj, ier   

      print *, xj(1:nj)
      print *, fk0(1:ms)

      call dirft1d2(nj,xj,cj0,iflag, ms,fk0,ier)
      print *, "DFT"
      print *, cj0(ms)
      do i=1,ms-1
         print *, cj0(i)
      enddo

      call nufft1d2f90(nj,xj,cj1,iflag, eps, ms,fk0,ier)
      print *, nj
      print *, cj1(ms)
      print *, "NFFT"
      print *, cj1(ms)
      do i=1,ms-1
         print *, cj1(i)
      enddo

      print *, "Their FFT"

      print *, fw(0)+fw(2*ms), fw(1)+fw(2*ms+1)
      do j=1,ms-1
         print *, fw(2*j), fw(2*j+1)
      enddo


      call errcomp(cj0,cj1,nj,err)
      print *,' ier = ',ier
      print *,' type 2 error = ',err
      stop
      end
c
c
c
c
c
      subroutine errcomp(fk0,fk1,n,err)
      implicit none
      integer k,n
      complex*16 fk0(n), fk1(n)
      real *8 salg,ealg,err
c
      ealg = 0d0
      salg = 0d0
      do k = 1, n
         ealg = ealg + cdabs(fk1(k)-fk0(k))**2
         salg = salg + cdabs(fk0(k))**2
      enddo
      err =sqrt(ealg/salg)
      return
      end
