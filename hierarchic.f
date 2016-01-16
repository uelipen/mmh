c hierachic.f
c written Dec 13, 1999, by Ue-Li Pen, pen@cita.utoronto.ca
c fill in an array of numbers in a hierarchical fashion.
c
c allows 
c 
      subroutine parallelrand(def,iseed)
      implicit none
      include 'relaxgl.fi'
      real def(ng1,ng2,ng3)
      integer iseed(4)
c locals      
      integer nstrip,nblock
      parameter (nblock=64**3,nstrip=ng1*ng2*ng3/nblock)
      integer iseedl(4)
      real rseed(4,nstrip)
      integer i,idist,j

c     gaussian      
      idist=3
      
      if (nstrip .lt. 1) then
         call DLARNV(IDIST, ISEED, ng1*ng2*ng3, def )
      else
         if (nblock*nstrip .ne. ng1*ng2*ng3) then
            write(*,*) 'prand: nblock,nstrip,n=',nblock,nstrip
     &           ,ng1*ng2*ng3
            stop
         endif
c generate uniformly distributed offsets         
         call dlarnv(1,iseed,4*nstrip,rseed)
         write(*,*) 'prand nstrip=',nstrip
c$omp parallel do default(none) shared(def,rseed), private(i,j,iseedl)
c$omp&   ,firstprivate(iseed,idist)
         do i=1,nstrip
            do j=1,4
               iseedl(j)=iseed(j)+rseed(j,i)*4096
               iseedl(j)=mod(iseedl(j),4096)
            enddo
            iseedl(4)=iseedl(4)-mod(iseedl(4),2)+1
            iseedl(4)=mod(iseedl(4),4096)
            call dlarnv(idist,iseedl,nblock,def((i-1)*nblock+1,1,1))
         enddo
          write(*,*) 'done parallel dlarnv'
      endif
      return
      end

      subroutine hierarch(rho1,rho2,n1,n2,n3)
      implicit none
      integer n1,n2,n3
      real rho1(n1,n2,n3),rho2(n1,n2,n3)
c locals
      integer m1,m2,m3,i,j,k,ic,it,jt,kt, ic1, ic0

! not known to work for unequal n123
      m1=1
      m2=1
      m3=1

      ic1=0

 100  continue
      ic0=ic1
      m1=min(m1,n1)
      m2=min(m2,n2)
      m3=min(m3,n3)
      write(*,*) 'unpack, m123=',m1,m2,m3
c$omp parallel default(none) firstprivate(m1,m2,m3,n1,n2,n3,ic0)
c$omp& shared(rho1,rho2) reduction(+:ic1) private(i,j,k,it,jt,kt,ic)
      do k=1,m3
c$omp do
         do j=1,m2
            do i=1,m1
               if (i .gt. m1/2 .or. j .gt. m2/2 .or. k .gt. m3/2) then
                  ic1=ic1+1
                  if (k.gt.m3/2) then
                     ic=i+(j-1)*m1+(k-1)*m1*m2-m1*m2*m3/8
                  else if (j .gt. m2/2) then
                     ic=i+(j-1)*m1+(k-1)*m1*m2-m1*m2*k/4
                  else
                     ic=i-m1/2+(j-1)*m1/2+(k-1)*(m1*m2/4)*3
                  endif
		  ic=ic+ic0
c 		  if (ic .ne. ic1) write(*,'(9i8)') ic,ic1,i,j,k,m1,m2,m3
c the following line causes problems with SGI 7.3
!                  rho2(i,j,k)=rho1(ic,1,1)
                  it=mod(ic-1,n1)+1
                  jt=mod((ic-it)/n1,n2)+1
                  kt=mod((ic-it-(jt-1)*n1)/(n1*n2),n3)+1
                  rho2(i,j,k)=rho1(it,jt,kt)
               endif
            enddo
         enddo
      enddo
c$omp end parallel
      m1=m1*2
      m2=m2*2
      m3=m3*2
      if (m1 .le. n1 .or. m2 .le. n2 .or. m3 .lt. n3) goto 100
      write(*,*) 'done unpack'
c$omp parallel default(private) firstprivate(n1,n2,n3)
c$omp& shared(rho1,rho2)
      do k=1,n3
c$omp do
         do j=1,n2
            do i=1,n1
               rho1(i,j,k)=rho2(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end


      subroutine sinc(rho1,rho2,n1,n2,n3)
      implicit none
      integer n1,n2,n3
      real rho1(n1,n2,n3),rho2(n1,n2,n3)
c locals
      integer m1,m2,m3,i,j,k,ic,it
      

      m1=n1/2
 100  continue
c$omp parallel private(i,j,k,it), firstprivate(m1,n2,n3)
c$omp& shared(rho2,rho1)
      do k=1,n3
c$omp do
         do j=1,n2
            do i=1,m1
               it=2*(i-1)+1
               rho2(i,j,k)=rho1(it,j,k)+rho1(it+1,j,k)
               rho2(i+m1/2,j,k)=rho1(it,j,k)-rho1(it+1,j,k)
            enddo
         enddo
      enddo
c$omp end parallel
      m1=m1/2
      if (m1 .gt. 1) goto 100
      
c not done yet

      end

      subroutine copy(rho1,rho2,n1,n2,n3,m1,m2,m3)
      implicit none
      integer m1,m2,m3,n1,n2,n3
      real rho1(n1,n2,n3),rho2(n1,n2,n3)
c locals
      integer i,j,k
c$omp parallel private(i,j,k) firstprivate(m1,m2,m3)
c$omp& shared(rho1,rho2)
      do k=1,m3
c$omp do
         do j=1,m2
            do i=1,m1
               rho1(i,j,k)=rho2(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end
      
      subroutine isinc(rho1,rho2,n1,n2,n3)
c inverse sinc transform
      implicit none
      integer n1,n2,n3
      real rho1(n1,n2,n3),rho2(n1,n2,n3)
c locals
      integer m1,m2,m3,i,j,k,it

      m1=1
 100  continue
c$omp parallel default(private) firstprivate(m1,n2,n3)
c$omp& shared(rho1,rho2)
      do k=1,n3
c$omp do
         do j=1,n2
            do i=1,m1
               it=2*(i-1)+1
               rho2(it,j,k)=(rho1(i,j,k)+rho1(i+m1,j,k))/sqrt(2.)
               rho2(it+1,j,k)=(rho1(i,j,k)-rho1(i+m1,j,k))/sqrt(2.)
            enddo
         enddo
      enddo
c$omp end parallel
      m1=m1*2
      call copy(rho1,rho2,n1,n2,n3,m1,n2,n3)
      if (m1 .le. n1/2) goto 100
      m2=1
 200  continue
c$omp parallel default(private) firstprivate(m2,n1,n3)
c$omp& shared(rho1,rho2)
      do k=1,n3
c$omp do
         do j=1,m2
            do i=1,n1
               it=2*(j-1)+1
               rho2(i,it,k)=(rho1(i,j,k)+rho1(i,j+m2,k))/sqrt(2.)
               rho2(i,it+1,k)=(rho1(i,j,k)-rho1(i,j+m2,k))/sqrt(2.)
            enddo
         enddo
      enddo
c$omp end parallel
      m2=m2*2
      call copy(rho1,rho2,n1,n2,n3,n1,m2,n3)
      if (m2 .le. n2/2) goto 200
      m3=1
 300  continue
c$omp parallel default(private) firstprivate(m3,n1,n2)
c$omp& shared(rho1,rho2)
      do k=1,m3
c$omp do
         do j=1,n2
            do i=1,n1
               it=2*(k-1)+1
               rho2(i,j,it)=(rho1(i,j,k)+rho1(i,j,k+m3))/sqrt(2.)
               rho2(i,j,it+1)=(rho1(i,j,k)-rho1(i,j,k+m3))/sqrt(2.)
            enddo
         enddo
      enddo
c$omp end parallel
      m3=m3*2
      call copy(rho1,rho2,n1,n2,n3,n1,n2,m3)
      if (m3 .le. n3/2) goto 300
      return
      end

      
      subroutine testh
      parameter( n=4)
      real rho1(n,n,n),rho2(n,n,n)
      do i=1,n**3
         rho1(i,1,1)=0
      enddo
      rho1(1,1,1)=1
      call hierarch(rho1,rho2,n,n,n)
      write(*,'(4f5.0)') rho1
      call isinc(rho1,rho2,n,n,n)
      write(*,'(4f10.5)') rho1
      end






















