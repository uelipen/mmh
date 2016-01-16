c     these routines rely on arrays being linearly mapped
c     and that compilations is not optimized
c     so that array compression may be performed in-place
c      
cccccc
cccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc
c  Jacket routine for dxml fft3r call.
      subroutine fft3r1(a,b,c,d,n1,n2,n3)
c real-to-complex
      integer n1,n2,n3,i,j,k
      real a(n1,n2,n3+2)
      real b(n1+2,n2,n3)
      complex c(n1/2+1,n2,n3),d(n1/2,n2,n3+2)
      include 'relaxgl.fi'
      complex temp(ng2,ng3)

      if (ng2 .ne. n2 .or. ng3 .ne. n3) then
	write(*,*) 'generic fft: dimensions dont match',ng2,n2,ng3,n3
        stop
      endif
      chk=b(n1/2,n2/2,n3/2)
      do k=1,n3
         do j=1,n2
            do i=1,n1
               a(i,j,k)=b(i,j,k)
            enddo
         enddo
      enddo
      if (chk .ne. a(n1/2,n2/2,n3/2)) then
	write(*,*) 'remap failure in genericfft.f'
      write(*,*)'make sure that fft3r1 and fft3rinv1 are not optimized'
        stop
      endif

      call fft3rk(a,temp,n1,n2,n3)
      do k=n3,1,-1
         do j=n2,1,-1
            do i=n1/2,1,-1
               c(i,j,k)=d(i,j,k)
            enddo
         enddo
      enddo
      do k=1,n3
         do j=1,n2
            c(n1/2+1,j,k)=temp(j,k)
         enddo
      enddo
      return
      end
      subroutine fft3rinv1(a,b,c,d,n1,n2,n3)
      integer n1,n2,n3,i,j,k
      real a(n1,n2,n3+2),xn3
      real b(n1+2,n2,n3)
      complex c(n1/2+1,n2,n3),d(n1/2,n2,n3+2)
c locals
      include 'relaxgl.fi'
      complex temp(ng2,ng3)

      xn3=1./(n1*n2*n3)
c
      if (ng2 .ne. n2 .or. ng3 .ne. n3) then
	write(*,*) 'generic fft: dimensions dont match',ng2,n2,ng3,n3
        stop
      endif

      do k=1,n3
         do j=1,n2
            temp(j,k)=c(n1/2+1,j,k)
         enddo
      enddo

      do k=1,n3
         do j=1,n2
            do i=1,n1/2
               d(i,j,k)=c(i,j,k)
            enddo
         enddo
      enddo
      
      call fft3kr(d,temp,n1,n2,n3)

      do k=n3,1,-1
         do j=n2,1,-1
            do i=n1,1,-1
               b(i,j,k)=a(i,j,k)*xn3
            enddo
         enddo
      enddo
      do k=1,n3
         do j=1,n2
            do i=n1+1,n1+2
               b(i,j,k)=0
            enddo
         enddo
      enddo
      return
      end
