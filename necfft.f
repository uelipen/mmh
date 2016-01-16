c NEC version:
      subroutine fft3(a,n,isign)
c 3-d real-to-complex fft with normalization at the inverse transform
c isign=1 means forward transform
c isign=-1 means inverse transform
c
      complex a(n/2+1,n,n)

      if (isign .eq. 1) then
         call fft3r(a,n)
      else if ( isign .eq. -1) then
         call fft3rinv(a,n)
      else
         write(*,*) 'fft3: unknown argument: isign=',isign
         stop
      endif
      return
      end

**********************************************************************

c MLIB Man Pages
c 
c crc3ft, zrc3ft - real-to-complex three-dimensional FFT - full storage mode
c 
c INTEGER l1, l2, l3, ldz, mdz, iopt, ier
c COMPLEX z(ldz, mdz, l3)
c CALL CRC3FT (z, l1, l2, l3, ldz, mdz, iopt, ier)
c 
c INTEGER*4       l1, l2, l3, ldz, mdz, iopt, ier
c COMPLEX*16      z(ldz, mdz, l3)
c CALL ZRC3FT (z, l1, l2, l3, ldz, mdz, iopt, ier)
c 
c DESCRIPTION
c These subprograms compute either the forward real-to-complex or the inverse
c complex-to-real three-dimensional discrete Fourier transform using a radix
c 2-3-5 fast Fourier transform (FFT) algorithm optimized for real input or
c output.
c A pair of companion subprograms, SRC3FT and DRC3FT, performs a similar
c operation, but in a space-conserving manner.
c 
c Option flag:
c iopt    >=      0       Compute forward real-to-complex transform.
c iopt    <       0       Compute inverse complex-to-real transform.
 
**********************************************************************

      subroutine fft3r(a,n)

      integer n,i,j,k,ier
      real a(n+2,n,n),x(2,n,n,n)
c$omp parallel default(private) shared(n,x,a)
      do k=1,n
c$omp do
        do j=1,n
          do i=1,n
            x(1,i,j,k)=a(i,j,k)
            x(2,i,j,k)=0
          enddo
        enddo
      enddo
c$omp end parallel
      call crc3ft(x,n,n,n,n,n,1,ier)
c$omp parallel default(private) shared(n,x,a)
      do k=1,n
c$omp do
        do j=1,n
          do i=1,n/2+1
            a(2*i-1,j,k)=x(1,i,j,k)
            a(2*i  ,j,k)=x(2,i,j,k)
          enddo
        enddo
      enddo
c$omp end parallel
      return
      end

**********************************************************************

      subroutine fft3rinv(a,n)

      integer n,i,j,k,ier
      real a(2,n/2+1,n,n),x(2,n,n,n)
c$omp parallel default(private) shared(n,x,a)
      do k=1,n
c$omp do
        do j=1,n
          do i=1,n/2+1
            x(1,i,j,k)= a(1,i,j,k)
            x(2,i,j,k)= a(2,i,j,k)
          enddo
c$omp end parallel
          do i=n/2+1+1,n
            x(1,i,j,k)= a(1,n+2-i,j,k)
            x(2,i,j,k)=-a(2,n+2-i,j,k)
          enddo
        enddo
      enddo
      call crc3ft(x,n,n,n,n,n,-1,ier)
c$omp parallel default(private) shared(n,x,a)
      do k=1,n
c$omp do
        do j=1,n
          do i=1,n/2
            a(1,i,j,k)=x(1,2*i-1,j,k)
            a(2,i,j,k)=x(1,2*i  ,j,k)
          enddo
          a(1,n/2+1,j,k)=0
          a(2,n/2+1,j,k)=0
        enddo
      enddo
c$omp end parallel
      return
      end









