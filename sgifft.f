c SGi version:
      subroutine fft3(a,n,isign)
c 3-d real-to-complex fft with normalization at the inverse transform
c isign=1 means forward transform
c isign=-1 means inverse transform
      complex*16 a(n/2+1,n,n)
      real table(100000),workspace(100000),scale
      integer omp_get_max_threads,iomaxt
c$    external omp_get_max_threads

        
      iomaxt=omp_get_max_threads()
      if (iomaxt .gt. 64) then
          call omp_set_num_threads(64)
      endif
      isys=0
      if (isign .eq. 1) then
	scale=1
      call dzfft3d(0, n, n, n, scale, a, n+2, n, 
     &       a,n/2+1,n ,table,work,isys)
        call dzfft3d(-isign, n, n, n, scale, a, n+2, n, 
     &       a,n/2+1,n ,table,work,isys)
      else
        scale=1./n**3
        call zdfft3d(0, n, n, n, scale, a, n/2+1, n, 
     &       a,n+2,n ,table,work,isys)
        call zdfft3d(-isign, n, n, n, scale, a, n/2+1, n, 
     &       a,n+2,n ,table,work,isys)
      endif

      if (iomaxt .gt. 64) then
          call omp_set_num_threads(iomaxt)
      endif
      return
      end


      subroutine sfft3(a,n,isign)
c 3-d real-to-complex fft with normalization at the inverse transform
c isign=1 means forward transform
c isign=-1 means inverse transform
      complex*8 a(n/2+1,n,n)
      real*4 workspace(1000000),table(100000)


      isys=0
      if (isign .eq. 1) then
        call scfft3d(0, n, n, n, scale, a, n+2, n,
     &       a,n/2+1,n ,table,work,isys)
        scale=1
        call scfft3d(-isign, n, n, n, scale, a, n+2, n,
     &       a,n/2+1,n ,table,work,isys)
      else
        call csfft3d(0, n, n, n, scale, a, n/2+1, n,
     &       a,n+2,n ,table,work,isys)
        scale=1./n**3
        call csfft3d(-isign, n, n, n, scale, a, n/2+1, n,
     &       a,n+2,n ,table,work,isys)
      endif

        

      return
      end
