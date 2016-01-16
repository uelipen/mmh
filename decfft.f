c this file is not yet used.
c DEC DXML version:
      subroutine fft3(a,n,isign)
c 3-d real-to-complex fft with normalization at the inverse transform
c isign=1 means forward transform
c isign=-1 means inverse transform
      complex a(n/2+1,n,n)
      integer sfft_3d, dfft_3d
      external sfft_3d, dfft_3d
        
      if (isign .eq. 1) then
         ierr=dfft_3d('R','C','F',a,a,n,n,n,n+2,n,1,1,1)
      else
         ierr=dfft_3d('C','R','B',a,a,n,n,n,n+2,n,1,1,1)
      endif
c I think dxml normalizes on inverse transform
      if (ierr .ne. 0) then
         write(*,*) 'forward fft failed, ier=',ierr
         write(*,*) 'isign,n=',isign,n
         stop
      endif
      return
      end


      subroutine sfft3(a,n,isign)
c 3-d real-to-complex fft with normalization at the inverse transform
c isign=1 means forward transform
c isign=-1 means inverse transform
      complex*8 a(n/2+1,n,n)
      integer sfft_3d, dfft_3d
      external sfft_3d, dfft_3d
        
      if (isign .eq. 1) then
         ierr=sfft_3d('R','C','F',a,a,n,n,n,n+2,n,1,1,1)
      else
         ierr=sfft_3d('C','R','B',a,a,n,n,n,n+2,n,1,1,1)
      endif
c I think dxml normalizes on inverse transform
      if (ierr .ne. 0) then
         write(*,*) 'forward fft failed, ier=',ierr
         write(*,*) 'isign,n=',isign,n
         stop
      endif
      return
      end



