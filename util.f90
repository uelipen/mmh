subroutine warray(a,n,nf)
  implicit none
  integer n,nf
  real*4 a(n)
  integer*2 iout(n),icompress(n)
  ! locals
  integer imax,imin,nicmp
  real amax,amin, atmp, aminl, amaxl


  if (.true.) then
	write(nf) a
  	return
  endif
  imax=2**15-1
  amax=maxval(abs(a))
! now get rid of zeros
  atmp=2*amax
  where (a .eq. 0) a=atmp
  amin=minval(abs(a))
  amin=max(amin,amax/1000000)
  if (amin .eq. amax) amin=amin/2
  aminl=log(amin)
  amaxl=log(amax)
  where (abs(a) < amin .or. a .eq. atmp) a=sign(amin,a)
  where (a > 0)
     iout=nint((log(a)-aminl)/(amaxl-aminl)*imax)
  elsewhere
     iout=-nint((log(-a)-aminl)/(amaxl-aminl)*imax)
  end where
     
  iout(2:)=iout(2:)-iout(:n-1)
!  call fcompress(iout,icompress,nicmp,n*2)
!  if (nicmp.lt.n .and. nicmp .gt. 0 ) then
!    write(nf) amin,amax,nicmp
!    write(nf) icompress(:nicmp)
!  else
    write(nf) amin,amax,n
    write(nf) iout
!  endif
!  write(*,*) 'warray: compression = ',2*n,nicmp
  return
end subroutine warray
  
