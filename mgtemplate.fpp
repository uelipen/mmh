c-*- Fortran -*-	
c
c global constant NMG must be defined before inclusion of this file.
c

      subroutine mg NMG(u,rhs,defpot,ubig,dx,nu,iopt)
      implicit none
      integer iopt,nu
      real u(NMG,NMG,NMG),rhs(NMG,NMG,NMG),defpot(NMG,NMG,NMG)
     &,ubig(nu,NMG,NMG,NMG),dx
cdir$ shared *u(:block,:block,:),*defpot(:block,:block,:)
cdir$ shared *rhs(:block,:block,:)      
c locals
      integer isolver,imgalg,irelaxopt,nrelax,i
#define N2 NMG/2
      real uhalf(N2,N2,N2),rhshalf(N2,N2,N2)
     &,defhalf(N2,N2,N2),usmall(1,N2,N2,N2)
cdir$ shared uhalf(:block,:block,:),defhalf(:block,:block,:)
cdir$ shared rhshalf(:block,:block,:),usmall(:,:block,:block,:)
      real relaxerr,r1


      isolver=mod(iopt,8)

#ifdef _ALPHA
c first touch placement of the subarrays
ccdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(uhalf,defhalf,rhshalf,usmall)
      call szero(uhalf,N2,N2,N2)
      call szero(defhalf,N2,N2,N2)
      call szero(rhshalf,N2,N2,N2)
      call szero(usmall,N2,N2,N2)
#endif
      nrelax=4*NG/NMG
      irelaxopt=isolver*8+3
c calculate residual on half grid      
      call relax(u,rhshalf,rhs,defpot,ubig,r1,dx,NMG,NMG,NMG
     &     ,nu,nrelax,irelaxopt) 

c      goto 59
      call restrict(defhalf,defpot,NMG,NMG,NMG)
      call restrictu(usmall,ubig,NMG,NMG,NMG,nu)
      call szero(uhalf,N2,N2,N2)

#if NMG == 8
      call mg4(uhalf,rhshalf,defhalf,usmall,dx*2,1,iopt)
#endif
#if NMG == 16      
      call mg8(uhalf,rhshalf,defhalf,usmall,dx*2,1,iopt)
#endif
#if NMG == 32      
      call mg16(uhalf,rhshalf,defhalf,usmall,dx*2,1,iopt)
#endif
#if NMG == 64
      call mg32(uhalf,rhshalf,defhalf,usmall,dx*2,1,iopt)
#endif
#if NMG == 128    
      call mg64(uhalf,rhshalf,defhalf,usmall,dx*2,1,iopt)
#endif
#if NMG == 256     
      call mg128(uhalf,rhshalf,defhalf,usmall,dx*2,1,iopt)
#endif
#if NMG == 512     
      call mg256(uhalf,rhshalf,defhalf,usmall,dx*2,1,iopt)
#endif
      if (NMG .gt. 512) then
         write(*,*) 'mg: invalid dimension, NMG=',NMG
         stop
      endif
      
      call inject(u,uhalf,NMG,NMG,NMG)
59    continue
      irelaxopt=isolver*8+2
      call relax(u,rhshalf,rhs,defpot,ubig,relaxerr,dx,NMG,NMG,NMG
     &     ,nu,nrelax,irelaxopt) 
      if (NMG .eq. NG) then
         write(*,*) 'relaxerr=',relaxerr
      endif
      if (NMG .ne. NG .and. relaxerr .gt. r1 ) then
        write(*,73) NMG,r1,relaxerr
73      format(' mg: relaxation diverged, level=',i3,' r1,2=',2e12.4)
c prolong resets u: forget about the current level relaxation!
c        call prolong(u,uhalf,NMG,NMG,NMG)
        call szero(u,NMG,NMG,NMG)
      else
        call zerosum(u,NMG,NMG,NMG)
      endif
      return
      end
