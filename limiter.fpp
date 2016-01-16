c -*- Fortran -*- 77 file limiter.fpp
c**************************************************************************
c File: limiter.fpp
c extracted from relaxing.fpp 11/2/97 
c to implement a better homogeneous grid limiter.

c
c the fundamental adjustable parameters are the number of times
c that nsmooth is called.  The first ambiguity is how accurately we
c want to track the density field.  \rho\rg is smoothed, if it
c smoothes too little, the limiters start failing.
c the second is the number of neighboring cells that are subject to the
c limiter, i.e. the number of times that routine spreadmax is called.
c the third is the nsmooth of defp, which causes the grid to only respond
c to long wavelength perturbations.

c
c*********************************************************************
      subroutine calcdefp(defp,tmp,tmp2,def,u,dtold1,dtaumesh,nfluid)
c assume that u(1,:,:,:) is the density field that we wish to keep
c constant.
      implicit none
      include 'relaxgl.fi'
      integer nfluid
      real defp(ng1,ng2,ng3), def(ng1,ng2,ng3),u(nfluid,ng1,ng2,ng3)
     &      ,tmp2(ng1+2,ng2,ng3)
     &      ,tmp(ng1,ng2,ng3), dtold1, rhomin, rhomax,rms,dtaumesh

c we want to recycle the tmp array to be of two different sizes
c (not any more, that caused too much grief, better to be clean)
c but note that algorithmically, tmp and tmp2 might be equivalenced.      
       call calcdefp1(defp,tmp,tmp2,def,u,dtold1,dtaumesh,nfluid)
c the distribution of the second tmp can be problematic.
      return
      end
      
      subroutine calcdefp1(defp,tmp,defplim,def,u,dtold1,dtaumesh
     &     ,nfluid)
c assume that u(1,:,:,:) is the density field that we wish to keep
c constant.
      implicit none
      include 'relaxgl.fi'
      integer nfluid
      real defp(ng1,ng2,ng3), def(ng1,ng2,ng3),u(nfluid,ng1,ng2,ng3)
     &      ,tmp(ng1,ng2,ng3), dtold1, rhomin, rhomax,rms,dtaumesh
      real defplim(ng1+2,ng2,ng3)
c
c we need a full sized array as local storage to solve the poisson
c equation.
c
c locals
      integer i,j,k,ip,im,jp,jm,kp,km,kmm,ko,koo,kl,ncompress,nexpand
     &     ,ncold,ii
      real umean, baj(2,7,ng1,ng2), phixx,phiyy,phizz,phixy,phixz
     &     ,phiyz,a11,a12,a13,a23,a33,b11,b12,b13,b23,b33,b22,a22
     &     ,det,s1m,s2m,s3m,s1p,s2p,s3p, compression, eflux, cmax
     &     ,c11,c12,c13,c22,c23,c33, cmpmax, detmax, detmin,pi
     &     ,d11,d22,a,b,am,theta,trace,r1,r2,r3, rhol, rhonn, dtold
     &		,rho2,xlum,temp,temp1,t1,xlow,xhigh, cmpmaxglobal
     &	, expansion, fluxmax, fluxmin, tlim, det1, det2,xk,yk,zk
     &     ,akernel,threshold
cmydist baj(*,*,*,block)
      parameter(pi=3.14159265358979323846264338328)
      logical firsttime
      data firsttime /.true./
      include 'globalpa.fi'
c for adaptive initial conditions:
      real*4 densinit(ng1,ng2,ng3)
cmydist densinit(*,block,*)
      common /densi/densinit
c each component is a real number specifying the initial \rho \rg in each
c cell

#ifdef COLD
      include 'cold.fi'
#else
      real cold
! to trick OMP 
#endif
#ifdef GMETRIC
      include 'gmetric.fi'
#endif

      dtold=dtold1
      cmax=0
      
c the maximal linear compression allowed:
      cmpmaxglobal=compressmax
      if (cmpmaxglobal .lt. 0) then
         tmp(1,1,1)=0
         tmp(2,1,1)=0
         tmp(3,1,1)=0
         do i=1,ng1*ng2*ng3
            tmp(1,1,1)=tmp(1,1,1)+u(1,i,1,1)**2
         temp1=max(abs(u(1,i,1,1)),0.001)
            temp1=
     & abs(u(5,i,1,1)-(u(2,i,1,1)**2+u(3,i,1,1)**2+u(4,i,1,1)**2
     &         )/temp1/2)/temp1
               tmp(3,1,1)=tmp(3,1,1)+temp1
               xlum=u(1,i,1,1)**2*sqrt(temp1)
            tmp(2,1,1)=tmp(2,1,1)+xlum
         enddo
         return
      endif

      if (firsttime) then
         firsttime = .false.
         write(*,*) 'calcdefp: cmpmax=',cmpmaxglobal
         if ( densinit(1,1,1) .eq. -4321 ) then
            write(*,*) 'density init field not defined, setting to 1'
           
c$omp parallel default(private) shared(densinit)
            do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
               do j=1,ng2
                  do i=1,ng1
                     densinit(i,j,k)=1
                  enddo
               enddo
            enddo
c$omp end parallel
         endif
      endif
      if (dtold .lt. 1.e-8) then
         write(*,*) 'calcdefp: warning: dtold=',dtold,' setting to 1'
         dtold=1
      endif
c 1. approach the constant smoothed density field after ca 20 steps
c*$* assert do(serial)
c$omp parallel default(private) shared(tmp,densinit,u,dtold) 
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               tmp(i,j,k)=(densinit(i,j,k)-u(1,i,j,k))
            enddo
         enddo
      enddo
c short note on the constant density algorithm:
c in the absence of velocities, the continuity equation should be
c \dot{\rho\rg} + \div \rho\rg e \grad \defp = 0
c but we dropped the \rho\rg inside the divergence.  If that was a
c slowly varying quantity, we could just say
c \dot{\rho\rg} + \rho\rg \div e \grad \defp = 0
c and the change is now dimensionless \dot{\rho\rg}/(\rho\rg)=\div e ...
c the sign of the change is thus binary.
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               t1=tmp(i,j,k)
               xlow=abs(u(1,i,j,k))*0.2
               xlow=0
               xhigh=2
c if we set xhigh too large, the tidal fields of the deformation can
c lead to deformations which the grid limiter cannot catch.
               if (abs(t1) .lt. xlow) t1=0
               if (abs(t1) .gt. xhigh) t1=sign(xhigh,t1)
               tmp(i,j,k)=t1/30/dtold
            enddo
         enddo
      enddo
c$omp end parallel

      do i=1,3
         call nsmooth(tmp)
      enddo
    
      koo=1
cfpp$ skip
      do k=1,2
         koo=3-koo
         ko=3-koo
         kp=mod(k,ng3)+1
         km=mod(k+ng3-2,ng3)+1
         kmm=mod(km+ng3-2,ng3)+1
c$omp parallel default(private) firstprivate(k,koo,ko,kp,km,kmm)
c$omp& shared(def,baj)  
c$dir prefer_parallel_ext
c$omp do
         do j=1,ng2
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
               phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
               phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
               phixy=(def(ip,jp,k)-def(im,jp,k)
     &              -def(ip,jm,k)+def(im,jm,k))/4
               phiyz=(def(i,jp,kp)-def(i,jp,km)
     &              -def(i,jm,kp)+def(i,jm,km))/4
               phixz=(def(ip,j,kp)-def(im,j,kp)
     &              -def(ip,j,km)+def(im,j,km))/4
               a11=1+phixx
               a12=phixy
               a13=phixz
               a22=1+phiyy
               a23=phiyz
               a33=1+phizz
               det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &                 -a12**2*a33
               b11=(a22*a33-a23**2)/det
               b12=(a13*a23-a12*a33)/det
               b13=(a12*a23-a13*a22)/det
               b22=(a11*a33-a13**2)/det
               b23=(a12*a13-a11*a23)/det
               b33=(a11*a22-a12**2)/det
               baj(koo,1,i,j)=b11
               baj(koo,2,i,j)=b12
               baj(koo,3,i,j)=b13
               baj(koo,4,i,j)=b22
               baj(koo,5,i,j)=b23
               baj(koo,6,i,j)=b33
               baj(koo,7,i,j)=det
            enddo
         enddo
c$omp end parallel
      enddo
      ncompress=0
      nexpand=0
      detmax=0
      detmin=1
      rho2=0
      temp=0
      xlum=0
      rms=0
      fluxmax=0
      fluxmin=0
cfpp$ nodepchk r
c$omp parallel default(private) firstprivate(koo,dtold,dtaumesh
c$omp& ,cmpmaxglobal)
c$omp& shared(def,u,baj,densinit,tmp,cold) reduction(max:fluxmax,detmax) 
c$omp& reduction(min:fluxmin,detmin)  
c$omp& reduction(+:nexpand,xlum,temp,rho2)
      do kl=3,ng3+2
         koo=3-koo
         ko=3-koo
         k=mod(kl+ng3-1,ng3)+1
         kp=mod(k,ng3)+1
         km=mod(k+ng3-2,ng3)+1
         kmm=mod(km+ng3-2,ng3)+1
c$dir no_recurrence, force_parallel_ext
c*$* assert no recurrence(baj)         
c 7.2 f77
c*$* assert do(concurrent)
!cdir nodep
c$omp do
         do j=1,ng2
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
c$dir no_recurrence, force_vector 
!cdir nodep
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
#ifdef GMETRIC
               det=gbaj(7,i,j,k)
               b11=gbaj(1,i,j,k)
               b12=gbaj(2,i,j,k)
               b13=gbaj(3,i,j,k)
               b22=gbaj(4,i,j,k)
               b23=gbaj(5,i,j,k)
               b33=gbaj(6,i,j,k)
#else
               phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
               phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
               phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
               phixy=(def(ip,jp,k)-def(im,jp,k)
     &              -def(ip,jm,k)+def(im,jm,k))/4
               phiyz=(def(i,jp,kp)-def(i,jp,km)
     &              -def(i,jm,kp)+def(i,jm,km))/4
               phixz=(def(ip,j,kp)-def(im,j,kp)
     &              -def(ip,j,km)+def(im,j,km))/4
               a11=1+phixx
               a12=phixy
               a13=phixz
               a22=1+phiyy
               a23=phiyz
               a33=1+phizz
               det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &                 -a12**2*a33
               b11=(a22*a33-a23**2)/det
               b12=(a13*a23-a12*a33)/det
               b13=(a12*a23-a13*a22)/det
               b22=(a11*a33-a13**2)/det
               b23=(a12*a13-a11*a23)/det
               b33=(a11*a22-a12**2)/det
#endif
               temp1=max(abs(u(1,i,1,1)),0.001)
               temp1=
     & abs(u(5,i,j,km)-(u(2,i,j,km)**2+u(3,i,j,km)**2+u(4,i,j,km)**2
     &         )/temp1/2)/temp1
               temp=temp+temp1*det
               xlum=xlum+u(1,i,j,km)**2/det*sqrt(temp1)
               rho2=rho2+u(1,i,j,km)**2/det
#ifdef COLD
               rhol= 1.5*u(1,i,j,km)/baj(ko,7,i,j)
c if any of the six neighboring cells is more than twice as dense,
c revert to TVD               
               rhonn=max(u(1,i,j,kmm)/baj(koo,7,i,j), u(1,i,j,k)/det
     &         ,u(1,im,j,km)/baj(ko,7,im,j),u(1,ip,j,km)/baj(ko,7,ip,j)
     &         ,u(1,i,jm,km)/baj(ko,7,i,jm),u(1,i,jp,km)/baj(ko,7,i,jp)
     &              )
               if ( rhonn .gt. rhol ) then
                  cold(i,j,km)=.false.
               endif
#endif
               
c while the exact expression asks for \rho\rg e v, we ignore the \rho\rg,
c which for potential flow causes no errors.
               s1p=(baj(ko,1,ip,j)*u(2,ip,j,km)
     &              +baj(ko,2,ip,j)*u(3,ip,j,km)
     &              +baj(ko,3,ip,j)*u(4,ip,j,km))/u(1,ip,j,km)
               s2p=(baj(ko,2,i,jp)*u(2,i,jp,km)
     &              +baj(ko,4,i,jp)*u(3,i,jp,km)
     &              +baj(ko,5,i,jp)*u(4,i,jp,km))/u(1,i,jp,km)
               s3p=(b13*u(2,i,j,k)
     &              +b23*u(3,i,j,k)
     &              +b33*u(4,i,j,k))/u(1,i,j,k)
               s1m=(baj(ko,1,im,j)*u(2,im,j,km)
     &              +baj(ko,2,im,j)*u(3,im,j,km)
     &              +baj(ko,3,im,j)*u(4,im,j,km))/u(1,im,j,km)
               s2m=(baj(ko,2,i,jm)*u(2,i,jm,km)
     &              +baj(ko,4,i,jm)*u(3,i,jm,km)
     &              +baj(ko,5,i,jm)*u(4,i,jm,km))/u(1,i,jm,km)
               s3m=(baj(koo,3,i,j)*u(2,i,j,kmm)
     &              +baj(koo,5,i,j)*u(3,i,j,kmm)
     &              +baj(koo,6,i,j)*u(4,i,j,kmm))/u(1,i,j,kmm)
               eflux=(s1p+s2p+s3p-s1m-s2m-s3m)/2
               eflux=eflux+tmp(i,j,km)

               baj(koo,1,i,j)=b11
               baj(koo,2,i,j)=b12
               baj(koo,3,i,j)=b13
               baj(koo,4,i,j)=b22
               baj(koo,5,i,j)=b23
               baj(koo,6,i,j)=b33
               baj(koo,7,i,j)=det
               if (det .gt. 10*densinit(i,j,km)) then
                  eflux=-2*abs(eflux)
                  nexpand=nexpand+1
               endif
               c11=1+phixx
               c12=phixy
               c13=phixz
               c22=1+phiyy
               c23=phiyz
               c33=1+phizz
c     now apply compression and expansion limiters
               trace=c11+c22+c33
               d11=c11-trace/3
               d22=c22-trace/3
               a=-(d11**2+c12**2+c13**2+d11*d22+d22**2+c23**2)
               b= - d11*c12**2 + d11**2*d22 - c12**2*d22 + c13**2*d22 + 
     &              d11*d22**2 - 2*c12*c13*c23 + d11*c23**2
               if (abs(a) .lt. 1.e-20) a=1.e-20
               am=2*sqrt(abs(-a)/3)
               theta=acos(max(-1.,min(1.,3*b/a/am)))/3
               r1=am*cos(theta)
               r2=am*cos(theta+2*pi/3.)
               r3=am*cos(theta+4*pi/3.)
               compression=1/min(r1+trace/3,r2+trace/3,r3+trace/3)
               expansion=1/max(r1+trace/3,r2+trace/3,r3+trace/3)
c     we compress to some multiple of the original density
               cmpmax=cmpmaxglobal*densinit(i,j,k)**(-1./3.)
c     dont allow low density regions to compress: they screw up the
c     grid regularity.
c**** it is not clear that this is necessary or usefull
c               cmpmax=max(1.5,min(cmpmax, cmpmax*expansion
c     &              *densinit(i,j,k)**(-1./3.)  ))
c****
c     cmpmax is 1/dx
               tlim=3*(1/cmpmax-1/compression)/dtold
               if (tlim .gt. 0) eflux=max(0.,eflux)

               detmax=max(detmax,det)
               detmin=min(detmin,det)
               eflux=min(1./dtaumesh,eflux)
               eflux=max(-1./dtaumesh,eflux)
               tmp(i,j,km)=eflux
               fluxmax=max(fluxmax,eflux*dtold)
               fluxmin=min(fluxmin,eflux*dtold)
            enddo
         enddo
      enddo
c$omp end parallel
c the limiters cause the mean to be non-zero      

#if (DEBUG>0) 
      write(*,100) fluxmax, fluxmin,dtold
 100  format('calcdefp: fluxmax,min=',2e10.2,' dtold=',e10.2)
#endif

      do i=1,1
         call nsmooth(tmp)
      enddo

c$omp parallel default(private) shared(defp)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               defp(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel
c multigrid will call rgzerosum...
c
c the +32 indicates to use the current value as initial guess.
c 2 indicates flow solver.
      call multigrid(defp,tmp,def,u,1.,ng1,ng2,ng3,nfluid,1,2)
      do i=1,0
      call multigrid(defp,tmp,def,u,1.,ng1,ng2,ng3,nfluid,1,2+32)
      enddo

c we are done with tmp for now, its storage will be recycled for defplim

c smoothing with only 1 can crash a 32^3 run on modi4
c
c  the limiter algorithm assumes that the grid is actually smooth.
c if it is not, the grid gradients can lead to big errors.
c The simplest fix is to smooth the grid sufficiently.
c
      do i=1,2
         call nsmooth(defp)
      enddo
 
c$omp parallel default(private) firstprivate(dtold,cmpmaxglobal)
c$omp& shared(defp,densinit,def,defplim) reduction(+:ncompress)
c$omp&  reduction(max:cmax) 
      do k=1,ng3
         kp=mod(k,ng3)+1
         km=mod(k+ng3-2,ng3)+1
c$dir prefer_parallel_ext
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               phixx=(defp(ip,j,k)-2*defp(i,j,k)+defp(im,j,k))*dtold
               phiyy=(defp(i,jp,k)-2*defp(i,j,k)+defp(i,jm,k))*dtold
               phizz=(defp(i,j,kp)-2*defp(i,j,k)+defp(i,j,km))*dtold
               phixy=(defp(ip,jp,k)-defp(im,jp,k)
     &              -defp(ip,jm,k)+defp(im,jm,k))/4*dtold
               phiyz=(defp(i,jp,kp)-defp(i,jp,km)
     &              -defp(i,jm,kp)+defp(i,jm,km))/4*dtold
               phixz=(defp(ip,j,kp)-defp(im,j,kp)
     &              -defp(ip,j,km)+defp(im,j,km))/4*dtold
               phixx=phixx+(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
               phiyy=phiyy+(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
               phizz=phizz+(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
               phixy=phixy+(def(ip,jp,k)-def(im,jp,k)
     &              -def(ip,jm,k)+def(im,jm,k))/4
               phiyz=phiyz+(def(i,jp,kp)-def(i,jp,km)
     &                 -def(i,jm,kp)+def(i,jm,km))/4
               phixz=phixz+(def(ip,j,kp)-def(im,j,kp)
     &                 -def(ip,j,km)+def(im,j,km))/4
               c11=1+phixx
               c12=phixy
               c13=phixz
               c22=1+phiyy
               c23=phiyz
               c33=1+phizz
c     now apply compression and expansion limiters
               trace=c11+c22+c33
               d11=c11-trace/3
               d22=c22-trace/3
               a=-(d11**2+c12**2+c13**2+d11*d22+d22**2+c23**2)
               b= - d11*c12**2 + d11**2*d22 - c12**2*d22 + c13**2*d22 + 
     &              d11*d22**2 - 2*c12*c13*c23 + d11*c23**2
               if (abs(a) .lt. 1.e-20) a=1.e-20
               am=2*sqrt(abs(-a)/3)
               theta=acos(max(-1.,min(1.,3*b/a/am)))/3
               r1=am*cos(theta)
               r2=am*cos(theta+2*pi/3.)
               r3=am*cos(theta+4*pi/3.)
               compression=1/min(r1+trace/3,r2+trace/3,r3+trace/3)
               expansion=1/max(r1+trace/3,r2+trace/3,r3+trace/3)
c     we compress to some multiple of the original density
               cmpmax=cmpmaxglobal*densinit(i,j,k)**(-1./3.)
c     dont allow low density regions to compress: they screw up the
c     grid regularity.
c**** it is not clear that this is necessary or usefull
c               cmpmax=max(1.5,min(cmpmax, cmpmax*expansion
c     &              *densinit(i,j,k)**(-1./3.)  ))
c****
c     cmpmax is 1/dx
               tlim=3*(1/cmpmax-1/compression)/dtold
c     0.5 means that we want to do the correction over two time steps.
c     the three is because the other two dimensions do expand as well.
               threshold=6
#if (DEBUG > 0) 
               threshold=1.5
               phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
               phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
               phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
               phixy=(def(ip,jp,k)-def(im,jp,k)
     &              -def(ip,jm,k)+def(im,jm,k))/4
               phiyz=(def(i,jp,kp)-def(i,jp,km)
     &              -def(i,jm,kp)+def(i,jm,km))/4
               phixz=(def(ip,j,kp)-def(im,j,kp)
     &              -def(ip,j,km)+def(im,j,km))/4
               c11=1+phixx
               c12=phixy
               c13=phixz
               c22=1+phiyy
               c23=phiyz
               c33=1+phizz
c     now apply compression and expansion limiters
               trace=c11+c22+c33
               d11=c11-trace/3
               d22=c22-trace/3
               a=-(d11**2+c12**2+c13**2+d11*d22+d22**2+c23**2)
               b= - d11*c12**2 + d11**2*d22 - c12**2*d22 + c13**2*d22 + 
     &              d11*d22**2 - 2*c12*c13*c23 + d11*c23**2
               if (abs(a) .lt. 1.e-20) a=1.e-20
               am=2*sqrt(abs(-a)/3)
               theta=acos(max(-1.,min(1.,3*b/a/am)))/3
               r1=am*cos(theta)
               r2=am*cos(theta+2*pi/3.)
               r3=am*cos(theta+4*pi/3.)
               compression=1/min(r1+trace/3,r2+trace/3,r3+trace/3)
               expansion=1/max(r1+trace/3,r2+trace/3,r3+trace/3)
#endif
               cmax=max(cmax,compression)
#if (DEBUG >1)
               if (compression/cmpmax .gt. threshold) then
            write(*,30)i,j,k,compression,expansion,temp1,tmp(i,j,k),tlim
     &              ,u(1,i,j,k)
 30               format('i=',3i3,'cmn,rc,fl,tl',6e9.2)
               endif
#endif
               defplim(i,j,k)=0.
               if (tlim .gt. 0) then
                  ncompress=ncompress+1
c     we have 1+phixx=1/compression so phixx=1/compression-1
c     if dphixx = 1/cmpmax-1/compression,
                  defplim(i,j,k)=tlim
               endif
            enddo
         enddo
      enddo
c$omp end parallel

      if (ncompress .gt. 0) then
c$omp parallel default(private) shared(defplim)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=ng1+1,ng1+2
               defplim(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel

c we want n steps to fall off to 1/3:
c x^n=1/3 => x=(1/3)^(1/n)
      j=4
      temp1=3.**(-1./j)
      do i=1,j
c spreadmax seems to lead to anisotropic configurations...
c         call spreadmax(defplim,temp1,ng1+2)
      enddo


c we handle the limiter separately using FFTs:
      call fft3(defplim,ng1,1)
c$omp parallel default(private) shared(defplim)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1+2
               ii=(i-1)/2
               xk=2*pi*ii/ng1
               yk=2*pi*(j-1.)/ng2
               zk=2*pi*(k-1.)/ng3
               akernel=2*cos(xk)+2*cos(yk)+2*cos(zk)-6
               if (abs(akernel).lt. 1e-10) akernel=1
               defplim(i,j,k)=defplim(i,j,k)/akernel
            enddo
         enddo
      enddo
c$omp end parallel
      defplim(1,1,1)=0
      call fft3(defplim,ng1,-1)



c one more possibility is to compute the RHS of defplim
c according to what actually comes out of \nabla defp

c$omp parallel default(private)
c$omp& shared(defp,defplim)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               defp(i,j,k)=defp(i,j,k)+defplim(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel
      endif  ! if tlim>0

      rhomin=u(1,1,1,1)
      rhomax=rhomin
      ncold=0
c*$* assert do(serial)
c$omp parallel default(private)
c$omp& shared(u,tmp,cold) reduction(+:rms) reduction(max:rhomax)
c$omp& reduction(min:rhomin)
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               rhomin=min(rhomin,u(1,i,j,k))
               rhomax=max(rhomax,u(1,i,j,k))
               rms=rms+(u(1,i,j,k)-1)**2
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               tmp(i,j,k)=-1
#ifdef COLD
c this apparently complicated expression means the following:
c if both neighboring cells along any one axis are hot, set cell to hot.
               if ( .not. ((cold(i,j,km) .or. cold(i,j,kp)) .and.
     &              (cold(i,jp,k) .or. cold(i,jm,k)) .and.
     &              (cold(im,j,k) .or. cold(ip,j,k))) ) 
     &              tmp(i,j,k)=1
#endif
               
            enddo
         enddo
      enddo
c$omp end parallel
#ifdef COLD
c*$* assert do(serial)
c$omp parallel default(private) shared(tmp,cold) reduction(+:ncold)
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               if (tmp(i,j,k) .gt. 0) cold(i,j,k)=.false.
               if (cold(i,j,k)) ncold=ncold+1
            enddo
         enddo
      enddo
c$omp end parallel
               
      if (ncold .gt. 0) then
         write(*,91) float(ncold)/(ng1*ng2*ng3)
      endif
 91   format(' cold fraction: ',f20.10)
#endif     
      rms=sqrt(rms/(ng1*ng2*ng3))
      write(*,10)cmax,detmin,detmax,rhomin,rhomax,rms
10    format('cmax=',f8.5,' rgmin,max=',e12.5,f6.2,' rhomnx=',2f7.2
     &           , ' rms=',f5.2)


      if (ncompress .gt. 0) then
         if (nexpand .gt. 0) then
            write(*,*)'compression limiter x '
     &                ,ncompress,' expansion x ',nexpand
         else
            write(*,*)'compression limiter applied ',ncompress,' times'
         endif
      endif

      tmp(1,1,1)=rho2
      tmp(2,1,1)=xlum
      tmp(3,1,1)=temp
      return
      end


c*********************************************************************
      subroutine nsmooth(a)
      implicit none
      include 'relaxgl.fi'
      real a(ng1,ng2,ng3)
c locals
      real v(ng1,ng2,2), topa(ng1,ng2)
cmydist v(*,block,*),topa(*,block)
      integer i,j,k,ip,im,jp,jm,kp,km,kv,ko

c$omp parallel default(private) shared(a,v,topa)
      kv=1
cfpp$ skip
c*$* assert do(serial)      
      do k=1,2
         ko=kv
         kv=mod(k,2)+1
c$dir force_parallel
c*$* assert do(concurrent)
c$omp do      
         do j=1,ng2
c$dir force_vector
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               v(i,j,kv)=a(i,j,k)/8+(a(ip,j,k)+a(im,j,k)
     &              +a(i,jp,k)+a(i,jm,k)+a(i,j,kp)+a(i,j,km))/16
     &  +(a(ip,jp,k)+a(ip,jm,k)+a(im,jp,k)+a(im,jm,k)+a(ip,j,kp)
     &   +a(ip,j,km)+a(im,j,kp)+a(im,j,km)+a(i,jp,kp)+a(i,jp,km)
     &   +a(i,jm,kp)+a(i,jm,km))/32
     &  +(a(ip,jp,kp)+a(ip,jp,km)+a(ip,jm,kp)+a(ip,jm,km)
     &   +a(im,jp,kp)+a(im,jp,km)+a(im,jm,kp)+a(im,jm,km))/64
            enddo
         enddo
      enddo


c$dir force_parallel
c$omp do
      do j=1,ng2
c$dir force_vector
         do i=1,ng1
            topa(i,j)=v(i,j,ko)
         enddo
      enddo


c$dir scalar      
c*$* assert do(serial)      
      do k=3,ng3
         ko=kv
         kv=mod(k,2)+1
c$dir force_parallel
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
c$dir force_vector
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               v(i,j,kv)=a(i,j,k)/8+(a(ip,j,k)+a(im,j,k)
     &  +a(i,jp,k)+a(i,jm,k)+a(i,j,kp)+a(i,j,km))/16
     &  +(a(ip,jp,k)+a(ip,jm,k)+a(im,jp,k)+a(im,jm,k)+a(ip,j,kp)
     &   +a(ip,j,km)+a(im,j,kp)+a(im,j,km)+a(i,jp,kp)+a(i,jp,km)
     &   +a(i,jm,kp)+a(i,jm,km))/32
     &  +(a(ip,jp,kp)+a(ip,jp,km)+a(ip,jm,kp)+a(ip,jm,km)
     &   +a(im,jp,kp)+a(im,jp,km)+a(im,jm,kp)+a(im,jm,km))/64
            enddo
         enddo
c$dir force_parallel
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
c$dir force_vector
            do i=1,ng1
c               a(i,j,km)=v(i,j,ko)
c the above line results invokes a very subtle bug on the SGI compiler,
c resulting in wrong answers.
               a(i,j,km)=v(i,j,2-mod(k,2))
            enddo
         enddo
      enddo


c$dir force_parallel
c$omp do
      do j=1,ng2
c$dir force_vector
         do i=1,ng1
            a(i,j,ng3)=v(i,j,kv)
            a(i,j,1)=topa(i,j)
         enddo
      enddo
c$omp end parallel

      return
      end

c*********************************************************************
      subroutine spreadmax(a,scale,ld1)
c goal: to expand the compression limiter regions 
c it replaces a with the maximum of a and scale*max(neighbors)
      implicit none
      include 'relaxgl.fi'
      integer ld1
      real a(ld1,ng2,ng3),scale
c locals
      real v(ng1,ng2,0:1), a1,anmax
cmydist v(*,block,*)
      integer i,j,k,ip,im,jp,jm,kp,km,ko


cfpp$ skip
c*$* assert do(serial)      
c$omp parallel default(private) firstprivate(scale)
c$omp& shared(a,v)
      do k=1,ng3
c$dir force_parallel
c*$* assert do(concurrent)
c$omp do      
         do j=1,ng2
c$dir force_vector
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               a1=a(i,j,k)
               anmax=max(a(ip,j,k),a(im,j,k),a(i,jp,k),a(i,jm,k)
     &              ,a(i,j,kp),a(i,j,km))
               v(i,j,mod(k,2))=max(a1,scale*anmax)
            enddo
         enddo
         if (k .gt. 1) then
c*$* assert do(concurrent)
c$omp do      
            do j=1,ng2
c$dir force_vector
               do i=1,ng1
                  a(i,j,km)=v(i,j,1-mod(k,2))
               enddo
            enddo
         endif
      enddo
      k=1
      km=ng3
c*$* assert do(concurrent)
c$omp do      
      do j=1,ng2
c$dir force_vector
         do i=1,ng1
            a(i,j,km)=v(i,j,1-mod(k,2))
         enddo
      enddo
c$omp end parallel

      return
      end

























