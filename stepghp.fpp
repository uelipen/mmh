c -*- Fortran -*- file stepghp.f
c written May 11, 1995 by Ue-Li Pen, upen@astro.princeton.edu
c
c update particle positions and calculate gravitational field.
c
c rewritten June 1, 1995 to use more efficient chaining lists.
c
c modified June 20 to keep all particles on the common block,
c which gives more flexibility to have a shorter size for
c them.      
#ifdef NBODY
      
      subroutine pcalcrho(rho,def)
c on exit, rho will be modified to include the contribution
c from particles.      
c dts is the effective 2*dt*a^2 in units of \tau      
      implicit none
      include 'relaxgl.fi'
      real rho(ng1,ng2,ng3), def(ng1,ng2,ng3)
      
c locals
      logical firsttime
      data firsttime /.true./
      save firsttime
      integer i,j,k
      include 'globalpa.fi'

      include 'nbody.fi'
#ifdef COLD
      include 'cold.fi'
#endif
      real fpig, bomega
      parameter(fpig=2./3.)



      bomega=omegab/omega0
      if (firsttime) then
         firsttime = .false.
         write(*,*) 'Nbody solver: omega_b/omega0=',bomega
#ifdef P3M
         write(*,*) 'using CRUDE P3M'
#endif 
      endif


#ifdef COLD
c$omp parallel default(private) shared(rho,cold) firstprivate(bomega)
#else
c$omp parallel default(private) shared(rho) firstprivate(bomega)
#endif
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               rho(i,j,k)=rho(i,j,k)*bomega/(1-bomega)
#ifdef COLD
               if (cold(i,j,k))  rho(i,j,k)=0
#endif
            enddo
         enddo
      enddo
c$omp end parallel
c      call makechain
      call makechainparallel
      call calcrho(rho)
#ifdef COLD
c$omp parallel default(private) shared(rho,cold) firstprivate(bomega)
#else
c$omp parallel default(private) shared(rho) firstprivate(bomega)
#endif
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               rho(i,j,k)=fpig*rho(i,j,k)*(1-bomega)
#ifdef COLD
               if (cold(i,j,k)) rho(i,j,k)=rho(i,j,k)/(1-bomega)
#endif
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end



      subroutine limitx
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'

c compiler bugs in the SGI:  it doesn't like common block objects
c more than 2GB in size.      
#ifdef _SGI_SOURCE

      call limitx1(xv,npart)
      return
      end

      subroutine limitx1(xv,npart)
      implicit none
      include 'relaxgl.fi'
      integer npart
#ifndef XVREAL8
      real*4 xv(6,npart)
#else
      real xv(6,npart)
#endif      
#endif      
c local
      integer i
      
c$omp parallel do default(shared) private(i)
      do i=1,npart
         if (xv(1,i) .lt. 1) xv(1,i)=xv(1,i)+ng1
c  it is possible for a particle to
c go several times around the grid in one step, which one needs to test
c for separately.
         if (xv(2,i) .lt. 1) xv(2,i)=xv(2,i)+ng2
         if (xv(3,i) .lt. 1) xv(3,i)=xv(3,i)+ng3
         if (xv(1,i) .ge. ng1+1) xv(1,i)=xv(1,i)-ng1
         if (xv(2,i) .ge. ng2+1) xv(2,i)=xv(2,i)-ng2
         if (xv(3,i) .ge. ng3+1) xv(3,i)=xv(3,i)-ng3
      enddo
      return
      end


      
      subroutine matchhydro(tmp,phi,u,def)
c match the gas field u for density and momentum      
      implicit none
      include 'relaxgl.fi'
      real u(5,ng1,ng2,ng3),tmp(ng1,ng2,ng3),phi(ng1,ng2,ng3)
     &     ,def(ng1,ng2,ng3)
      include 'nbody.fi'

#ifdef ZELDOVICH
c displaces the particles using the Zeldovich approximation
c
      integer i,j,k,ip,jp,kp,ix,n1,n2,n3
      real rms,rmean,rmom,xoff, drho

#if defined(_SGI_SOURCE) || defined(_ALPHA)
      integer ii
      integer ifoldy(ng2),ifoldz(ng1)
      real pi,xk,yk,zk,akernel
      parameter(pi=3.14159265358979323846264338328)
#endif


      n3=ng3*nclinear
      n2=ng2*nclinear
      if (ng2 .eq. 1 ) n2=1
      n1=ng1*nclinear
      if (ng1 .eq. 1) n1=1
      rms=0
c$omp parallel default(private) reduction(+:rms) shared(u)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               rms=rms+(u(1,i,j,k)-1)**2
            enddo
         enddo
      enddo
c$omp end parallel
      rms=sqrt(rms/(ng1*ng2*ng3))
      write(*,*)'matchrho: initializing particles, rms=',rms
       
#if defined(_SGI_SOURCE) || defined(_ALPHA)
      write(*,*) 'matchhydro:n123=',n1,n2,n3
      write(*,*) 'matchhydro: npartcell,nclinear=',npartcell,nclinear
      
c$omp parallel default(private) shared(u,rhostartsmall)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
c transpose so the extra dimension is on NG3, makes life easier
c for 1-D simulations               
               drho=u(1,i,j,k)-1
               rhostartsmall(k,j,i)=1/(1+drho)
            enddo
         enddo
      enddo
c$omp end parallel

      call sfft3(rhostartsmall,ng1,1)
      
      do i=1,(ng2+1)/2
         ifoldy(i)=i
      enddo
      do i=(ng2+1)/2+1,ng2
         ifoldy(i)=n2-ng2+i
      enddo
      do i=1,(ng1+1)/2
         ifoldz(i)=i
      enddo
      do i=(ng1+1)/2+1,ng1
         ifoldz(i)=n1-ng1+i
      enddo
c$omp parallel default(private) shared(rhostartup,ifoldy,ifoldz,rhostartsmall)
c$omp& firstprivate(n1,n2,n3) 
      do k=1,n1
c$omp do
         do j=1,n2
            do i=1,n3+2
               rhostartup(i,j,k)=0
            enddo
         enddo
      enddo
c invert the laplacian kernel      
c*$* assert do (concurrent)
      do k=1,ng1
c$omp do
         do j=1,ng2
            do i=1,ng3+2
               ii=(i-1)/2
               xk=2*pi*ii/n3
               yk=2*pi*(ifoldy(j)-1.)/n2
               zk=2*pi*(ifoldz(k)-1.)/n1
               akernel=2*cos(xk)+2*cos(yk)+2*cos(zk)-6
               if (abs(akernel).lt. 1e-10) akernel=1
               rhostartup(i,ifoldy(j),ifoldz(k))=rhostartsmall(i,j,k)
     &             /akernel /nclinear
c we divide by nclinear because the gradient operation
c will need the unit conversions               
            enddo
         enddo
      enddo
c$omp end parallel


      call sfft3(rhostartup,n1,-1)
c      write(*,*) (rhostartup(i,1,1),i=1,n3)

c$omp parallel private(i,j) shared(xv)
      do i=1,npart
c$omp do
         do j=1,6
            xv(j,i)=0
         enddo
      enddo
c$omp end parallel
c      xv=0
#else
      if (npartcell .ne. 1) then
         write(*,*) 'matchrho: need FFT to match uneven particles'
         pause
      endif
c
c$omp parallel default(private) shared(tmp,u,phi)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               tmp(i,j,k)=-u(1,i,j,k)
               phi(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel
      
      call multigrid(phi,tmp,def,1.,ng1,ng2,ng3,1,1)
      do i=1,10
         call multigrid(phi,tmp,def,1.,ng1,ng2,ng3,1,1+32)
      enddo
c$omp parallel default(private) shared(rhostartup,phi)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
c transpose so the extra dimension is on NG3, makes life easier
c for 1-D simulations               
               rhostartup(k,j,i)=phi(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel

#endif
      xoff=1.+0.5/nclinear
c*$* assert do (concurrent)
c$omp parallel default(private) firstprivate(n1,n2,n3,xoff)
c$omp& shared(pmass,xv,rhostartup,u)
      do k=1,n1
         kp=mod(k,n1)+1
c$omp do
         do j=1,n2
            jp=mod(j,n2)+1
            do i=1,n3
               ip=mod(i,n3)+1
               ix=k+(j-1)*n1+(i-1)*n1*n2
               pmass(ix)=1.
               xv(3,ix)=(i-1.)/nclinear+(rhostartup(ip,j,k)
     &              -rhostartup(i,j,k)
     &              + rhostartup(ip,jp,k)-rhostartup(i,jp,k)
     &              + rhostartup(ip,j,kp)-rhostartup(i,j,kp)
     &              + rhostartup(ip,jp,kp)-rhostartup(i,jp,kp))/4
     &              + xoff
               xv(2,ix)=(j-1.)/nclinear+(rhostartup(i,jp,k)
     &              -rhostartup(i,j,k)
     &              + rhostartup(i,jp,kp)-rhostartup(i,j,kp)
     &              + rhostartup(ip,jp,k)-rhostartup(ip,j,k)
     &              + rhostartup(ip,jp,kp)-rhostartup(ip,j,kp))/4
     &              + xoff
               xv(1,ix)=(k-1.)/nclinear+(rhostartup(i,j,kp)
     &              -rhostartup(i,j,k)
     &              + rhostartup(i,jp,kp)-rhostartup(i,jp,k)
     &              + rhostartup(ip,j,kp)-rhostartup(ip,j,k)
     &              + rhostartup(ip,jp,kp)-rhostartup(ip,jp,k))/4
     &              + xoff
c
#if ! (defined(_SGI_SOURCE) || defined(_ALPHA))
                rmean=(u(1,i,j,k)+u(1,i,j,kp)+u(1,i,jp,k)+u(1,i,jp,kp)
     &   +u(1,ip,j,k)+u(1,ip,j,kp)+u(1,ip,jp,k)+u(1,ip,jp,kp))/8
               rmom=(u(2,i,j,k)+u(2,i,j,kp)+u(2,i,jp,k)+u(2,i,jp,kp)
     &     +u(2,ip,j,k)+u(2,ip,j,kp)+u(2,ip,jp,k)+u(2,ip,jp,kp))/8
               xv(4,ix)=rmom/rmean
               rmom=(u(3,i,j,k)+u(3,i,j,kp)+u(3,i,jp,k)+u(3,i,jp,kp)
     &     +u(3,ip,j,k)+u(3,ip,j,kp)+u(3,ip,jp,k)+u(3,ip,jp,kp))/8
               xv(5,ix)=rmom/rmean
               rmom=(u(4,i,j,k)+u(4,i,j,kp)+u(4,i,jp,k)+u(4,i,jp,kp)
     &     +u(4,ip,j,k)+u(4,ip,j,kp)+u(4,ip,jp,k)+u(4,ip,jp,kp))/8
               xv(6,ix)=rmom/rmean
#endif               
            enddo
         enddo
      enddo
c$omp end parallel

c      call makechain

c (end of ZELDOVICH)
#else
c particles have no displacement, the fluctuations from from mass fluctuations
c alone.
      real xoff,wx,wy,wz,w000,w001,w010,w011,w100,w101,w110,w111,p1
     &     ,rho000,rho001,rho010,rho011,rho100,rho101,rho110,rho111,rms
     &     ,drho
      integer i,j,k,n1,n2,n3,i1,j1,k1,ixp,iyp,izp,ix, idist, iseed(4)


      rms=0
c$omp parallel default(private) shared(u) reduction(+:rms)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               rms=rms+(u(1,i,j,k)-1)**2
            enddo
         enddo
      enddo
c$omp end parallel
      rms=sqrt(rms/(ng1*ng2*ng3))
      write(*,*)'matchrho(no zeldovich): initializing particles, rms='
     &    ,rms
      rms=0
c$omp parallel default(private) shared(def) reduction(+:rms)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               rms=rms+def(i,j,k)**2
            enddo
         enddo
      enddo
c$omp end parallel
      rms=sqrt(rms/(ng1*ng2*ng3))
      if (rms .gt. 0) then
         write(*,*) 'initial grid already deformed, centering particles'
         xoff=1
      else
         xoff=1.+0.5/nclinear
#ifdef P3M
         write(*,*) 'P3M: using centered particles'
         xoff=1.
#endif
      endif
      n3=ng3*nclinear
      n2=ng2*nclinear
      if (ng2 .eq. 1 ) n2=1
      n1=ng1*nclinear
      if (ng1 .eq. 1) n1=1
c*$* assert do (concurrent)
c$omp parallel default(private) firstprivate(n1,n2,n3,xoff)
c$omp& shared(xv,u,pmass)
      do k=1,n1
c$omp do
         do j=1,n2
!cdir nodep
            do i=1,n3
               ix=k+(j-1)*n1+(i-1)*n1*n2
               xv(3,ix)=(i-1.)/nclinear
     &              + xoff
               xv(2,ix)=(j-1.)/nclinear
     &              + xoff
               xv(1,ix)=(k-1.)/nclinear
     &              + xoff
               i1=xv(1,ix)
               j1=xv(2,ix)
               k1=xv(3,ix)
               wx=xv(1,ix)-i1
               wy=xv(2,ix)-j1
               wz=xv(3,ix)-k1
               w111=wx*wy*wz
               w110=wx*wy*(1-wz)
               w101=wx*(1-wy)*wz
               w100=wx*(1-wy)*(1-wz)
               w011=(1-wx)*wy*wz
               w010=(1-wx)*wy*(1-wz)
               w001=(1-wx)*(1-wy)*wz
               w000=(1-wx)*(1-wy)*(1-wz)
               ixp=mod(i1,ng1)+1
               iyp=mod(j1,ng2)+1
               izp=mod(k1,ng3)+1
               rho000=u(1,i1,j1,k1)
               rho001=u(1,i1,j1,izp)
               rho010=u(1,i1,iyp,k1)
               rho011=u(1,i1,iyp,izp)
               rho100=u(1,ixp,j1,k1)
               rho101=u(1,ixp,j1,izp)
               rho110=u(1,ixp,iyp,k1)
               rho111=u(1,ixp,iyp,izp)
               p1=rho000*w000+rho001*w001+rho010*w010
     &           +rho011*w011+rho100*w100+rho101*w101
     &           +rho110*w110+rho111*w111
               pmass(ix)=p1
               xv(4,ix)=0
               xv(5,ix)=0
               xv(6,ix)=0
            enddo
         enddo
      enddo
c$omp end parallel
#endif

c if the initial grid is already deformed, we need to take special
c care that the particles match the density EXACTLY.  Otherwise we 
c generate  spurious growing modes, which can be deadly.

c$omp parallel default(private)shared(tmp)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               tmp(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel
      call limitx
      call makechainparallel
      if (rms .gt. 0 .and. nclinear .gt. 1) then
c we need to adjust the particle densities to match perfectly.
         call calcrho(tmp)
         rms=0
c$omp parallel default(private) reduction(+:rms) firstprivate(n1,n2)
c$omp& shared(u,tmp,pmass)
         do k=1,ng3
c$omp do
            do j=1,ng2
!cdir nodep
               do i=1,ng1
                  drho=(u(1,i,j,k)-tmp(i,j,k))
                  rms=rms+drho**2
                  ix=(i-1)*nclinear+1+(j-1)*nclinear*n1
     &                 +(k-1)*n1*n2*nclinear
               pmass(ix)=pmass(ix)+nclinear**3*(u(1,i,j,k)-tmp(i,j,k))
                  tmp(i,j,k)=0
               enddo
            enddo
         enddo
c$omp end parallel
         rms=sqrt(rms/(ng1*ng2*ng3))
         write(*,*) 'curvilinear error=',rms
      endif

      if (nmassless .gt. 0 ) then
         if (nmassless .ne. ng1*ng2*ng3) then
            write(*,*) 'matchhydro:'
            write(*,*) 'dont know what to do with nmassless=',nmassless
            write(*,*) 'try setting it equal to ng1*ng2*ng3'
            stop
         endif
         idist=1
c uniform distribution [0,1)
         iseed(1)=43
         iseed(2)=12
         iseed(3)=95
         iseed(4)=19
         call slarnv(idist, iseed, 6*nmassless, xv(1,npmassive+1))
c$omp parallel default(private) shared(xv,pmass)
         do k=1,ng3
c$omp do
            do j=1,ng2
               do i=1,ng1
                  ix=k+(j-1)*ng1+(i-1)*ng1*ng2+npmassive
                  xv(1,ix)=xv(1,ix)+i
                  xv(2,ix)=xv(2,ix)+j
                  xv(3,ix)=xv(3,ix)+k
                  pmass(ix)=0
                  xv(4,ix)=0
                  xv(5,ix)=0
                  xv(6,ix)=0
               enddo
            enddo
         enddo
c$omp end parallel
      endif
      call limitx
      call makechainparallel
      call calcrho(tmp)
      rms=0
c$omp parallel default(private) shared(u,tmp) reduction(+:rms)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               drho=(u(1,i,j,k)-tmp(i,j,k))
               rms=rms+drho**2
            enddo
         enddo
      enddo
c$omp end parallel
      rms=sqrt(rms/(ng1*ng2*ng3))
      write(*,*)'bottom of matchrho: rms particle-gas=',rms
      return
      end

c end of MATCHHYDRO

#ifndef _ALPHA
      
      subroutine makechain
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      
      
c locals
      integer i,j,k,ix,iy,iz,missed,nhead
c$omp parallel default(private) shared(hoc)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               hoc(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel
c*$* assert do(concurrent)      
      do i=1,npmassive
         ix=xv(1,i)
         iy=xv(2,i)
         iz=xv(3,i)
         ll(i)=-hoc(ix,iy,iz)
         hoc(ix,iy,iz)=i
      enddo

c for parallel processing, we need to check if the procedure
c resulted in data conflicts.
c
      nhead=0
c I really dont know why I did that.
c*$* assert do(concurrent)      
      do i=1,npmassive
         j=abs(ll(i))
         if (j .gt. 0) then
	    ll(j)=-ll(j)
         else
             nhead=nhead+1
         endif
      enddo

c*$* assert do(concurrent)      
c zpj: compiler do not allow nhead to be reduction. may not parallizable
cc$omp parallel default(private) shared(hoc,ll) reduction(-:nhead)
      do k=1,ng3
cc$omp do
         do j=1,ng2
            do i=1,ng1
               ix=hoc(i,j,k)
               if (ix .gt. 0) then
                  ll(ix)=-ll(ix)
                  nhead=nhead-1
               endif
            enddo
         enddo
      enddo
cc$omp end parallel
      missed=nhead
c$doacross local(i) reduction(missed)
      do i=1,npmassive
         if( ll(i) .lt. 0 ) missed=missed+1
      enddo
      if ( missed .gt. 0 ) then
         write(*,*)' makechain: redoing scalar. nhead,missed='
     &	      ,nhead,missed
         call makechainscalar
      endif
      return
      end


      subroutine makechainscalar
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      
c locals
      integer i,j,k,ix,iy,iz

c$omp parallel default(private) shared(hoc)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               hoc(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel
c it works as follows:
c if hoc(i,j,k) = 0, that means there are no particles in that cell
c else, hoc(i,j,k) is the first particle in that cell.
c the next particle is given by ll(hoc(i,j,k)), etc.
         
      do i=1,npmassive
         ix=xv(1,i)
         iy=xv(2,i)
         iz=xv(3,i)
         ll(i)=hoc(ix,iy,iz)
c each hit puts the new particle at the head-of-chain.
c if ll(i)=0, there are no more particles
         hoc(ix,iy,iz)=i
      enddo
      return
      end

#endif

      subroutine makechainparallel
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      
c locals
      integer nparallel
#ifdef _SX5
      parameter(nparallel=1024)
#else
#ifdef _SGI_SOURCE
      parameter(nparallel=256)
#else
      parameter(nparallel=64)
#endif
#endif
      integer i,j,k,ix,iy,iz, hocpar(3,ng3,nparallel), nchunk, jl, jh
     &     ,inext,j1,ii
c hoc(1,,) stores the Head-Of-Chain, hoc(2,,) the length of the chain,
c and hoc(3,,) stores the Tail-Of-Chain.      
      logical firsttime
      data firsttime /.true./


#if ! (defined(_SGI_SOURCE) || defined(_ALPHA) || defined(_SX5) )
      call makechainscalar
      return
      
#endif
      
      if (firsttime) then
         firsttime=.false.
         write(*,*) ' using parallel makechain, nparallel=',nparallel
      endif

c the new twist is that we first construct a 1-D chain using only the z-axis
c$omp parallel default(private) shared(hocpar)
      do i=1,nparallel
c$omp do
         do j=1,ng3
            hocpar(1,j,i)=0
            hocpar(2,j,i)=0
            hocpar(3,j,i)=0
         enddo
      enddo
c$omp end parallel
      nchunk=npmassive/nparallel
      if (nchunk*nparallel .ne. npmassive) then
         write(*,*) 'parallel makechain: need npmassive divisible by ',
     &        'nparallel',npmassive,nparallel
         stop
      endif
c$omp parallel default(private) firstprivate(nchunk)
c$omp& shared(xv,hocpar,ll)
#ifdef _SX5
c on a vector machine, we want the vector loop to be inside      
      do j1=0,nchunk-1
!cdir pardo for
         do i=1,nparallel
            j=j1*nparallel+i
#else  
c zpj: I made the following lines comment
cc$doacross  shared (HOCPAR,XV,LL) local (J,IZ,II,jl,jh)
cc$&  MP_SCHEDTYPE=GSS
c zpj: end
c*KAP*parallel region  shared (HOCPAR,XV,LL) local (J,IZ,II,jl,jh)
c*KAP*parallel do
c zpj: this is one of the two do loops which does not parallize j loop
c$omp do
      do ii=1,nparallel
         jl=nchunk*(ii-1)+1
         jh=nchunk*ii
         do j=jl,jh
#endif
            iz=xv(3,j)
c the first hit is the tail-of-chain
            if (hocpar(1,iz,ii) .eq. 0) hocpar(3,iz,ii)=j
            ll(j)=hocpar(1,iz,ii)
            hocpar(1,iz,ii)=j
            hocpar(2,iz,ii)=hocpar(2,iz,ii)+1
         enddo
      enddo
c$omp end parallel
c$omp parallel default(private) shared(hocpar,ll)
#ifdef _SX5
c$omp do
      do j=2,nparallel
!cdir nodep
         do i=1,ng3
#else            
c*KAP*end parallel region
c*$* assert do(concurrent)      
c$omp do
      do i=1,ng3
         do j=2,nparallel
#endif
c now we merge the whole structure into hoc(1,:,1)
            if (hocpar(2,i,j) .gt. 0) then
               ll(hocpar(3,i,j))=hocpar(1,i,1)
               hocpar(1,i,1)=hocpar(1,i,j)
               hocpar(2,i,1)=hocpar(2,i,1)+hocpar(2,i,j)
            endif
         enddo
      enddo
c$omp end parallel
#ifdef _SX5
      icountz=hocpar(2,:,1)
      koffset(1)=0
      do i=2,ng3
         koffset(i)=koffset(i-1)+icountz(i-1)
      enddo
      i=maxval(icountz)
c$omp parallel default(private) firstprivate(i,icountz,koffsset)
c$omp& shared(ll2,hocpar,ll) 
c$omp do
      do j=1,i
         where (j .le. icountz)
            ll2(koffset+j)=hocpar(1,:,1)
            hocpar(1,:,1)=ll(hocpar(1,:,1))
         endwhere
      enddo
c$omp end parallel
c this isnt optimal if nparallel>>ng3
c in that case, one should make sure that particles are not arranged
c in ascending z, and do the loop parallel over hocpar(,,:)      

      
c on a vector machine, we build the linked list in calcrho      
      return
#endif           
      
c$omp parallel default(private) shared(hoc,hocpar,ll,xv)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               hoc(i,j,k)=0
            enddo
         enddo
      enddo

cc*$* assert do(concurrent)  MP_SCHEDTYPE=GSS
c zpj: I made the following comments
cc$ doacross local(i,iz,inext,ix,iy,j,k) shared(hocpar,ll,xv)
cc$&  MP_SCHEDTYPE=GSS
c zpj: end
c*KAP*parallel region local(i,iz,inext,ix,iy,j,k) 
c*KAP*& shared(hoc,hocpar,ll,xv)
c*KAP*parallel do
c$omp do
      do j=1,ng3
c now rebuild the whole structure, with the advantage that everything
c is concurrent: we can do each z-layer in parallel.
c Dont be confused that j actual refers to the z-index
         i=hocpar(1,j,1)
         iz=j
c now we proceed fearlessly, knowing that the whole chain lies on our
c z-layer.
         do k=1,hocpar(2,j,1)
            inext=ll(i)
            ix=xv(1,i)
            iy=xv(2,i)
            ll(i)=hoc(ix,iy,iz)
            hoc(ix,iy,iz)=i
            i=inext
         enddo
      enddo
c*KAP*end parallel region
c$omp end parallel
      return
      end




      
      subroutine calcrho1(rho)
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      real rho(ng1,ng2,ng3)
c 27 M 23 A
c locals      
      real wx,wy,wz,p1,w000,w001,w010,w011,w100,w101,w110,w111
     &     , amass,rho000,rho001,rho010,rho011,rho100,rho101,rho110
     &     , rho111
      integer i,j,k,ix,iy,iz,ixp,iyp,izp,ipart
      logical fexception

      amass=0
      do k=1,ng3
         amass=amass-rho(1,1,k)
      enddo
      i=1
      j=1
      ix=1
      iy=1
      do k=1,ng3
         iz=k
         izp=mod(iz,ng3)+1
         rho000=rho(ix,iy,iz)
         rho001=rho(ix,iy,izp)
         
         ipart=hoc(i,j,k)
         do while (ipart .ne. 0)
            wz=xv(3,ipart)-iz
            p1=pmass(ipart)
            w001=wz
            w000=(1-wz)

            rho000=rho000+w000*p1
            rho001=rho001+w001*p1
            ipart=ll(ipart)
         enddo
         rho(ix,iy,iz)=rho000
         rho(ix,iy,izp)=rho001
      enddo
      do k=1,ng3
         amass=amass+rho(i,j,k)
      enddo
      p1=1./npartcell
      if (abs(amass-npmassive*p1) .gt. 0.01) then
         write(*,*)'calcrho: missed mass=',npart*p1-amass
      endif
      return
      end

#ifdef _SX5
c on a vector machine we need a more fine grained splitup
c but we can sacrifice cache usage
c      
      subroutine calcrho(rho)
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      real rho(ng1,ng2,ng3)
c 27 M 23 A
c locals
      integer nveclen,nv
      parameter(nveclen=1024,nv=32)
      integer nhoc(nv,nv,nv,nveclen,3), nhoc1(nv**3,nveclen,3)
      equivalence (nhoc,nhoc1)
      real wx,wy,wz,p1,w000,w001,w010,w011,w100,w101,w110,w111
      integer i,j,k,ix,iy,iz,ixp,iyp,izp,ipart,ko,ii,jj,io,jo,ko
     &     , nchunk,j1,jjmax

      if (ng1*ng2 .eq. 1) then
         call calcrho1(rho)
         return
      endif

      nhoc=0
      nchunk=npmassive/nveclen
      if (nchunk*nveclen .ne. npmassive) then
         write(*,*) 'vector makechain: need npmassive divisible by ',
     &        'nveclen',npmassive,nveclen
         stop
      endif
      do i=0,nchunk-1
!cdir pardo for
         do j1=1,nveclen
            j=i*nveclen+j1
            ix=(xv(1,j)-1)*nv/ng1+1
            iy=(xv(2,j)-1)*nv/ng2+1
            iz=(xv(3,j)-1)*nv/ng3+1
c the first hit is the tail-of-chain
            if (nhoc(ix,iy,iz,j1,1) .eq. 0) nhoc(ix,iy,iz,j1,3)=j
            ll(j)=nhoc(ix,iy,iz,j1,1)
            nhoc(ix,iy,iz,j1,1)=j
            nhoc(ix,iy,iz,j1,2)=nhoc(ix,iy,iz,j1,2)+1
         enddo
      enddo
      do j=2,nveclen
c now we merge the whole structure into hoc(1,:,1)
         where (nhoc1(:,j,2) .gt. 0)
            ll(nhoc1(:,j,3))=nhoc1(:,1,1)
            nhoc1(:,1,1)=nhoc1(:,j,1)
            nhoc1(:,1,2)=nhoc1(:,1,2)+nhoc1(:,j,2)
         endwhere
      enddo

      do ko=0,1
         do jo=0,1
            do io=0,2
               jjmax=maxval(nhoc(io::2,jo::2,ko::2,1,2))
               do jj=0,jjmax
!cdir pardo for
!cdir nodep
                  do ii=0,nv**3/8-1
                  i=mod(ii,nv/2)*2+io
                  j=mod(ii*2/nv,nv/2)*2+jo
                  k=ii*4/nv**2*2+ko  ! note: int div is not associative
                  ipart=nhoc(i,j,k,1,1)
                  if (ipart > 0) then
                     nhoc(i,j,k,1,1)=ll(ipart)
                     ix=xv(1,ipart)
                     iy=xv(2,ipart)
                     iz=xv(3,ipart)
                     ixp=mod(ix,ng1)+1
                     iyp=mod(iy,ng2)+1
                     izp=mod(iz,ng3)+1
                     wx=xv(1,ipart)-ix
                     wy=xv(2,ipart)-iy
                     wz=xv(3,ipart)-iz
                     
                     p1=pmass(ipart)/npartcell
                     w111=wx*wy*wz
                     w110=wx*wy*(1-wz)
                     w101=wx*(1-wy)*wz
                     w100=wx*(1-wy)*(1-wz)
                     w011=(1-wx)*wy*wz
                     w010=(1-wx)*wy*(1-wz)
                     w001=(1-wx)*(1-wy)*wz
                     w000=(1-wx)*(1-wy)*(1-wz)
                     rho(ix,iy,iz)=rho(ix,iy,iz)+w000*p1
                     rho(ix,iy,izp)=rho(ix,iy,izp)+w001*p1
                     rho(ix,iyp,iz)=rho(ix,iyp,iz)+w010*p1
                     rho(ix,iyp,izp)=rho(ix,iyp,izp)+w011*p1
                     rho(ixp,iy,iz)=rho(ixp,iy,iz)+w100*p1
                     rho(ixp,iy,izp)=rho(ixp,iy,izp)+w101*p1
                     rho(ixp,iyp,iz)=rho(ixp,iyp,iz)+w110*p1
                     rho(ixp,iyp,izp)=rho(ixp,iyp,izp)+w111*p1
                  endif
               enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
      
#else      
      subroutine calcrho(rho)
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      real rho(ng1,ng2,ng3)
c 27 M 23 A
c locals      
      real wx,wy,wz,p1,w000,w001,w010,w011,w100,w101,w110,w111
     &     , rho000,rho001,rho010,rho011,rho100,rho101,rho110
     &     , rho111,asum,bsum
      integer i,j,k,ix,iy,iz,ixp,iyp,izp,ipart,ko
      logical fexception
      if (ng1*ng2 .eq. 1) then
         call calcrho1(rho)
         return
      endif
      fexception=.false.
c      call usum(asum,rho,ng1,ng2,ng3)
      
c$omp parallel default(private) shared(rho,hoc,xv,pmass,ll)
      do ko=1,2
c this should be perfectly safe to parallelize.  The compiler must
c put a synchronization barrier at the end of the k loop.   
c zpj: I made the following comments      
cc$doacross local(i,j,k,ix,iy,iz,ixp,iyp,izp,ipart,wx,wy,wz,
cc$&  w000,w001,w010,w011,w100,w101,w110,w111,p1,
cc$&  rho000,rho001,rho010,rho011,rho100,rho101,rho110,rho111)
cc$&  shared(xv,hoc,ll,rho,pmass, ng1,ng2,ng3,ko,npartcell)
c zpj: end
c*KAP*parallel region local(i,j,k,ix,iy,iz,ixp,iyp,izp,ipart,wx,wy,wz,
c*KAP*&  w000,w001,w010,w011,w100,w101,w110,w111,p1,
c*KAP*&  rho000,rho001,rho010,rho011,rho100,rho101,rho110,rho111)
c*KAP*&  shared(xv,hoc,ll,rho,pmass,ko)
c*KAP*parallel do
c zpj: the second do loop which does not parallize on j loop
c$omp do
      do k=ko,ng3,2
         do j=1,ng2
            do i=1,ng1
               ix=i
               iy=j
               iz=k
               ixp=mod(ix,ng1)+1
               iyp=mod(iy,ng2)+1
               izp=mod(iz,ng3)+1
               rho000=rho(ix,iy,iz)
               rho001=rho(ix,iy,izp)
               rho010=rho(ix,iyp,iz)
               rho011=rho(ix,iyp,izp)
               rho100=rho(ixp,iy,iz)
               rho101=rho(ixp,iy,izp)
               rho110=rho(ixp,iyp,iz)
               rho111=rho(ixp,iyp,izp)

               ipart=hoc(i,j,k)
#ifndef _ALPHA
               do while (ipart .ne. 0)
#else
 20               continue
                  if (ipart .eq. 0) goto 10
#endif
               wx=xv(1,ipart)-ix
               wy=xv(2,ipart)-iy
               wz=xv(3,ipart)-iz

               p1=pmass(ipart)/npartcell
c 3 S
c               
c               if (wx .lt. 0 .or. wx .gt. 1 .or. wy .lt. 0 .or. 
c     &             wy .gt. 1 .or. wz .lt. 0 .or. wz .gt. 1) then
c                  write(*,*) 'calcrho: error, wxyz=',wx,wy,wz
c                   fexception = fexception .or. .true.
c               endif
               w111=wx*wy*wz
               w110=wx*wy*(1-wz)
               w101=wx*(1-wy)*wz
               w100=wx*(1-wy)*(1-wz)
               w011=(1-wx)*wy*wz
               w010=(1-wx)*wy*(1-wz)
               w001=(1-wx)*(1-wy)*wz
               w000=(1-wx)*(1-wy)*(1-wz)
c 16 M 12 S

               rho000=rho000+w000*p1
               rho001=rho001+w001*p1
               rho010=rho010+w010*p1
               rho011=rho011+w011*p1
               rho100=rho100+w100*p1
               rho101=rho101+w101*p1
               rho110=rho110+w110*p1
               rho111=rho111+w111*p1
c 8 M 8 A               
               ipart=ll(ipart)
#ifndef _ALPHA
               enddo
#else
               goto 20
 10            continue
#endif
               rho(ix,iy,iz)=rho000
               rho(ix,iy,izp)=rho001
               rho(ix,iyp,iz)=rho010
               rho(ix,iyp,izp)=rho011
               rho(ixp,iy,iz)=rho100
               rho(ixp,iy,izp)=rho101
               rho(ixp,iyp,iz)=rho110
               rho(ixp,iyp,izp)=rho111
            enddo
         enddo
      enddo
c*KAP*end parallel region
c close ko loop      
      enddo
c$omp end parallel
c      call usum(bsum,rho,ng1,ng2,ng3)
c      write(*,*) 'calcrho:sum=',bsum-asum
      return
      end
#endif
      
      subroutine stepxv(phi,def,defn,tdt,odt,a)
c external entry routine
c dt is the current step size to step forward to.
c odt is the last time steps dt
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      real def(ng1,ng2,ng3)
     &     ,defn(ng1,ng2,ng3),phi(ng1,ng2,ng3)
     &     ,tdt, odt, a

#ifdef P3M
c      call p3m(def,tdt/2,a)
c I think we can do two steps at once by cheating a little,
c still needs to be checked if it is really 2nd order.
      call p3m(def,(odt+tdt)/2,a)
#endif
      call stepxv1(phi,def,defn,tdt,odt,a)
#ifdef P3M
c      call makechainparallel
c      call p3m(def,tdt/2,a)      
#endif
      
      end
      
      subroutine stepxv1(phi,def,defn,tdt,odt,a)
c external entry routine
c dt is the current step size to step forward to.
c odt is the last time steps dt
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      real def(ng1,ng2,ng3)
     &     ,defn(ng1,ng2,ng3),phi(ng1,ng2,ng3)
     &     ,tdt, odt, a
      
c locals
      integer k, iflip, ik0,ik1,kp
      real xvdot(2,12,ng1,ng2),dt


      if (tdt .lt. 1.e-10) then
         write(*,*) 'stepxv: dt=',tdt,' setting to 1E-10'
      endif
      dt=max(tdt,1.e-10)

      kp=1
      iflip=1
      call calcxvdot(xvdot,phi,def,defn,dt,a,kp,iflip)
               
      do k=1,ng3
         kp=mod(k,ng3)+1
         iflip=3-iflip
         call calcxvdot(xvdot,phi,def,defn,dt,a,kp,iflip)
         ik0=3-iflip
         ik1=iflip
         call pushparticle(xvdot,odt,dt,k,ik0,ik1)
      enddo

      call limitx
      return
      end

      subroutine pushparticle(xvdot,odt,tdt,k,ik0,ik1)
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
c we put the k index on xvdot first for cache purposes
      real xvdot(2,12,ng1,ng2),odt,tdt
      integer ik0,ik1,k
      
c locals
      real wx,wy,wz,w000,w001,w010,w011,w100,w101,w110
     &     ,w111,dv1,dv2,dv3,v1,v2,v3,dx1,dx2,dx3,vm11,vm12,vm13
     &     ,vm22,vm23,vm33
     &     ,xvdot1000,xvdot1001,xvdot1010,xvdot1011
     &     ,xvdot1100,xvdot1101,xvdot1110,xvdot1111
     &     ,xvdot2000,xvdot2001,xvdot2010,xvdot2011
     &     ,xvdot2100,xvdot2101,xvdot2110,xvdot2111
     &     ,xvdot3000,xvdot3001,xvdot3010,xvdot3011
     &     ,xvdot3100,xvdot3101,xvdot3110,xvdot3111
     &     ,xvdot4000,xvdot4001,xvdot4010,xvdot4011
     &     ,xvdot4100,xvdot4101,xvdot4110,xvdot4111
     &     ,xvdot5000,xvdot5001,xvdot5010,xvdot5011
     &     ,xvdot5100,xvdot5101,xvdot5110,xvdot5111
     &     ,xvdot6000,xvdot6001,xvdot6010,xvdot6011
     &     ,xvdot6100,xvdot6101,xvdot6110,xvdot6111
     &     ,xvdot7000,xvdot7001,xvdot7010,xvdot7011
     &     ,xvdot7100,xvdot7101,xvdot7110,xvdot7111
     &     ,xvdot8000,xvdot8001,xvdot8010,xvdot8011
     &     ,xvdot8100,xvdot8101,xvdot8110,xvdot8111
     &     ,xvdot9000,xvdot9001,xvdot9010,xvdot9011
     &     ,xvdot9100,xvdot9101,xvdot9110,xvdot9111
     &     ,xvdota000,xvdota001,xvdota010,xvdota011
     &     ,xvdota100,xvdota101,xvdota110,xvdota111
     &     ,xvdotb000,xvdotb001,xvdotb010,xvdotb011
     &     ,xvdotb100,xvdotb101,xvdotb110,xvdotb111
     &     ,xvdotc000,xvdotc001,xvdotc010,xvdotc011
     &     ,xvdotc100,xvdotc101,xvdotc110,xvdotc111
      integer i,j,ix,iy,iz,ixp,iyp,ipart
      logical fexception
     
c Grand total: 127 M, 117 A  -> 244 FLOP/particle

      fexception=.false.

      iz=k
#ifdef _SX5
!cdir pardo for
!cdir nodep      
      do i=1,icountz(k)
         ipart=ll2(koffset(k)+i)
         ix=xv(1,ipart)
         iy=xv(2,ipart)
#else
c zpj: I made the follwoing lines comments
cc$doacross local(i,j,ix,iy,ixp,iyp,ipart,wx,wy,wz,w000,w001,w010
cc$&,w011,w100,w101,w110,w111,dv1,dv2,dv3,v1,v2,v3,dx1,dx2,dx3,vm11,vm12
cc$&,vm13,vm22,vm23,vm33
cc$&     ,xvdot1000,xvdot1001,xvdot1010,xvdot1011
cc$&     ,xvdot1100,xvdot1101,xvdot1110,xvdot1111
cc$&     ,xvdot2000,xvdot2001,xvdot2010,xvdot2011
cc$&     ,xvdot2100,xvdot2101,xvdot2110,xvdot2111
cc$&     ,xvdot3000,xvdot3001,xvdot3010,xvdot3011
cc$&     ,xvdot3100,xvdot3101,xvdot3110,xvdot3111
cc$&     ,xvdot4000,xvdot4001,xvdot4010,xvdot4011
cc$&     ,xvdot4100,xvdot4101,xvdot4110,xvdot4111
cc$&     ,xvdot5000,xvdot5001,xvdot5010,xvdot5011
cc$&     ,xvdot5100,xvdot5101,xvdot5110,xvdot5111
cc$&     ,xvdot6000,xvdot6001,xvdot6010,xvdot6011
cc$&     ,xvdot6100,xvdot6101,xvdot6110,xvdot6111
cc$&     ,xvdot7000,xvdot7001,xvdot7010,xvdot7011
cc$&     ,xvdot7100,xvdot7101,xvdot7110,xvdot7111
cc$&     ,xvdot8000,xvdot8001,xvdot8010,xvdot8011
cc$&     ,xvdot8100,xvdot8101,xvdot8110,xvdot8111
cc$&     ,xvdot9000,xvdot9001,xvdot9010,xvdot9011
cc$&     ,xvdot9100,xvdot9101,xvdot9110,xvdot9111
cc$&     ,xvdota000,xvdota001,xvdota010,xvdota011
cc$&     ,xvdota100,xvdota101,xvdota110,xvdota111
cc$&     ,xvdotb000,xvdotb001,xvdotb010,xvdotb011
cc$&     ,xvdotb100,xvdotb101,xvdotb110,xvdotb111
cc$&     ,xvdotc000,xvdotc001,xvdotc010,xvdotc011
cc$&     ,xvdotc100,xvdotc101,xvdotc110,xvdotc111
cc$&) shared(ik0,ik1,iz)
cc$&  MP_SCHEDTYPE=GSS
c zpj:end
c*KAP*parallel region local(i,j,ix,iy,ixp,iyp,ipart,wx,wy,wz,w000
c*KAP*&,w001,w010,w011,w100,w101,w110,w111,dv1,dv2,dv3,v1,v2,v3
c*KAP*&,dx1,dx2,dx3,vm11,vm12,vm13,vm22,vm23,vm33
c*KAP*&     ,xvdot1000,xvdot1001,xvdot1010,xvdot1011
c*KAP*&     ,xvdot1100,xvdot1101,xvdot1110,xvdot1111
c*KAP*&     ,xvdot2000,xvdot2001,xvdot2010,xvdot2011
c*KAP*&     ,xvdot2100,xvdot2101,xvdot2110,xvdot2111
c*KAP*&     ,xvdot3000,xvdot3001,xvdot3010,xvdot3011
c*KAP*&     ,xvdot3100,xvdot3101,xvdot3110,xvdot3111
c*KAP*&     ,xvdot4000,xvdot4001,xvdot4010,xvdot4011
c*KAP*&     ,xvdot4100,xvdot4101,xvdot4110,xvdot4111
c*KAP*&     ,xvdot5000,xvdot5001,xvdot5010,xvdot5011
c*KAP*&     ,xvdot5100,xvdot5101,xvdot5110,xvdot5111
c*KAP*&     ,xvdot6000,xvdot6001,xvdot6010,xvdot6011
c*KAP*&     ,xvdot6100,xvdot6101,xvdot6110,xvdot6111
c*KAP*&     ,xvdot7000,xvdot7001,xvdot7010,xvdot7011
c*KAP*&     ,xvdot7100,xvdot7101,xvdot7110,xvdot7111
c*KAP*&     ,xvdot8000,xvdot8001,xvdot8010,xvdot8011
c*KAP*&     ,xvdot8100,xvdot8101,xvdot8110,xvdot8111
c*KAP*&     ,xvdot9000,xvdot9001,xvdot9010,xvdot9011
c*KAP*&     ,xvdot9100,xvdot9101,xvdot9110,xvdot9111
c*KAP*&     ,xvdota000,xvdota001,xvdota010,xvdota011
c*KAP*&     ,xvdota100,xvdota101,xvdota110,xvdota111
c*KAP*&     ,xvdotb000,xvdotb001,xvdotb010,xvdotb011
c*KAP*&     ,xvdotb100,xvdotb101,xvdotb110,xvdotb111
c*KAP*&     ,xvdotc000,xvdotc001,xvdotc010,xvdotc011
c*KAP*&     ,xvdotc100,xvdotc101,xvdotc110,xvdotc111
c*KAP*&) shared(ik0,ik1,iz,xv)
c*KAP*parallel do

c$omp parallel default(private) firstprivate(ik0,ik1,odt,tdt,k,iz)
c$omp& shared(hoc,xvdot,xv,ll)
c$omp do
      do j=1,ng2
         do i=1,ng1
            ix=i
            iy=j
            ipart=hoc(i,j,k)
#endif
            ixp=mod(ix,ng1)+1
            iyp=mod(iy,ng2)+1
c the compiler doesnt seem to understand while loops, so we lift
c invariant expressions by hand...            
            xvdot1000=xvdot(ik0,1,ix,iy)
            xvdot2000=xvdot(ik0,2,ix,iy)
            xvdot3000=xvdot(ik0,3,ix,iy)
            xvdot4000=xvdot(ik0,4,ix,iy)
            xvdot5000=xvdot(ik0,5,ix,iy)
            xvdot6000=xvdot(ik0,6,ix,iy)
            xvdot7000=xvdot(ik0,7,ix,iy)
            xvdot8000=xvdot(ik0,8,ix,iy)
            xvdot9000=xvdot(ik0,9,ix,iy)
            xvdota000=xvdot(ik0,10,ix,iy)
            xvdotb000=xvdot(ik0,11,ix,iy)
            xvdotc000=xvdot(ik0,12,ix,iy)
            
            xvdot1001=xvdot(ik1,1,ix,iy)
            xvdot2001=xvdot(ik1,2,ix,iy)
            xvdot3001=xvdot(ik1,3,ix,iy)
            xvdot4001=xvdot(ik1,4,ix,iy)
            xvdot5001=xvdot(ik1,5,ix,iy)
            xvdot6001=xvdot(ik1,6,ix,iy)
            xvdot7001=xvdot(ik1,7,ix,iy)
            xvdot8001=xvdot(ik1,8,ix,iy)
            xvdot9001=xvdot(ik1,9,ix,iy)
            xvdota001=xvdot(ik1,10,ix,iy)
            xvdotb001=xvdot(ik1,11,ix,iy)
            xvdotc001=xvdot(ik1,12,ix,iy)

            xvdot1010=xvdot(ik0,1,ix,iyp)
            xvdot2010=xvdot(ik0,2,ix,iyp)
            xvdot3010=xvdot(ik0,3,ix,iyp)
            xvdot4010=xvdot(ik0,4,ix,iyp)
            xvdot5010=xvdot(ik0,5,ix,iyp)
            xvdot6010=xvdot(ik0,6,ix,iyp)
            xvdot7010=xvdot(ik0,7,ix,iyp)
            xvdot8010=xvdot(ik0,8,ix,iyp)
            xvdot9010=xvdot(ik0,9,ix,iyp)
            xvdota010=xvdot(ik0,10,ix,iyp)
            xvdotb010=xvdot(ik0,11,ix,iyp)
            xvdotc010=xvdot(ik0,12,ix,iyp)
            
            xvdot1011=xvdot(ik1,1,ix,iyp)
            xvdot2011=xvdot(ik1,2,ix,iyp)
            xvdot3011=xvdot(ik1,3,ix,iyp)
            xvdot4011=xvdot(ik1,4,ix,iyp)
            xvdot5011=xvdot(ik1,5,ix,iyp)
            xvdot6011=xvdot(ik1,6,ix,iyp)
            xvdot7011=xvdot(ik1,7,ix,iyp)
            xvdot8011=xvdot(ik1,8,ix,iyp)
            xvdot9011=xvdot(ik1,9,ix,iyp)
            xvdota011=xvdot(ik1,10,ix,iyp)
            xvdotb011=xvdot(ik1,11,ix,iyp)
            xvdotc011=xvdot(ik1,12,ix,iyp)

            xvdot1100=xvdot(ik0,1,ixp,iy)
            xvdot2100=xvdot(ik0,2,ixp,iy)
            xvdot3100=xvdot(ik0,3,ixp,iy)
            xvdot4100=xvdot(ik0,4,ixp,iy)
            xvdot5100=xvdot(ik0,5,ixp,iy)
            xvdot6100=xvdot(ik0,6,ixp,iy)
            xvdot7100=xvdot(ik0,7,ixp,iy)
            xvdot8100=xvdot(ik0,8,ixp,iy)
            xvdot9100=xvdot(ik0,9,ixp,iy)
            xvdota100=xvdot(ik0,10,ixp,iy)
            xvdotb100=xvdot(ik0,11,ixp,iy)
            xvdotc100=xvdot(ik0,12,ixp,iy)
            
            xvdot1101=xvdot(ik1,1,ixp,iy)
            xvdot2101=xvdot(ik1,2,ixp,iy)
            xvdot3101=xvdot(ik1,3,ixp,iy)
            xvdot4101=xvdot(ik1,4,ixp,iy)
            xvdot5101=xvdot(ik1,5,ixp,iy)
            xvdot6101=xvdot(ik1,6,ixp,iy)
            xvdot7101=xvdot(ik1,7,ixp,iy)
            xvdot8101=xvdot(ik1,8,ixp,iy)
            xvdot9101=xvdot(ik1,9,ixp,iy)
            xvdota101=xvdot(ik1,10,ixp,iy)
            xvdotb101=xvdot(ik1,11,ixp,iy)
            xvdotc101=xvdot(ik1,12,ixp,iy)

            xvdot1110=xvdot(ik0,1,ixp,iyp)
            xvdot2110=xvdot(ik0,2,ixp,iyp)
            xvdot3110=xvdot(ik0,3,ixp,iyp)
            xvdot4110=xvdot(ik0,4,ixp,iyp)
            xvdot5110=xvdot(ik0,5,ixp,iyp)
            xvdot6110=xvdot(ik0,6,ixp,iyp)
            xvdot7110=xvdot(ik0,7,ixp,iyp)
            xvdot8110=xvdot(ik0,8,ixp,iyp)
            xvdot9110=xvdot(ik0,9,ixp,iyp)
            xvdota110=xvdot(ik0,10,ixp,iyp)
            xvdotb110=xvdot(ik0,11,ixp,iyp)
            xvdotc110=xvdot(ik0,12,ixp,iyp)
            
            xvdot1111=xvdot(ik1,1,ixp,iyp)
            xvdot2111=xvdot(ik1,2,ixp,iyp)
            xvdot3111=xvdot(ik1,3,ixp,iyp)
            xvdot4111=xvdot(ik1,4,ixp,iyp)
            xvdot5111=xvdot(ik1,5,ixp,iyp)
            xvdot6111=xvdot(ik1,6,ixp,iyp)
            xvdot7111=xvdot(ik1,7,ixp,iyp)
            xvdot8111=xvdot(ik1,8,ixp,iyp)
            xvdot9111=xvdot(ik1,9,ixp,iyp)
            xvdota111=xvdot(ik1,10,ixp,iyp)
            xvdotb111=xvdot(ik1,11,ixp,iyp)
            xvdotc111=xvdot(ik1,12,ixp,iyp)

#ifndef _SX5             
cSGI            nii=0
#ifndef _ALPHA
            do while (ipart .ne. 0)
cSGI            nii=nii+1
cSGI            ipart=ll(ipart)
cSGI            enddo
cSGI            ipart=hoc(i,j,k)
c the SGI compiler does not know how to pipeline while loops
c but the unrolled loop is still too complicated
cSGI            do ii=1,nii
#else
 20         continue
            if (ipart .eq. 0) goto 10
#endif
#endif
c in our convention we require that the particle positions are in the
c interval [1,ng+1)
               
            wx=xv(1,ipart)-ix        ! the distance to that grid point
            wy=xv(2,ipart)-iy
            wz=xv(3,ipart)-iz
c 3 S
c            if (wx .lt. 0 .or. wx .gt. 1 .or. wy .lt. 0 .or. wy .gt. 1
c     &           .or. wz .lt. 0 .or. wz .gt. 1) then
c               fexception = .true.
c               write(*,*) 'pushparticle: error, wxyz=',wx,wy,wz
c            endif
            w111=wx*wy*wz
            w110=wx*wy*(1-wz)
            w101=wx*(1-wy)*wz
            w100=wx*(1-wy)*(1-wz)
            w011=(1-wx)*wy*wz
            w010=(1-wx)*wy*(1-wz)
            w001=(1-wx)*(1-wy)*wz
            w000=(1-wx)*(1-wy)*(1-wz)
c 16 M 12 S

            dv1=xvdot4000*w000+xvdot4001*w001+xvdot4010*w010
     &           +xvdot4011*w011+xvdot4100*w100+xvdot4101*w101
     &           +xvdot4110*w110+xvdot4111*w111
            dv2=xvdot5000*w000+xvdot5001*w001+xvdot5010*w010
     &           +xvdot5011*w011+xvdot5100*w100+xvdot5101*w101
     &           +xvdot5110*w110+xvdot5111*w111
            dv3=xvdot6000*w000+xvdot6001*w001+xvdot6010*w010
     &              +xvdot6011*w011+xvdot6100*w100+xvdot6101*w101
     &              +xvdot6110*w110+xvdot6111*w111
c 24 M, 21 A
            v1=xv(4,ipart)+dv1*odt/2
            v2=xv(5,ipart)+dv2*odt/2
            v3=xv(6,ipart)+dv3*odt/2
c 3 M 3 A
            dx1=xvdot1000*w000+xvdot1001*w001+xvdot1010*w010
     &           +xvdot1011*w011+xvdot1100*w100+xvdot1101*w101
     &           +xvdot1110*w110+xvdot1111*w111
            dx2=xvdot2000*w000+xvdot2001*w001+xvdot2010*w010
     &           +xvdot2011*w011+xvdot2100*w100+xvdot2101*w101
     &           +xvdot2110*w110+xvdot2111*w111
            dx3=xvdot3000*w000+xvdot3001*w001+xvdot3010*w010
     &           +xvdot3011*w011+xvdot3100*w100+xvdot3101*w101
     &           +xvdot3110*w110+xvdot3111*w111
            vm11=xvdot7000*w000+xvdot7001*w001+xvdot7010*w010
     &           +xvdot7011*w011+xvdot7100*w100+xvdot7101*w101
     &           +xvdot7110*w110+xvdot7111*w111
            vm12=xvdot8000*w000+xvdot8001*w001+xvdot8010*w010
     &           +xvdot8011*w011+xvdot8100*w100+xvdot8101*w101
     &           +xvdot8110*w110+xvdot8111*w111
            vm13=xvdot9000*w000+xvdot9001*w001+xvdot9010*w010
     &           +xvdot9011*w011+xvdot9100*w100+xvdot9101*w101
     &           +xvdot9110*w110+xvdot9111*w111
            vm22=xvdota000*w000+xvdota001*w001+xvdota010*w010
     &           +xvdota011*w011+xvdota100*w100+xvdota101*w101
     &           +xvdota110*w110+xvdota111*w111
            vm23=xvdotb000*w000+xvdotb001*w001+xvdotb010*w010
     &           +xvdotb011*w011+xvdotb100*w100+xvdotb101*w101
     &           +xvdotb110*w110+xvdotb111*w111
            vm33=xvdotc000*w000+xvdotc001*w001+xvdotc010*w010
     &           +xvdotc011*w011+xvdotc100*w100+xvdotc101*w101
     &           +xvdotc110*w110+xvdotc111*w111
c 72 M 63 A            

            dx1=dx1+vm11*v1+vm12*v2+vm13*v3
            dx2=dx2+vm12*v1+vm22*v2+vm23*v3
            dx3=dx3+vm13*v1+vm23*v2+vm33*v3
c 9 M 9 A
            xv(1,ipart)=xv(1,ipart)+dx1
            xv(2,ipart)=xv(2,ipart)+dx2
            xv(3,ipart)=xv(3,ipart)+dx3
            xv(4,ipart)=v1+dv1*tdt/2
            xv(5,ipart)=v2+dv2*tdt/2
            xv(6,ipart)=v3+dv3*tdt/2
c 3 M 6 A            
c             
#ifndef _SX5
            ipart=ll(ipart)
#ifndef _ALPHA
            enddo
#else
            goto 20
 10         continue
#endif
         enddo
#endif
      enddo
c*KAP*end parallel region
c$omp end parallel
      if (fexception) then
         write(*,*) 'pushparticle error, wxyz'
         pause
      endif
      return
      end

               
         
      
      subroutine calcxvdot(xvdot,phi,def,defn,dt,a
     &     ,k,iflip)
c defn(,,) is the deformation potential at the next time step      
      implicit none
      include 'relaxgl.fi'
      integer k,iflip
      real xvdot(2,12,ng1,ng2),phi(ng1,ng2,ng3),def(ng1,ng2,ng3)
     &     ,defn(ng1,ng2,ng3),dt,a
c dt is the time to step forward,
      
c locals
      integer i,j,ip,im,jp,jm,kp,km
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &     ,a23, a33, det, b11, b12, b13, b22, b23, b33, a1,a2,a3
     &     ,gradx,grady,gradz,dx1,dx2,dx3,bdot11,bdot12,bdot13,bdot21
     &     ,bdot22,bdot23,bdot31,bdot32,bdot33


c$omp parallel default(private) firstprivate(k,a,iflip,dt)
c$omp&  shared(def,phi,xvdot,defn)        
      kp=mod(k,ng3)+1 
      km=mod(k+ng3-2,ng3)+1
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
     &           -def(ip,jm,k)+def(im,jm,k))/4
            phiyz=(def(i,jp,kp)-def(i,jp,km)
     &           -def(i,jm,kp)+def(i,jm,km))/4
            phixz=(def(ip,j,kp)-def(im,j,kp)
     &           -def(ip,j,km)+def(im,j,km))/4
            a11=1+phixx
            a12=phixy
            a13=phixz
            a22=1+phiyy
            a23=phiyz
            a33=1+phizz
            det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &           -a12**2*a33
            b11=(a22*a33-a23**2)/det
            b12=(a13*a23-a12*a33)/det
            b13=(a12*a23-a13*a22)/det
            b22=(a11*a33-a13**2)/det
            b23=(a12*a13-a11*a23)/det
            b33=(a11*a22-a12**2)/det

            gradx=-(phi(ip,j,k)-phi(im,j,k))/2
            grady=-(phi(i,jp,k)-phi(i,jm,k))/2
            gradz=-(phi(i,j,kp)-phi(i,j,km))/2

            a1=(b11*gradx+b12*grady+b13*gradz)*a
            a2=(b12*gradx+b22*grady+b23*gradz)*a
            a3=(b13*gradx+b23*grady+b33*gradz)*a

            xvdot(iflip,4,i,j)=a1
            xvdot(iflip,5,i,j)=a2
            xvdot(iflip,6,i,j)=a3

            dx1=(b11*a1+b12*a2+b13*a3)*dt**2/2
            dx2=(b12*a1+b22*a2+b23*a3)*dt**2/2
            dx3=(b13*a1+b23*a2+b33*a3)*dt**2/2

c
c     we are trying to solve \dot{\xi^a}=e^a_i (v^i-\dphi^i)
c which to second order reads \xi^+ - \xi = (e^a_i + \dot{e}^a_i\dt/2) v^i \dt
c     + e^a_i \dot{v}^i\dt^2/2 - e^a_i \dphi^i \dt - (\dot{e}^a_i \dphi
c        + e^a_i \dot{\dphi})\dt^2/2
c            


            phixx=(defn(ip,j,k)-2*defn(i,j,k)+defn(im,j,k)
     &           -phixx)/dt
            phiyy=(defn(i,jp,k)-2*defn(i,j,k)+defn(i,jm,k)
     &           -phiyy)/dt
            phizz=(defn(i,j,kp)-2*defn(i,j,k)+defn(i,j,km)
     &           -phizz)/dt
            phixy=((defn(ip,jp,k)-defn(im,jp,k)
     &           -defn(ip,jm,k)+defn(im,jm,k))/4-phixy)/dt
            phiyz=((defn(i,jp,kp)-defn(i,jp,km)
     &           -defn(i,jm,kp)+defn(i,jm,km))/4-phiyz)/dt
            phixz=((defn(ip,j,kp)-defn(im,j,kp)
     &           -defn(ip,j,km)+defn(im,j,km))/4-phixz)/dt
            a11=phixx
            a12=phixy
            a13=phixz
            a22=phiyy
            a23=phiyz
            a33=phizz
c \dot{A^-1}  =  - A^-1 \dot{A} A^-1
c We note that ABA is symmetric if A and B are symmetric
            bdot11=-(b11*a11+b12*a12+b13*a13)
            bdot12=-(b11*a12+b12*a22+b13*a23)
            bdot13=-(b11*a13+b12*a23+b13*a33)
            bdot21=-(b12*a11+b22*a12+b23*a13)
            bdot22=-(b12*a12+b22*a22+b23*a23)
            bdot23=-(b12*a13+b22*a23+b23*a33)
            bdot31=-(b13*a11+b23*a12+b33*a13)
            bdot32=-(b13*a12+b23*a22+b33*a23)
            bdot33=-(b13*a13+b23*a23+b33*a33)

            bdot11=bdot11*b11+bdot12*b12+bdot13*b13
            bdot12=bdot11*b12+bdot12*b22+bdot13*b23
            bdot13=bdot11*b13+bdot12*b23+bdot13*b33
            bdot22=bdot21*b12+bdot22*b22+bdot23*b23
            bdot23=bdot21*b13+bdot22*b23+bdot23*b33
            bdot33=bdot31*b13+bdot32*b23+bdot33*b33

            
            xvdot(iflip,7,i,j)=b11*dt+bdot11*dt**2/2
            xvdot(iflip,8,i,j)=b12*dt+bdot12*dt**2/2
            xvdot(iflip,9,i,j)=b13*dt+bdot13*dt**2/2
            xvdot(iflip,10,i,j)=b22*dt+bdot22*dt**2/2
            xvdot(iflip,11,i,j)=b23*dt+bdot23*dt**2/2
            xvdot(iflip,12,i,j)=b33*dt+bdot33*dt**2/2

c store -\dot{\phi^i}
            gradx=(defn(im,j,k)-defn(ip,j,k)
     &           -def(im,j,k)+def(ip,j,k))/dt/2
            grady=(defn(i,jm,k)-defn(i,jp,k)
     &           -def(i,jm,k)+def(i,jp,k))/dt/2
            gradz=(defn(i,j,km)-defn(i,j,kp)
     &           -def(i,j,km)+def(i,j,kp))/dt/2

            a1=b11*gradx+b12*grady+b13*gradz
            a2=b12*gradx+b22*grady+b23*gradz
            a3=b13*gradx+b23*grady+b33*gradz
            
            a1=a1*dt+(bdot11*gradx+bdot12*grady+bdot13*gradz
     &           )*dt**2/2
            a2=a2*dt+(bdot12*gradx+bdot22*grady+bdot23*gradz
     &           )*dt**2/2
            a3=a3*dt+(bdot13*gradx+bdot23*grady+bdot33*gradz
     &           )*dt**2/2

            xvdot(iflip,1,i,j)=a1+dx1
            xvdot(iflip,2,i,j)=a2+dx2
            xvdot(iflip,3,i,j)=a3+dx3
         enddo
      enddo
c$omp end parallel
      return
      end

#ifdef P3M
      subroutine p3m(def,dt,a)
c
c
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      real def(ng1,ng2,ng3)
      real dt,a
! locals
      real ro, rmin, fc
      parameter (fc=2./3./4./3.141592653,ro=2.1,rmin=0.05)
      integer ix,iy,iz,ixmin,ixmax,iymin,iymax,izmin,izmax,k1,j1,i1,
     &     k2,j2,i2,ipart,ipart2,i,j,k,ko,kt
c for 512^3 particles, icount overflows at 32 bits.
      integer*8 icount
c havent done a careful estimate for how small a value of nro we can
c get away with.
      integer nro
      parameter (nro=ro/2.+2)
      real ebaj(6,ng1,ng2,nro)
      real x2,y2,z2,x,y,z,dx,dy,dz,dxc,dyc,dzc,rc,rcurv,s,fx,fy,fz,
     &     e11,e12,e13,e22,e23,e33, amax, dk, wz
      integer kmap(ng3)
      
      do k=1,nro-1
         kmap(k)=k
         call calcebaj(ebaj,def,kmap,k)
      enddo
      icount=0
      amax=0
      do k=1,ng3
         ko=mod(k+nro-2+ng3,ng3)+1
         kt=mod(ko-1,nro)+1
         kmap(ko)=kt
         call calcebaj(ebaj,def,kmap,kt)
c this OMP directive isnt really safe -- one needs to space
c processors at least 2*ro+1 iterations apart.         
cc$OMP parallel do default(private)
cc$OMP& shared(hoc,ebaj,xv,a,k,kmap) reduction(+:icount,max:amax)
cc$OMP& OMP_SCHEDTYPE=GSS         
c$omp parallel default(private) firstprivate(a,k,dt)
c$omp& shared(hoc,ebaj,xv,kmap,ll,pmass) reduction(+:icount)
c$omp&         reduction(max:amax) 
cc$omp& OMP_SCHEDTYPE=GSS    
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               ipart=hoc(i,j,k)
!               do while (ipart .ne. 0) 
 100              continue
                  if (ipart .eq. 0) goto 500 !exit
                  x=xv(1,ipart)
                  y=xv(2,ipart)
                  z=xv(3,ipart)
                  izmin=z   ! we only search downwards
                  wz=z-izmin
                  izmax=z+ro
                  do k1=izmin,izmax
                     dk=sqrt(ro**2-(max(0.,k1-k-wz))**2)
                     iymin=y-dk+ng2
                     iymax=y+dk+ng2
                     do j1=iymin,iymax
c could do a similar trick to dk using j1  
                        ixmin=x-dk+ng1 ! make sure we round down.
                        ixmax=x+dk+ng1
                        do i1=ixmin,ixmax
                           k2=mod(k1-1,ng3)+1
                           j2=mod(j1-1,ng2)+1
                           i2=mod(i1-1,ng1)+1
                           ipart2=hoc(i2,j2,k2)
!                          do
 200                       continue
                           if (ipart2 .eq. 0) goto 400 ! exit
                           z2=xv(3,ipart2)
                           dz=z-z2
                           if (dz .gt. ng3/2) dz=dz-ng3
                           if (dz .lt. -ng3/2) dz=dz+ng3
c I think the above test is never true
                           if (dz .gt. 0.) goto 300 ! like cycle
                           x2=xv(1,ipart2)
                           y2=xv(2,ipart2)
                           dx=x-x2
                           if (dx .gt. ng1/2) dx=dx-ng1
                           if (dx .lt. -ng1/2) dx=dx+ng1
                           dy=y-y2
                           if (dy .gt. ng2/2) dy=dy-ng2
                           if (dy .lt. -ng2/2) dy=dy+ng2
                           if (dz .eq. 0.) then
                              if (ipart .eq. ipart2) goto 300
                              if (dy .gt. 0.) goto 300
                              if (dy .eq. 0. .and. dx .gt. 0) goto 300
                           endif
c the above line may be wrong on the first time step
c                           if (ipart2 .eq. ipart) goto 300 
c this test should be redundant from the previous test
                           rcurv=sqrt(dx**2+dy**2+dz**2)
                           if (rcurv .gt. ro) goto 300
                           ix=x2+dx/2-0.5+ng1
                           ix=mod(ix,ng1)+1
                           iy=y2+dy/2-0.5+ng2
                           iy=mod(iy,ng2)+1
                           iz=z2+dz/2-0.5+ng3
                           iz=mod(iz,ng3)+1
c     for laziness and speed, we use NGP for metric
                           e11=ebaj(1,ix,iy,kmap(iz))
                           e12=ebaj(2,ix,iy,kmap(iz))
                           e13=ebaj(3,ix,iy,kmap(iz))
                           e22=ebaj(4,ix,iy,kmap(iz))
                           e23=ebaj(5,ix,iy,kmap(iz))
                           e33=ebaj(6,ix,iy,kmap(iz))
c suffix 'c' stands for cartesian comoving
                           dxc=dx*e11+dy*e12+dz*e13
                           dyc=dx*e12+dy*e22+dz*e23
                           dzc=dx*e13+dy*e23+dz*e33
                           rc=sqrt(dxc**2+dyc**2+dzc**2)
c hard force cutoff, corresponds to empty shell particle
                           if (rc .lt. rmin) goto 300 ! like cycle
                           s=a*fc*(1-rcurv**3/(2.8+rcurv**5)**(3./5.))
c for PP, we assume gas traces dark matter
                           fx=-s*dxc/rc**3
                           fy=-s*dyc/rc**3
                           fz=-s*dzc/rc**3
c in operator splitting, we can set \dot{x}=0 and still be 2nd order
                           xv(4,ipart)=xv(4,ipart)+fx*dt*pmass(ipart2)
                           xv(5,ipart)=xv(5,ipart)+fy*dt*pmass(ipart2)
                           xv(6,ipart)=xv(6,ipart)+fz*dt*pmass(ipart2)
                           xv(4,ipart2)=xv(4,ipart2)-fx*dt*pmass(ipart)
                           xv(5,ipart2)=xv(5,ipart2)-fy*dt*pmass(ipart)
                           xv(6,ipart2)=xv(6,ipart2)-fz*dt*pmass(ipart)
                           icount=icount+1
                           amax=max(amax,s/rc**3)
 300                       continue
                           ipart2=ll(ipart2)
c     enddo
                           goto 200
 400                       continue
                        enddo
                     enddo
                  enddo
                  ipart=ll(ipart)
c               enddo
                  goto 100
 500              continue
            enddo
         enddo
c$omp end parallel
      enddo
      write(*,*) 'cost,pcfl=',icount*1./npmassive,dt*sqrt(amax/2)
c the last number is something like an N-body CFL number, should be < 1
      return
      end


c**********************************************************************
      subroutine calcebaj(baj,def,indx,k)
      implicit none
      include 'relaxgl.fi'
      include 'gmetric.fi'
      integer k,indx(ng3)
      real baj(6,ng1,ng2,nstore),def(ng1,ng2,ng3)
cmydist baj(*,*,block,*),def(*,block,*)
c
c locals
      integer i,j,ip,jp,im,jm,kp,km
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &      , a23, a33, det, b11,b12,b13,b22,b23,b33



c$dir prefer_parallel_ext            
cc$doacross local(i,j)
c$omp parallel default(private) firstprivate(k)
c$omp& shared(def,baj,indx)
c$omp do
      do j=1,ng2
         do i=1,ng1
            ip=mod(i,ng1)+1
            im=mod(i+ng1-2,ng1)+1
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            kp=mod(k,ng3)+1
            km=mod(k+ng3-2,ng3)+1
            phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
            phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
            phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
            phixy=(def(ip,jp,k)-def(im,jp,k)
     &           -def(ip,jm,k)+def(im,jm,k))/4
            phiyz=(def(i,jp,kp)-def(i,jp,km)
     &           -def(i,jm,kp)+def(i,jm,km))/4
            phixz=(def(ip,j,kp)-def(im,j,kp)
     &           -def(ip,j,km)+def(im,j,km))/4
            a11=1+phixx
            a12=phixy
            a13=phixz
            a22=1+phiyy
            a23=phiyz
            a33=1+phizz
            baj(1,i,j,indx(k))=a11
            baj(2,i,j,indx(k))=a12
            baj(3,i,j,indx(k))=a13
            baj(4,i,j,indx(k))=a22
            baj(5,i,j,indx(k))=a23
            baj(6,i,j,indx(k))=a33
         enddo
      enddo
c$omp end parallel
      return
      end


      
#endif
      subroutine writepmoment(amap)
      implicit none
      include 'relaxgl.fi'
      include 'nbody.fi'
      real*4 amap(10,ng1,ng2,ng3)
      
c locals      
      real vx,vy,vz,vxx,vxy,vxz,vyy,vyz,vzz
      integer i,j,k,ipart,imass


c zpj: I made the following lines comments
cc$doacross local(i,j,k,ipart,vx,vy,vz,vxx,vxy,vxz,vyy,vyz,vzz,imass)
cc$& , shared(xv,hoc,ll, ng1,ng2,ng3)
c zpj: end
c$omp parallel default(private) shared(hoc,xv,amap,ll)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               vx=0
               vy=0
               vz=0
               vxx=0
               vxy=0
               vxz=0
               vyy=0
               vyz=0
               vzz=0
               imass=0
               ipart=hoc(i,j,k)
c     do while (ipart .ne. 0)
 20            continue
               if (ipart .eq. 0) goto 10

               imass=imass+1
               vx=vx+xv(4,ipart)
               vy=vy+xv(5,ipart)
               vz=vz+xv(6,ipart)
               vxx=vxx+xv(4,ipart)**2
               vyy=vyy+xv(5,ipart)**2
               vzz=vzz+xv(6,ipart)**2
               vxy=vxy+xv(4,ipart)*xv(5,ipart)
               vxz=vxz+xv(4,ipart)*xv(6,ipart)
               vyz=vyz+xv(5,ipart)*xv(6,ipart)
               
               ipart=ll(ipart)
c     enddo
               goto 20
 10            continue

               amap(1,i,j,k)=vx
               amap(2,i,j,k)=vy
               amap(3,i,j,k)=vz
               amap(4,i,j,k)=vxx
               amap(5,i,j,k)=vxy
               amap(6,i,j,k)=vxz
               amap(7,i,j,k)=vyy
               amap(8,i,j,k)=vyz
               amap(9,i,j,k)=vzz
               amap(10,i,j,k)=imass
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end
                  
#endif




























