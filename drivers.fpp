c -*- Fortran -*-
c*********************************************************************
c 11/21/97: reduced this file to a single driver routine.
c
c
c*********************************************************************
#define BBKS      
      call initialize
      call gendriver
      end



#ifdef BBKS
c***************************************************************************
c generate initial powerspectrum according to BBKS prescription.
c credits: based on routines written by Guohong Xu
c
c      
#ifdef _ALPHA
      subroutine gendriver
      implicit none
      integer nbuf
      parameter(nbuf=300)
      include 'relaxgl.fi'
      real flat(ng1*ng2*ng3+300+ng1*ng2*2/9,10)

c locals
      integer iwordsize,ioff
cdec$ alias get_page_size, "getpagesize"
      integer get_page_size,ipagesize
      external get_page_size

      ipagesize=get_page_size()
      write(*,*) 'pagesize=',ipagesize
      iwordsize=%loc(flat(2,1))-%loc(flat(1,1))
      write(*,*) 'wordsize=',iwordsize
      write(*,*) 'flat at ', mod(%loc(flat(1,1)),ipagesize),
     &     mod(%loc(flat(2,1)),ipagesize)
      ioff=(ipagesize-mod(%loc(flat(1,1)),ipagesize))/iwordsize+1
      write(*,*) 'ioff=',ioff
      write(*,*) 'flat ioff ', mod(%loc(flat(ioff,1)),ipagesize),
     &     mod(%loc(flat(ioff+1,1)),ipagesize)
      call gendriver1(flat(ioff,1),flat(ioff,6),flat(ioff,7),
     &     flat(ioff,8),flat(ioff,9),flat(ioff,10))
      end
      subroutine gendriver1(u,phi,tmp,def,defp,tmp2)
      
#else
      subroutine gendriver
#endif
      implicit none
c everything is local
      include 'relaxgl.fi'
      include 'globalpa.fi'
      integer nfluidcmp,i,n,ip,im,j,k,loopni, nzout, ndimzout
      parameter(nfluidcmp=5,ndimzout=5)
      real pi,xk,yk,zk,sum,akernel
      parameter(pi=3.14159265358979311599796346854)
c modified 3/30/00 to use dimensions consistently, at the cost of some
c extra memory.  With the current design of limiter.fpp, one can not
c equivalance tmp and tmp2.      
c previously, tmp(2+,,) was also used as tmp()      
      real u(5,ng1,ng2,ng3), phi(ng1,ng2,ng3)
      real tmp(ng1,ng2,ng3),def(ng1,ng2,ng3),defp(ng1,ng2,ng3)
      real tmp2(ng1+2,ng2,ng3)
      real ai,af,cfl,ampl,anorm,pspect, afreq, rmin, rmax,qk
     &      ,pk,x, cosm_a, grow0,grow1,rms,zout(ndimzout), a 
cmydist u(*,*,block,*),phi(*,block,*),def(*,block,*)
cmydist defp(*,block,*),tmp(*,block,*),tmp2(*,block,*)
c these directives were for the SGI origin 2000 compiler.
c the tmp array is generally used as a normal array:
      character*80 scmd
c      data zout /32.,16.,1.,0.5,0./
      data zout /4.,2.,1.,0.5,0./
c      data zout /0.,-.1,-.1,-.1/
      integer  ii, nwave, iopt, time, is3, io, jo, ko, i1, i2, i3
#ifdef _ALPHA
      integer magic1
cdec$ alias get_page_size, "getpagesize"
      integer get_page_size,ipagesize
      external get_page_size

ccdec$ alias migrate_next_touch,"_OtsMigrateNextTouch"
c      external migrate_next_touch

cdec$ alias addr_to_rad, "_OtsAddrToRad"
      integer  addr_to_rad
      external addr_to_rad

cdec$ alias getradid, "_OtsGetRADID"
      integer  getradid
      external getradid

      integer omp_get_thread_num,cpuid,cpu_get_rad
      external omp_get_thread_num,cpuid,cpu_get_rad


#endif
      real*4 densinit(ng1,ng2,ng3)
cmydist densinit(*,block,*)
      common /densi/densinit
      real*4 cr4
#ifndef _ALPHA_KAP_BUG
      cr4(x)=x
#else
      external cr4
#endif

#if defined(_SGI_SOURCE) || defined(_ALPHA)

ccdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(u,phi,def,defp,tmp,tmp2)
c$omp parallel default(private) shared(u,phi,def,defp,tmp,tmp2) 
c make sure data is local using first touch placement:
c*$* assert do (serial)
      do k=1,ng3
c*$* assert do (concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               u(1,i,j,k)=0
               u(2,i,j,k)=0
               u(3,i,j,k)=0
               u(4,i,j,k)=0
               u(5,i,j,k)=0
               phi(i,j,k)=0
               def(i,j,k)=0
               defp(i,j,k)=0
               tmp(i,j,k)=0
               tmp2(i,j,k)=0
               densinit(i,j,k)=0
            enddo
         enddo 
c$omp enddo nowait
      enddo
c$omp end parallel

      write(*,*) 'done first touch'
      loopni=1
      call loadopt
      nzout=0
      write(*,*) 'dump u at z=',(zout(i),i=1,ndimzout)
      open(10,file='state_chk.dat',form='unformatted'
     &     ,status='old',err=75)
      read(10) loopni,a
      close(10)
      do i=1,ndimzout-1
         if (a .ge. .9999/(1+zout(i))) then
            nzout=i
         endif
      enddo
      ai=a
      write(*,*) 'BBKS: using checkpoint data'
      goto 105
75    continue

      write(*,*) 'generating IC from BBKS'

      write (*,100) sigma8, omega0, hubbleh0,
     &              boxsize, redshiftzi
  100 format (1x,'CDM simulation:'
     &      /,7x,'Sigma8=',f10.5,
     &      /,7x,'Omega0= ',f10.3, 3x,'Hubble=',f10.5,
     &      /,7x,'Boxsize=',f10.3, 3x,'z_i=   ',f10.5)

c normal distribution with mean=0 sigma=1
      call parallelrand(def,iseed)
      call hierarch(def,phi,ng1,ng2,ng3)
      call isinc(def,phi,ng1,ng2,ng3)

c$omp parallel default(private) shared(tmp2,def)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=ng1+1,ng1+2
               tmp2(i,j,k)=0
            enddo
         enddo
      enddo

      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               tmp2(i,j,k)=def(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel

c do the FFT to fix any potential problems with symmetries
      call fft3(tmp2,ng1,1)
      
      CALL pk_normalization(anorm,omega0,hubbleh0,sigma8)
      write(*,*) 'CDM normalization value =',anorm

      anorm=anorm*(2*pi/boxsize)**3
      if (redshiftzi .gt. 0) then
         cosm_a=1/(1+redshiftzi)
      else
         cosm_a=0.001
      endif
      if (omega0 .lt. 0.99) then
         call LCDM_grow(grow0, grow1, cosm_a, omega0,omegal)
      else
         write(*,*) 'approximating growth factor to \Omega_0=1'
         grow0=cosm_a
      endif

      write(*,*) 'growth factor relative to \Omega=1:',cosm_a/grow0
      
c anorm scales P(k), so it goes as \delta^2      
      anorm=anorm*grow0**2
      write(*,*) 'start amplitude computation'
cc*$* assert do(concurrent)
c zpj I make the following two line comments 
cc$doacross local(i,j,k,ii,xk,yk,zk,akernel,afreq,qk,x,pk)
cc$&  shared(phi) MP_SCHEDTYPE=INTERLEAVE
c*KAP* parallel region local(i,j,k,ii,xk,yk,zk,akernel,afreq,qk,x,pk)
c*KAP*parallel do

c$omp parallel default(private)
c$omp& firstprivate(boxsize,omega0,hubbleh0,anorm,pk,redshiftzi)
c$omp& shared(phi)
      do i=1,ng1/2+1
c$omp do
         do j=1,i
            do k=1,j
c at this point we only compute one octant, and copy them into
c the full array in the next set of loops
               ii=(i-1)
               xk=2*pi*ii/ng1
               yk=2*pi*(j-1.)/ng2
               zk=2*pi*(k-1.)/ng3
               if (yk .gt. pi) yk=2*pi-yk
               if (zk .gt. pi) zk=2*pi-zk
               akernel=sqrt(xk**2+yk**2+zk**2)
               afreq=akernel*ng3/2/pi
               if (afreq .gt. 0.5) then
                  qk=akernel*ng3/(boxsize*omega0*hubbleh0)
                  call spectrum(pk,qk)
c the factor of 2 is to normalize the sum of two amplitudes.
                  x=sqrt(anorm*pk/2)
#ifndef P3M
                  if (akernel .gt. pi/2 .and. redshiftzi .le. 0) x=0
#endif
                  phi(i,j,k)=x
               endif
            enddo
         enddo
      enddo
c$omp end parallel

c*KAP*end parallel region
      write(*,*)'done amplitude computation'
      phi(1,1,1)=0
!      write(*,'(i2,f20.6)') (i,tmp2(i,i,i),i=1,16)
!      write(*,'(i2,f20.6)') (i,phi(i,i,i),i=1,16)

c$omp parallel default(private) shared(phi,tmp2)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1+2
               io=(i-1)/2+1
               ko=k
               jo=j
               if (ko .gt. ng3/2) ko=ng3+2-k
               if (jo .gt. ng2/2) jo=ng2+2-j
               i1=max(io,jo,ko)
               i3=min(io,jo,ko)
               i2=io+jo+ko-i1-i3
               pk=phi(i1,i2,i3)
               x=tmp2(i,j,k)
c allow for a frequency cutoff               
               tmp2(i,j,k)=x*pk*sqrt(2.*ng1*ng2*ng3)
            enddo
         enddo
      enddo
c$omp end parallel

      tmp2(1,1,1)=ng1*ng2*ng3
      write(*,*) 'done initial seeding'
!      write(*,'(i2,f20.6)') (i,tmp2(i,i,i),i=1,16)
      call fft3(tmp2,ng1,-1)
      write(*,*) 'done FFT'
!      write(*,'(i2,f20.6)') (i,tmp2(i,i,i),i=1,16)

      x=0
c$omp parallel private(i,j,k) shared(tmp2,def)
c$omp& reduction(+:x)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               x=x+tmp2(i,j,k)
               def(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel

      if(abs(ng1*ng2*ng3/x-1).gt. 0.01)then
	write(*,*) 'fft normalization problem'
        write(*,*) 'tmp2(0)=',x,' should be ',ng1*ng2*ng3
        stop
      endif
      if (lrecut) call recut(def,tmp2,phi)
c the old recut subroutine had systematic problems with power aliasing
      if (lrecut .and. .false.) then
         write(*,*) 'using recut and disk data'
         open(10,file='def0.dat',form='unformatted',status='old')
         read(10) def
         read(10) phi
         close(10)
         open(10,file='deltadef.dat',form='unformatted',status='old')
         read(10) defp
         close(10)

c$omp parallel default(private) shared(defp,phi,tmp2,densinit)
         do k=1,ng3
c$omp do
            do j=1,ng2
               do i=1,ng1
                  tmp2(i,j,k)=(defp(i,j,k)+1)*phi(i,j,k)
                  densinit(i,j,k)=phi(i,j,k)
                  defp(i,j,k)=0
                  phi(i,j,k)=0
               enddo
            enddo
         enddo
c$omp end parallel
      endif


      if (redshiftzi .le. 0) then
         sum=0
         rms=0
         rmin=tmp2(1,1,1)
         rmax=rmin
c$omp parallel default(private)
c$omp& shared(tmp2) reduction(max:rmax) reduction(min:rmin)
         do k=1,ng3
c$omp do
            do j=1,ng2
               do i=1,ng1
                  rmax=max(rmax,tmp2(i,j,k))
                  rmin=min(rmin,tmp2(i,j,k))
               enddo
            enddo
         enddo
c$omp end parallel
         ai=0.8*cosm_a/(1-rmin)
         redshiftzi=1/ai-1
         write(*,*) 'using auto z_i=',redshiftzi
         write(*,*) 'Note: using smoothed initial density field'

c$omp parallel default(private) shared(tmp2,cosm_a,ai)
         do k=1,ng3
c$omp do
            do j=1,ng2
               do i=1,ng1
                  tmp2(i,j,k)=(tmp2(i,j,k)-1)*ai/cosm_a+1
               enddo
            enddo
         enddo
c$omp end parallel

         cosm_a=ai
         do i=1,ndimzout-1
            if (ai .ge. .9999/(1+zout(i))) then
               nzout=i
            endif
         enddo
      endif

      sum=0
      rms=0
      rmin=tmp2(1,1,1)
      rmax=rmin
c      open(13,file='rho0.dat',form='unformatted')
c      write(13) tmp2
c      close(13)
#ifdef _ALPHA
	write(*,*) 'magic=',magic
        magic1=magic
       ipagesize=get_page_size()
       write(*,*) 'pagesize=',ipagesize
       write(*,*) 'u(1) on=',addr_to_rad(u(1,1,1,1)),%loc(u(1,1,1,1))
     &  ,mod(%loc(u(1,1,1,1)),ipagesize)
#endif
      call szero(phi,ng1,ng2,ng3)
      call szero(defp,ng1,ng2,ng3)
      call szero(tmp,ng1,ng2,ng3)
c$omp parallel default(private) shared(tmp2,u,phi,defp,tmp)
c$omp&  reduction(+:sum,rms) reduction(max:rmax) reduction(min:rmin)
      do k=1,ng3
c$omp do
         do j=1,ng2
            if (k .eq. 1 .and. mod(j,ng2/32) .eq. 1) 
     &write(*,*) j,addr_to_rad(u(1,1,j,1)),addr_to_rad(tmp2(1,j,1)),
     &   cpu_get_rad()
            do i=1,ng1
               u(1,i,j,k)=tmp2(i,j,k)
               sum=sum+tmp2(i,j,k)
               rms=rms+(tmp2(i,j,k)-1)**2
               rmax=max(rmax,tmp2(i,j,k))
               rmin=min(rmin,tmp2(i,j,k))
               u(2,i,j,k)=0
               u(3,i,j,k)=0
               u(4,i,j,k)=0
               u(5,i,j,k)=1.e-20
               tmp2(i,j,k)=0
c               phi(i,j,k)=0
c               defp(i,j,k)=0
c               tmp(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel
#ifdef _ALPHA
	write(*,*) 'magic=',magic
c        magic=magic1
c      call loadopt
	write(*,*) 'magic=',magic
#endif
      write(*,35) sum/(ng1*ng2*ng3), rmin, rmax, sqrt(rms/(ng1*ng2*ng3))
 35   format('rhomean=',e15.5,' rhomin,max=',2f10.5,' rms=',f10.5)


c we start here on a checkpoint run
      ai=1./((1.+redshiftzi)*5./3.)
105   continue
      af=1.

      loopni=1

      
 135  continue
      nzout=nzout+1
      af=1./(1+zout(nzout))
      cfl=0.5
      if (compressmax .lt. 0) cfl=0.9
      open(12,file='essentials.dat')
      write(12,*) 'ai,af,cfl=',ai,af,cfl
      close(12)
      call runloop(u,phi,tmp,tmp2,def,defp,ai,af,cfl,loopni)
c following line suggested by Jason Prochaska:
      af=1./(1+zout(nzout))
      if (af .gt. 0.99) return
      ai=af
      loopni=-abs(loopni)
      if (nzout .lt. 10) then
      write(scmd,175,err=235) nzout,nzout
 175  format('mkdir Z.',i1
     &      ,'; cp state_chk* u_chk* def_chk* prho.dat Z.',
     &     i1,'/ &')
      else
      write(scmd,175,err=235) nzout,nzout
 185  format('mkdir Z.',i2
     &      ,'; cp state_chk* u_chk* def_chk* prho.dat Z.',
     &     i2,'/ &')
      endif

      call system(scmd)
      if (nzout .lt. ndimzout) goto 135

      return
 235  write(*,*) 'BBKS:i/o error: continuing anyways'
      if (nzout .lt. ndimzout) goto 135
      return
      end



c****************************************************************
      subroutine recut(def,tmp2,t1)
c      
c lay out a dense mesh over the specified region of the grid
c       
      implicit none
      include 'relaxgl.fi'
      real tmp2(ng1+2,ng2,ng3),def(ng1,ng2,ng3),t1(ng1,ng2,ng3)

c locals
      include 'globalpa.fi' 
      real pi,sum,wx,wy,wz,w000,w001,w010
     &     ,w011,w100,w101,w110,w111,sigma2, anorm, cosm_a, grow0, grow1
     &     ,xk,yk,zk,akernel,afreq,qk,pk,sumw, sum1,p1
     &     ,rho000,rho001,rho010,rho011,rho100,rho101,rho110,rho111
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &     ,a23, a33, det, x, y, z, cfact, correctionvolume, detmin
     &     ,weight
      integer icx,icy,icz,nbl,i,j,k,io,jo,ko, idist, ifold(ng3)
     &     ,ii, iopt,i1,i2,i3,ip,im,jp,jm,kp,km, ix, iy, iz
     &     , ixp, iyp, izp, iweight, nsample
      parameter(nsample=3)
      real*4 densinit(ng1,ng2,ng3)
cmydist densinit(*,block,*)
      common /densi/densinit

      
      pi=4*atan(1.)
c detmin is the threshold above which we do any density corrections
      detmin=.9
      
      icx=icrecut(1)
      icy=icrecut(2)
      icz=icrecut(3)
      write(*,91) icx,icy,icz
 91   format(' recut: centering on ',3i5)

      call definit(t1,rcut0,rcompress1,rcompress2)
c now move the regridded region onto the cluster:
      sum=0

c$omp parallel default(private) reduction(+:sum)
c$omp& shared(def,t1,tmp2,icx,icy,icz,densinit) 
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
c suffix "o" is source cell               
               io=i-ng1/2-icx
               jo=j-ng2/2-icy
               ko=k-ng3/2-icz
               io=mod(io+2*ng1-1,ng1)+1
               jo=mod(jo+2*ng2-1,ng2)+1
               ko=mod(ko+2*ng3-1,ng3)+1
               def(i,j,k)=t1(io,jo,ko)
               densinit(i,j,k)=0
               sum=sum+(tmp2(i,j,k)-1)**2
            enddo
         enddo
      enddo
c$omp end parallel

      write(*,*)'recut: intrinsic rms variation='
     &      ,sqrt(sum/(ng1*ng2*ng3))
      sum=0
      correctionvolume=0
c$omp parallel default(private) firstprivate(detmin)
c$omp& shared(def,tmp2,densinit,t1) reduction(+:sum,correctionvolume)
      do k=1,ng3
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
     &           -a12**2*a33
               x=i+(def(ip,j,k)-def(im,j,k))/2
               y=j+(def(i,jp,k)-def(i,jm,k))/2
               z=k+(def(i,j,kp)-def(i,j,km))/2
               if (x .lt. 1) x=x+ng1
               if (y .lt. 1) y=y+ng2
               if (z .lt. 1) z=z+ng3
               ix=x
               iy=y
               iz=z
               wx=x-ix
               wy=y-iy
               wz=z-iz
               ix=mod(ix-1,ng1)+1
               iy=mod(iy-1,ng2)+1
               iz=mod(iz-1,ng3)+1
               w111=wx*wy*wz
               w110=wx*wy*(1-wz)
               w101=wx*(1-wy)*wz
               w100=wx*(1-wy)*(1-wz)
               w011=(1-wx)*wy*wz
               w010=(1-wx)*wy*(1-wz)
               w001=(1-wx)*(1-wy)*wz
               w000=(1-wx)*(1-wy)*(1-wz)
               ixp=mod(ix,ng1)+1
               iyp=mod(iy,ng2)+1
               izp=mod(iz,ng3)+1
               rho000=tmp2(ix,iy,iz)
               rho001=tmp2(ix,iy,izp)
               rho010=tmp2(ix,iyp,iz)
               rho011=tmp2(ix,iyp,izp)
               rho100=tmp2(ixp,iy,iz)
               rho101=tmp2(ixp,iy,izp)
               rho110=tmp2(ixp,iyp,iz)
               rho111=tmp2(ixp,iyp,izp)
               p1=rho000*w000+rho001*w001+rho010*w010
     &           +rho011*w011+rho100*w100+rho101*w101
     &           +rho110*w110+rho111*w111
               t1(i,j,k)=det 
               sum=sum+p1*det
c we want rg0sum to result in zero net fluctuations if p1=0.
               if (det .gt. detmin) 
     &              correctionvolume=correctionvolume+(det-detmin)
      do ko=-nsample,nsample
         do jo=-nsample,nsample
            do io=-nsample,nsample
               iweight=(1+abs(ko/nsample))*(1+abs(jo/nsample))
     &                    *(1+abs(io/nsample))
               weight=1./(2*nsample)**3/iweight
c trapezoidal rule integration over the original domain
               x=i+(def(ip,j,k)-def(im,j,k))/2+io/(2.*nsample)
               y=j+(def(i,jp,k)-def(i,jm,k))/2+jo/(2.*nsample)
               z=k+(def(i,j,kp)-def(i,j,km))/2+ko/(2.*nsample)
               if (x .lt. 1) x=x+ng1
               if (y .lt. 1) y=y+ng2
               if (z .lt. 1) z=z+ng3
               ix=x
               iy=y
               iz=z
               wx=x-ix
               wy=y-iy
               wz=z-iz
               ix=mod(ix-1,ng1)+1
               iy=mod(iy-1,ng2)+1
               iz=mod(iz-1,ng3)+1
               w111=wx*wy*wz
               w110=wx*wy*(1-wz)
               w101=wx*(1-wy)*wz
               w100=wx*(1-wy)*(1-wz)
               w011=(1-wx)*wy*wz
               w010=(1-wx)*wy*(1-wz)
               w001=(1-wx)*(1-wy)*wz
               w000=(1-wx)*(1-wy)*(1-wz)
               ixp=mod(ix,ng1)+1
               iyp=mod(iy,ng2)+1
               izp=mod(iz,ng3)+1
               rho000=tmp2(ix,iy,iz)
               rho001=tmp2(ix,iy,izp)
               rho010=tmp2(ix,iyp,iz)
               rho011=tmp2(ix,iyp,izp)
               rho100=tmp2(ixp,iy,iz)
               rho101=tmp2(ixp,iy,izp)
               rho110=tmp2(ixp,iyp,iz)
               rho111=tmp2(ixp,iyp,izp)
               p1=rho000*w000+rho001*w001+rho010*w010
     &           +rho011*w011+rho100*w100+rho101*w101
     &           +rho110*w110+rho111*w111
c we subtract one to be less prone to rounding error.
               densinit(i,j,k)=densinit(i,j,k)+weight*(p1*det-1)
            enddo
         enddo
      enddo
            enddo
         enddo
      enddo
c$omp end parallel
      cfact=(sum-ng1*ng2*ng3)/correctionvolume
      sum=sum/(ng1*ng2*ng3)
      write(*,*) 'mean subgridded rho=',sum,'  fudge factor=',cfact
c now we use the underdense region to correct the mass excess/deficit
      sum=0

c$omp parallel default(private) firstprivate(detmin,cfact)
c$omp& shared(t1,densinit,tmp2) reduction(+:sum)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               det=t1(i,j,k)
               if ( det .gt. detmin) then
                  densinit(i,j,k) = densinit(i,j,k)-cfact*(det-detmin)
               endif
c this trading makes us a bit less prone to rounding error.  tmp
c is usually double precision, while densinit is single precision.
               tmp2(i,j,k)=densinit(i,j,k)+1
               densinit(i,j,k)=t1(i,j,k)
               sum=sum+tmp2(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel

      sum=sum/(ng1*ng2*ng3)
      write(*,*) 'corrected subgridded rho=',sum
c
      return
      end
      
      

C------------------------------------------------------------------------
      subroutine spectrum(pk, qk)
c 
c     Power-Spectrum function (BBKS):
c     CDM adiabatic fluctuation
c
c if I understand it correctly, qk is in Mpc^-1.
c      
      pk = log(1.0+2.34*qk)/(2.34*qk)
     &   *(1.0+3.89*qk+(16.1*qk)**2+(5.46*qk)**3+(6.71*qk)**4)**(-1./4.)
      pk = qk*pk*pk
      RETURN
      end
      


C------------------------------------------------------------------------
      subroutine pk_normalization(value,omega,hubble,sigma8)
c
c        using the Gaussian Integration, 16th order, to normalize PK
c
      integer modeltype
      real factor
      real ss, sum, xmin,step,h2, xn
      real qn,pk
      real w(0:15),x(0:15)
      data x/0.095012509837637440185d0, 0.281603550779258913230d0,
     &       0.458016777657227386342d0, 0.617876244402643748447d0,
     &       0.755404408355003033895d0, 0.865631202387831743880d0,
     &       0.944575023073232576078d0, 0.989400934991649932596d0,
     &      -0.095012509837637440185d0, -0.281603550779258913230d0,
     &      -0.458016777657227386342d0, -0.617876244402643748447d0,
     &      -0.755404408355003033895d0, -0.865631202387831743880d0,
     &      -0.944575023073232576078d0, -0.989400934991649932596d0/
      data w/0.189450610455068496285d0, 0.182603415044923588867d0,
     &       0.169156519395002538189d0, 0.149595988816576732081d0,
     &       0.124628971255533872052d0, 0.095158511682492784810d0,
     &       0.062253523938647892863d0, 0.027152459411754094852d0,
     &       0.189450610455068496285d0, 0.182603415044923588867d0,
     &       0.169156519395002538189d0, 0.149595988816576732081d0,
     &       0.124628971255533872052d0, 0.095158511682492784810d0,
     &       0.062253523938647892863d0, 0.027152459411754094852d0/
C     write(6,*) 'PK_NORM',omega,hubble
      factor = 4.0*3.1415926535/(8.**3)
      step = atan(1.0)*8.0
      h2 = step/2.0d0
      xmin = 0.0
      sum = 0.0
 10   ss = 0.0
      do i=0,127
          nx = MOD(i,16)
          nn = i/16
          xn = xmin+nn*step+h2*(1.0+x(nx))
          wthx = 3.0*(sin(xn)-xn*cos(xn))/xn**3
          qn = xn/(8.0*omega*hubble)
          call spectrum(pk, qn)
          func = pk*wthx**2*xn**2
          ss = ss + w(nx)*func
      enddo
      xmin = xmin + 8.0*step
      sum = sum + ss*h2
      if (ABS(ss*h2/sum).GE.1.0d-9) goto 10
      sum = sum*factor
      value = sigma8**2/sum
      return
      end



      subroutine lcdm_grow(grow0, grow1, cosm_a, omega0,lambda0)
      real cosm_a, grow0, grow1
      real omega0, lambda0
C-----------------------------------------------------------------------
C  Linear growth factor for CDM+LAMBDA models
C  Equations from Peebles (1980), equation (13.6)
C  INPUT: cosm_a, omega0,lambda0
C  OUTPUT: grow0 (growth rate), grow1 (D grow0 / D cosm_a)
C-----------------------------------------------------------------------
      real w(0:15),x(0:15)
      data x/0.095012509837637440185d0, 0.281603550779258913230d0,
     &       0.458016777657227386342d0, 0.617876244402643748447d0,
     &       0.755404408355003033895d0, 0.865631202387831743880d0,
     &       0.944575023073232576078d0, 0.989400934991649932596d0,
     &      -0.095012509837637440185d0, -0.281603550779258913230d0,
     &      -0.458016777657227386342d0, -0.617876244402643748447d0,
     &      -0.755404408355003033895d0, -0.865631202387831743880d0,
     &      -0.944575023073232576078d0, -0.989400934991649932596d0/
      data w/0.189450610455068496285d0, 0.182603415044923588867d0,
     &       0.169156519395002538189d0, 0.149595988816576732081d0,
     &       0.124628971255533872052d0, 0.095158511682492784810d0,
     &       0.062253523938647892863d0, 0.027152459411754094852d0,
     &       0.189450610455068496285d0, 0.182603415044923588867d0,
     &       0.169156519395002538189d0, 0.149595988816576732081d0,
     &       0.124628971255533872052d0, 0.095158511682492784810d0,
     &       0.062253523938647892863d0, 0.027152459411754094852d0/
C-----------------------------------------------------------------------
      LOGICAL INITIALIZE
      data INITIALIZE/.TRUE./
      SAVE INITIALIZE, grow_z0, a_e, eps



      if (abs(lambda0) .lt. 1.e-8) then
         write(*,*) 'lcdmgrow: using zero lambda, L=',lambda0
         x1 = (1./omega0-1.)
         grow0  = 1.+3./x1+3.*sqrt(1.+x1)/(x1**1.5)
     &        *log(sqrt(1.+x1)-sqrt(x1))
         x1 = (1./omega0-1.)*cosm_a
         grow  = 1.+3./x1+3.*sqrt(1.+x1)/(x1**1.5)
     &        *log(sqrt(1.+x1)-sqrt(x1))
         grow0=grow/grow0
         return
      endif
      
 100  CONTINUE
      IF (INITIALIZE) THEN
        a_e = (0.5*omega0/lambda0)**(1./3.)
        eps = (omega0+lambda0-1.)/(3.*lambda0*a_e*a_e)
        xa = 1./a_e
      ELSE
        xa = cosm_a/a_e
      ENDIF
      step = xa
      h2 = step/2.
      sum = 0.
      do i=0,15
         xn = h2*(1.0+x(i))
         func = xn**1.5*(xn**3-3.*eps*xn+2.)**(-1.5)
         sum = sum + w(i)*func
      enddo
      sum = sum*h2
      d2 = sqrt(xa**3-3*eps*xa+2.)/xa**1.5
      IF (INITIALIZE) THEN
        INITIALIZE = .FALSE.
        grow_z0 = d2*sum
        goto 100
      ELSE
        grow0 = d2*sum/grow_z0
        dd2 = -3.*(1-eps*xa)/(xa**2.5*sqrt(xa**3-3*eps*xa+2))
        func = xa**1.5*(xa**3-3.*eps*xa+2.)**(-1.5)
        grow1 = (dd2*sum + d2*func)/a_e/grow_z0
      ENDIF
      return
      end


#endif      


c the power spectrum P(x) is in our case only the transform of
      real function pspects(x)
      real x, andx

c andx is the power law index that is often associated with "n"
      andx=-1
      pspects=(x+1.e-10)**(andx/2)

      return
      end


      
c*************************************************************************
      subroutine loadopt
      implicit none
      include 'relaxgl.fi'
      include 'globalpa.fi'
#ifdef COLD
      include 'cold.fi'
      integer i,j,k
#endif      
      real*4 densinit(ng1,ng2,ng3)
cmydist densinit(*,block,*)
      common /densi/densinit
      logical checked
      common /crestart/ checked

      checked=.false.
      densinit(1,1,1)=-4321
      omega0=1
      omegab=0
      compressmax=20
      open(10,file='COSMOPAR.DAT',status='old',err=73)
      read(10,*) omega0
      read(10,*) omegab
      read(10,*) compressmax
      read(10,*) hubbleh0
      read(10,*) redshiftzi
      read(10,*) boxsize
      read(10,*) sigma8
      read(10,*) iseed
      read(10,*) omegal
      read(10,*,err=17) lrecut
      read(10,*) icrecut
      read(10,*) rcompress1, rcompress2
      read(10,*) rcut0
      close(10)
      goto 85
 17   continue
      lrecut=.false.
      rcompress1=1
      rcompress2=1
      rcut0=0.5
      goto 85
 73   continue
      write(*,*) 'no file COSMOPAR.DAT found, using default values'
      omega0=1
      omegab=0.048
      compressmax=10
      hubbleh0=0.5
      redshiftzi=100
      boxsize=80
      sigma8=1.05
      iseed(1)=1234
      iseed(2)=1234
      iseed(3)=1234
      iseed(4)=1
      omegal=0
      
 85   continue
#ifdef COLD      
      if (compressmax .lt. 0) then
         write(*,*) 'not using COLD flows without MM'

c$omp parallel default(private) shared(cold)
         do k=1,ng3
c$omp do
            do j=1,ng2
               do i=1,ng1
                  cold(i,j,k)=.false.
               enddo
            enddo
         enddo
c$omp end parallel

      else
         write(*,*) 'loadopt: using cold flow option'
      endif
#endif


      write(*,*) 'loadopt: global parameters are:'
      write(*,100) omega0,omegab,omegal,hubbleh0,sigma8,redshiftzi
 100  format('Omega_[0,b,L]=',3f7.4,'   H_0=',f5.3,'   sigma8=',f6.3
     &     ,'   z_i=',f6.1)
      write(*,110) compressmax, iseed, ng1,ng2,ng3
 110  format('compress lim=',f5.1,' iseed=',4i5,' gridsize=',3i4)
#ifdef NBODY
      write(*,*)' loadopt: using N-body'
#endif
      if (lrecut) then
         write(*,*) 'recutting with parameters:'
         write(*,120) icrecut,rcompress1,rcompress2,rcut0
 120     format(' center=',3i4,' comp1,2=',2f4.1,' r0=',f5.3)
      endif
      magic=-123456
      
      return
      end

      
      
c*************************************************************************
      subroutine runloop(u,phi,tmp,tmp2,def,defp,ai,af,cfl,loopnia)
      implicit none
      include 'relaxgl.fi'
      real u(5,ng1,ng2,ng3),phi(ng1,ng2,ng3),tmp(ng1,ng2,ng3)
     &     ,tmp2(ng1+2,ng2,ng3),def(ng1,ng2,ng3),defp(ng1,ng2,ng3)
cmydist u(*,*,block,*),phi(*,block,*),tmp(*,block,*)
cmydist def(*,block,*),defp(*,block,*),tmp2(*,block,*)
      real ai,af,cfl
      integer loopnia
c locals
      include 'globalpa.fi'
c *** locally declared array
#ifdef NBODY      
       real defold(ng1,ng2,ng3)
#endif      
       integer n,nfluidcmp,i,j,k, isopt, loopni

      parameter(nfluidcmp=5)
      real c,cd,a,tau,dtau,t, da, dt1, dascale, gcmpmax, dtauold
     &     , dtoldcdefp, cdoc, cfluid, aichk
c for use in conformal time
      real dtproj0, tprojold, dtproj, tprojnew, tauconf
      external dascale
      logical lastloop
c checkpoint/restart information: only shared by checkpoint
      common /rlchkpnt/ a,aichk,dtau,tau, tprojold,n
      real*4 cr4
      real x
#ifndef _ALPHA_KAP_BUG
      cr4(x)=x
#else
      external cr4
#endif
#ifdef _ALPHA
       character*24   fdate
       external       fdate
#endif
c loopni is negative on a reentry

      aichk=ai
      if (magic .ne. -123456) then
         write(*,*) 'runloop: global parameter block incorrectly '
     &        ,'initialized, magic=',magic
      stop
      endif
      loopni=abs(loopnia)  
      dtoldcdefp=dtau
      write(*,*) 'runloop: omega0,b, cmax=',omega0,omegab, compressmax
      gcmpmax=compressmax

      lastloop=.false.
c      taufinal=-3/sqrt(af)
      a=ai
      tprojold=tauconf(a,omega0,omegal)
c tell stepghp to initialize
      isopt=1
      if (loopnia .lt. 0) then
         isopt=2
         goto 24
      endif
      dtau=1
      dtoldcdefp=dtau
      tau=-3/sqrt(ai)
      open(10,file='state_chk.dat',status='old',form='unformatted'
     &  ,err=24)
      close(10)
      write(*,*)' using checkpointed data'
      call restart(u,def,defp,phi,tmp)
#ifdef NBODY
c$omp parallel default(private) shared(defold,tmp)
      do k=1,ng3
c$omp do
      do j=1,ng2
      do i=1,ng1
          defold(i,j,k)=tmp(i,j,k)
      enddo
      enddo
      enddo
c$omp end parallel
#endif
      dtoldcdefp=dtau
c setting dtau very small is dangerous in calcdefp, so we
c set it to some large value there.      
      if (a .gt. 0.99) then
         write(*,*) 'last checkpoint was at end of run'
         stop
      endif
      loopni=n+1
      isopt=0
 24   continue

#ifdef GMETRIC
      call gcalcbaj(def)
#endif
c project the simulation once every light crossing time
      dtproj0=boxsize*100/2.998e5

      do n=loopni,10000000
#ifndef _NT
#ifndef _ALPHA
          call system('date')
#else
	write(*,*) fdate()
#endif
#endif
c zpj: make the following line comment for check
c           if (n .ge. 3) stop 
c zpj: check end

         tprojnew=tauconf(a,omega0,omegal)
         dtproj=tprojnew-tprojold
#ifdef PROJOUT
         if (dtproj.ge.dtproj0) then
            call sz_gas(u,def,phi,a)
            call projdm(def)
c we want to have an average writeout of dtproj0
c            tprojold=tprojold+dtproj0
         endif
#endif

         if (gcmpmax .gt. 0) then
c compressmax is a common block variable shared with calcdefp
            compressmax=max(5.,min(gcmpmax, 1.1*gcmpmax*a))
         endif
         dtoldcdefp=min(dtoldcdefp,1.)


         call calcdefp1(defp,tmp,tmp2,def,u,dtoldcdefp,dtau,nfluidcmp)
         call eulercfl(c,cfluid,u,def,defp)
         c=max(1.e-10,c)
         call cfldefp(cd,def,defp)
         cd=max(1.e-10,cd)
         if (cd .gt. c .or. c .gt. 2*cfluid) then
            write(*,*)'cd/c=',cd/c,cd/cfluid
         endif
         c=max(c,cd)
         cdoc=c/cfluid
         dtau=cfl/sqrt(3.)/c
         dtau=min(dtau,1.)
         da=dascale(t,a,dtau)

c we want da/a < 0.1
         if (dtau .gt. abs(a*dtau/da/20)) 
     &       write(*,*)'expansion limited step'
         dtau=min(dtau,abs(a*dtau/da/20))
         da=dascale(t,a,2*dtau)
         dtauold=dtau
         dtoldcdefp=dtau*cdoc
         dtoldcdefp=min(dtoldcdefp,1.)
         if (a+da .gt. af) then
c use bisection rule to find optimal dtau
            dt1=dtau/2
            dtau=dtau-dt1
            do i=1,20
               dt1=dt1/2
               da=dascale(t,a,2*dtau)
               if (a+da .gt. af) then
                  dtau=dtau-dt1
               else
                  dtau=dtau+dt1
               endif
            enddo
            lastloop = .true.
         endif

         da=dascale(t,a,dtau)
         a=a+da
         if (gcmpmax .gt. 0) then
c compressmax is a common block variable shared with calcdefp
            compressmax=max(5.,min(gcmpmax, 1.1*gcmpmax*a))
         endif
c routine stepgh takes a sandwiched sequence: deformation update,
c hydro, gravity*2, hydro deformation.  This moves a total of 2dtau.
c
c         call bremsstrahlung(u,def,a,t,tau,dtau)

#ifdef NBODY
         call stepgh(u,phi,def,defp,defold,tmp2,isopt,a,dtau,dtoldcdefp)
#else
         call stepgh(u,phi,def,defp,tmp,tmp2,isopt,a,dtau,dtoldcdefp)
#endif    
c remember that for N-body integration, stepgh needs tmp to remain
c unchanged between calls.
c         call bremsstrahlung(u,def,a,t,tau,dtau)
#ifdef PROJOUT
         if (dtproj.ge.dtproj0) then
            call projphidot(phi,def,a,da)
c we want to have an average writeout of dtproj0
            tprojold=tprojold+dtproj0
         endif
#endif
         isopt=0
         
         da=dascale(t,a,dtau)
         a=a+da
         tau=tau+2*dtau
         write(*,10)n,t,tau,a,dtau,1/a-1
 10      format(' n=',i4,'  t=',e8.2,' tau=',f8.3,'  a=',f6.3
     &        ,' dtau=',e9.2,' z=',f8.3)
         compressmax=gcmpmax
         if (lastloop) goto 20
#ifdef _SGI_SOURCE
         call adjustthread
#endif

         if (mod(n,200) .eq. 0 .and. ng2 .ge. 64) then
#ifdef NBODY
c$omp parallel default(private) shared(tmp,defold)
             do k=1,ng3
c$omp do
             do j=1,ng2
             do i=1,ng1
                tmp(i,j,k)=defold(i,j,k)
             enddo
             enddo
             enddo
c$omp end parallel
#endif
            call  checkpoint(u,def,defp,phi,tmp)
         endif


         if (mod(n,11) .eq. 0 .and. ng2 .ge. 64) then
         endif
      enddo

 20   continue

c step the particles
      dtau=-1
#ifdef NBODY
       call stepgh(u,phi,def,defp,defold,tmp2,isopt,a,dtau,dtau)
#else
      call stepgh(u,phi,def,defp,tmp,tmp2,isopt,a,dtau,dtau)
#endif 
c we use dtau<0 in the checkpoint file as a signal that we have finished
c the run.      
      dtau=dtauold
c we need to be careful to store a representative dtau, and not
c the synchronization dtau, which could screw up calcdefp.  
#ifdef NBODY
c$omp parallel private(i,j,k) shared(tmp,defold)
       do k=1,ng3
c$omp do
          do j=1,ng2
             do i=1,ng1
                tmp(i,j,k)=defold(i,j,k)
             enddo
          enddo
       enddo
c$omp end parallel
#endif
    
      call checkpoint(u,def,defp,phi,tmp)
#ifdef NBODY    
c$omp parallel default(private) shared(tmp)  
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               tmp(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel

c#ifndef _SGI_SOURCE
c it is rumoured that the following may cause unpredictable segmentation
c faults on f77 7.3 irix 6.5 because tmp2 is dimensioned differently
c in the callee
      call calcrho(tmp)
      open(15,file='prho.dat',form='unformatted',err=235)
      write(15,err=235) (cr4(tmp(i,1,1)),i=1,ng1*ng2*ng3)
      close(15,err=235)
c#endif
#endif
      loopnia=n+1
#ifdef PROJOUT
c      call sz_gas(u,def,a)
#endif
      return
 235  write(*,*) 'runoop:i/o error: continuing anyways'
      return
      end




c*************************************************************************
      subroutine restart(u,def,defp,phi,tmp)
      implicit none
      include 'relaxgl.fi'
c this line is real*4 for conversion purposes      
      real*4 u(5,ng1,ng2,ng3), def(ng1,ng2,ng3), tmp(ng1,ng2,ng3)
     &     ,defp(ng1,ng2,ng3), phi(ng1,ng2,ng3)
     &     ,phiold4(ng1,ng2,ng3)
c the last line is locals^

      integer i,j
c  cumulative date from "runloop"
      real a,dtau, tau,tprojold, ai
      integer n
      common /rlchkpnt/ a,ai,dtau,tau, tprojold,n

c "stepgh"      
      real adotk, adotw, engy0, engy0k, gdtold, dtold
     &     , phiold(ng1,ng2,ng3)
c common block is only shared with checkpoint/restart and stepgh
      common /cstepgh/ phiold, adotk, adotw, engy0, engy0k
     &     , gdtold, dtold
      equivalence (phiold,phiold4)
#ifdef COLD
      include 'cold.fi'
#endif      
      real*4 densinit(ng1,ng2,ng3)
cmydist densinit(*,block,*)
      common /densi/densinit
      logical checked
      common /crestart/ checked
     
#ifdef _ALPHA
       character*24   fdate
       external       fdate
#endif


#ifdef NBODY
      include 'nbody.fi'
      
c from stepgh:
      
      real aold, defold(ng1,ng2,ng3)
c this line is real*4 for conversion purposes      
      real*4 defold4(ng1,ng2,ng3), xv1(6,npart), pmass1(npart)
c the common block is only shared with checkpoint and restart      
      common /defold/ aold
      equivalence(xv1,xv)
      equivalence(pmass1,pmass)
      
#ifdef NBODY
c$omp parallel do default(none) shared(xv) 
      do i=1,npart
         xv(1,i)=0
      enddo
c$omp end parallel do
#endif
#ifdef _ALPHA
      write(*,*) 'starting restart nbody ',fdate()
#endif
      open(15,file='rawpart_chk.dat',form='unformatted')
c      read(15) ((xv1(j,i),j=1,6),i=1,npart)
      read(15) xv1
      close(15)
#ifdef _ALPHA
      write(*,*) 'done reading xv ',fdate()
#endif
      open(15,file='pmass0.dat',form='unformatted',err=110)
      read(15) pmass1
      close(15)
      goto 120
 110  continue
      write(*,*) 'no pmass0.dat found, using pmass=1'
      do i=1,npart
         pmass1(i)=1
      enddo
 120  continue
      open(15,file='densinit0.dat',form='unformatted',err=130)
      read(15) densinit
      close(15)
      goto 140
 130  continue
      write(*,*) 'no densinit0.dat found, using 1'
      do i=1,ng1*ng2*ng3
         densinit(i,1,1)=1
      enddo
 140  continue

      if (xv1(1,1) .ne. xv(1,1) ) then
         write(*,*) 'restart: expanding particles to real*8'
         call expandr48(xv,xv1,6*npart)
         call expandr48(pmass,pmass1,npart)
      else
         write(*,*) 'restart: keeping real*4 particle data'
      endif
      open(15,file='defold_chk.dat',form='unformatted')
c      read(15) (defold4(i,1,1),i=1,ng1*ng2*ng3)
      read(15) defold4
      read(15) aold
      close(15)
      call expandr48(tmp,defold4,ng1*ng2*ng3)

      call makechainparallel
      
#endif

#ifdef _ALPHA
      write(*,*) 'starting restart gas',fdate()
#endif
      open(15,file='u_chk.dat',form='unformatted')
c      read(15) (u(i,1,1,1),i=1,5*ng1*ng2*ng3)
      read(15) u
      close(15)
      call expandr48(u,u,5*ng1*ng2*ng3)
      open(15,file='def_chk.dat',form='unformatted')
c      read(15) (def(i,1,1),i=1,ng1*ng2*ng3)
      read(15) def
      close(15)
      call expandr48(def,def,ng1*ng2*ng3)
      open(15,file='defp_chk.dat',form='unformatted')
c      read(15) (defp(i,1,1),i=1,ng1*ng2*ng3)
      read(15) defp
      close(15)
      call expandr48(defp,defp,ng1*ng2*ng3)
      open(15,file='phiold_chk.dat',form='unformatted')
c      read(15) (phiold4(i,1,1),i=1,ng1*ng2*ng3)
      read(15) phiold4
      close(15)
      call expandr48(phiold,phiold4,ng1*ng2*ng3)
      open(15,file='phi_chk.dat',form='unformatted')
c      read(15) (phi(i,1,1),i=1,ng1*ng2*ng3)
      read(15) phi
      close(15)
      call expandr48(phi,phi,ng1*ng2*ng3)
      open(15,file='state_chk.dat',form='unformatted')
      read(15) n,a,dtau,tau
      read(15) adotk, adotw, engy0, engy0k, gdtold, dtold, tprojold
      close(15)
#ifdef COLD
      open(15,file='cold_chk.dat',form='unformatted')
      read(15) (cold(i,1,1),i=1,ng1*ng2*ng3)
      close(15)
#endif      
      checked=.true.

#ifdef _ALPHA
      write(*,*) 'done restart',fdate()
#endif
      return
      end


      
      subroutine expandr48(b,a,n)
      implicit none
      integer n
      real*4 a(n)
      real b(n)
c local
      integer i

      
c*$* assert do(serial)
      do i=n,1,-1
         b(i)=a(i)
      enddo
      return
      end


#ifdef _ALPHA_KAP_BUG
c normally, cr4 is a statement function, but the alpha KF77 processor
c optimizes it away.

      real*4 function cr4(x)
      real x

      cr4=x
      return
      end
#endif
      

      
c*************************************************************************
      subroutine checkpoint(u,def,defp,phi,tmp)
c currently, tmp(,,) doesnt need to be saved      
      implicit none
      include 'relaxgl.fi'
      real u(5,ng1,ng2,ng3), def(ng1,ng2,ng3), tmp(ng1,ng2,ng3)
     &     ,defp(ng1,ng2,ng3),phi(ng1,ng2,ng3)

c locals
      integer i,j
      
c  cumulative data from "runloop"
      real a, tau,dtau,tprojold, ai
      integer n
      common /rlchkpnt/ a,ai,dtau,tau, tprojold,n
      logical firsttime, checked
      data firsttime /.true./
      common /crestart/ checked
c the crucial feature here is the calling order:  loadopt sets
c checked to false, and loadopt must be called before either checkpoint
c or restart.  If restart was called before checkpoint, we do not checkpoint
c the invariant data, i.e. the particle masses.
      real*4 densinit(ng1,ng2,ng3)
cmydist densinit(*,block,*)
      common /densi/densinit

c "stepgh"      
      real adotk, adotw, engy0, engy0k, gdtold, phiold(ng1,ng2,ng3)
     &     , dtold
c common block is only shared with checkpoint/restart      
      common /cstepgh/ phiold, adotk, adotw, engy0, engy0k
     &     , gdtold, dtold
     

      real x
c function to convert to 4 byte real from any size
#ifdef COLD
      include 'cold.fi'
#endif      

#ifdef NBODY
      include 'nbody.fi'
#ifdef PROJOUT
c      real*4 xve(6,npmassive)
#endif
      
c from stepgh:
      
      real aold,defold(ng1,ng2,ng3)
c the common block is only shared with checkpoint and restart      
      common /defold/ aold

#endif
      
      real*4 cr4
#ifndef _ALPHA_KAP_BUG
      cr4(x)=x
#else
      external cr4
#endif

      
      open(15,file='u_chk',form='unformatted',err=900)
      write(15,err=900) (cr4(u(i,1,1,1)),i=1,5*ng1*ng2*ng3)
      close(15,err=900)
      open(15,file='def_chk',form='unformatted',err=900)
      write(15,err=900) (cr4(def(i,1,1)),i=1,ng1*ng2*ng3)
      close(15,err=900)
      open(15,file='defp_chk',form='unformatted',err=900)
      write(15,err=900) (cr4(defp(i,1,1)),i=1,ng1*ng2*ng3)
      close(15,err=900)
      open(15,file='phi_chk',form='unformatted',err=900)
      write(15,err=900) (cr4(phi(i,1,1)),i=1,ng1*ng2*ng3)
      close(15,err=900)
      open(15,file='phiold_chk',form='unformatted',err=900)
      write(15,err=900) (cr4(phiold(i,1,1)),i=1,ng1*ng2*ng3)
      close(15,err=900)
      open(15,file='state_chk',form='unformatted',err=900)
      write(15,err=900) n,a,dtau,tau
      write(15,err=900) adotk, adotw, engy0, engy0k, gdtold, dtold
     &     , tprojold
      close(15,err=900)
#ifdef COLD
      open(15,file='cold_chk',form='unformatted',err=900)
      write(15,err=900) (cold(i,1,1),i=1,ng1*ng2*ng3)
      close(15,err=900)
#endif      
      if (.not. checked) then
         write(*,*) 'checking densinit field'
         open(15,file='densinit0.dat',form='unformatted',err=900)
         write(15,err=900) (densinit(i,1,1),i=1,ng1*ng2*ng3)
         close(15,err=900)
      endif
#ifdef NBODY
      if (.not. checked) then
         open(15,file='pmass0.dat',form='unformatted',err=900)
         write(15,err=900) (pmass(i),i=1,npart)
         write(15,err=900) ai
         close(15,err=900)
      endif
      open(15,file='rawpart_chk',form='unformatted',err=900)
c we usually keep particles in real*4, so no need to convert
      write(15,err=900) ((xv(j,i),j=1,6),i=1,npart)
      close(15,err=900)
      open(15,file='defold_chk',form='unformatted',err=900)
      write(15,err=900) (cr4(tmp(i,1,1)),i=1,ng1*ng2*ng3)
      write(15,err=900) aold
      close(15,err=900)
      
#ifdef PROJOUT
c      call xvmap(xv,xve,def,npmassive,ng1,a)
c      call sz_gas(u,def,a)
c      call sz_proj(xve,pmass,npmassive,a)
#endif
#endif

      checked=.true.
      call system('for n in *_chk ; do mv $n $n.dat ; done ')
      write(*,*) 'checkpoint: successful'
      
      return
900   write(*,*) 'checkpoint: I/O failure, continuing anyways.'
      close(15,err=910)
      return
 910  write(*,*) 'recovery close failed as well'
      return
      end



      subroutine ustat(u)
      implicit none
      include 'relaxgl.fi'
      real u(5,ng1,ng2,ng3)

      integer i,j,k,mni(ng3),mnj(ng3),mxi(ng3),mxj(ng3)
     &   ,mni1,mnj1,mnk1,mxi1,mxj1,mxk1
      real umin(ng3),umax(ng3),usum, urms, u1, umax1,umin1

      usum=0
      urms=0
c we have to manually do the MAXLOC reduction

cc$omp parallel default(private)
cc$omp& shared(u) reduction(+:usum,urms) reduction(max:umax,mni,mnj)
cc$omp& reduction(min:umin,mni,mnj)
      do k=1,ng3
         umin(k)=1000
         umax(k)=-100
cc$omp do
         do j=1,ng2
            do i=1,ng1
              u1=u(1,i,j,k)
               if (umax(k) .lt. u1) then
                  mxi(k)=i
                  mxj(k)=j
                  umax(k)=u1
               endif
               if (umin(k) .gt. u1) then
                  mni(k)=i
                  mnj(k)=j
                  umin(k)=u1
               endif
               usum=usum+u1
               urms=urms+(u1-1)**2
            enddo
         enddo
      enddo
cc$omp end parallel

      umax1=umax(1)-1
      umin1=umin(1)+1
      do k=1,ng3
         if (umax1 .lt. umax(k) ) then
            umax1=umax(k)
            mxi1=mxi(k)
            mxj1=mxj(k)
            mxk1=k

            mxj1=mxj(k)
            mxk1=k
         endif
         if (umin1 .gt. umin(k) ) then
            umin1=umin(k)
            mni1=mni(k)
            mnj1=mnj(k)
            mnk1=k
         endif
      enddo
      write(*,75) umin1,umax1,usum/(ng1*ng2*ng3)
     &             ,sqrt(urms/(ng1*ng2*ng3))
 75   format(' umin/max, mean, rms=',4e14.6)
      write(*,85) mni1,mnj1,mnk1,mxi1,mxj1,mxk1
 85   format(' min,max pos=',6i5)
      write(*,95) (u(i,mni1,mnj1,mnk1),i=1,5)
 95   format(' u(min)=',5e12.4)

      return
      end




















