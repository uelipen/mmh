c  Output xyz projected density field
c adapted for mmh by Ue-Li 11/13/97
        subroutine sz_proj(xv,amass,np,a)
        implicit none
        include 'globalpa.fi'
        include 'relaxgl.fi'
        integer np
        real*4 xv(6,np),amass(np)
        real*8 a

! locals
        integer n,nrhoim
        parameter (n=64,nrhoim=1)
        real*4 rho(n,n,n,7)
        common /sz/ rho
        real*4 dx,omegam,omegav,h0,epsilon,aout
        integer i,j,k,n1sz,n2sz,n3sz,i1,i2,i3,ni1,ni2,ni3,l
        parameter (n1sz=n,n2sz=n,n3sz=n)
        real*4 xy(n1sz,n2sz,nrhoim),yz(n2sz,n3sz,nrhoim)
     1  ,xz(n1sz,n3sz,nrhoim)
        real*4 xy1(n1sz,n2sz,nrhoim),yz1(n2sz,n3sz,nrhoim)
     1  ,xz1(n1sz,n3sz,nrhoim)
        real*4 xy2(n1sz,n2sz,nrhoim),yz2(n2sz,n3sz,nrhoim)
     1  ,xz2(n1sz,n3sz,nrhoim)
        real*4 xy3(n1sz,n2sz,nrhoim),yz3(n2sz,n3sz,nrhoim)
     1  ,xz3(n1sz,n3sz,nrhoim)
        real me,conksz,consz,yhe,ne,te,conx,dxscale,sig1,sig2,sig3
        real rhorhobar,xpart,vmean1,vmean2,vmean3,v2mean1,v2mean2
        real v2mean3,masstot,dd1,dd2,dd3,am
        logical firsttime
        save firsttime
        data firsttime /.true./
        

        aout=a
        dx=boxsize/hubbleh0/n
        if (firsttime) then
           firsttime=.false.
           h0=hubbleh0
           omegam=omega0
           omegav=0
           open(21,file='szdm.dat',status='unknown',form='unformatted')
           rewind 21
           write(21)n,n,n,nrhoim
           write(21)dx,omegam,omegav,h0
         open(22,file='kinszdm.dat',status='unknown',form='unformatted')
           rewind 22
           write(22)n,n,n,nrhoim
           write(22)dx,omegam,omegav,h0
         open(23,file='xraydm.dat',status='unknown',form='unformatted')
           rewind 23
           write(23)n,n,n,nrhoim
           write(23)dx,omegam,omegav,h0
         open(24,file='tempdm.dat',status='unknown',form='unformatted')
           rewind 24
           write(24)n,n,n,nrhoim
           write(24)dx,omegam,omegav,h0
        endif
           
c conversion factors
        yhe=0.24
c optical depth assuming He ionized
        ne=(1-yhe/2.0)*2.3e-5
c velocity in v/c: conversion factor for kinetic sz: v/c*tau
        conksz=1./2.998e5*ne
c electron T 0.6mp(v/c)**2 or (400km/s)**2=1/3keV: multiply with 3-d 
c velocity dispersion
        te=1./400.0**2/3.0*1000.0
        me=511000.0 
c Xray conversion factor in units of ergs/cm^2/s
        conx=5.9e-13*sqrt(te)
        dxscale=n*1./ng1
c SZ conversion factor: 2y=2taukT/me multiplied with dx
c (1+z)**2omega_bh^2 to be put in later
        consz=2*te/me*ne
c initialize
c$omp parallel do default(private)
c$omp&        shared(xy,xy1,xy2,xy3,yz,yz1,yz2,yz3,xz,xz1,xz2,xz3,rho
c$omp&       ,amass,xv, dxscale,te,consz,ni1,ni2,ni3)
        do j=1,n
           do i=1,n
              xy(i,j,ni3)=0.0
              xy1(i,j,ni3)=0.0
              xy2(i,j,ni3)=0.0
              xy3(i,j,ni3)=0.0
              yz(i,j,ni1)=0.0
              yz1(i,j,ni1)=0.0
              yz2(i,j,ni1)=0.0
              yz3(i,j,ni1)=0.0
              xz(i,j,ni2)=0.0
              xz1(i,j,ni2)=0.0
              xz2(i,j,ni2)=0.0
              xz3(i,j,ni2)=0.0
           enddo
        enddo
       do l=1,7
       do k=1,n
c$omp do
        do j=1,n
           do i=1,n
              rho(i,j,k,l)=0
           enddo
           enddo
        enddo
        enddo
       masstot=0.0
c$omp do
       do j=1,np
          dd1=xv(1,j)*dxscale
          dd2=xv(2,j)*dxscale
          dd3=xv(3,j)*dxscale
c  Nearest grid point, different than CIC or TSC.
          i1=dd1+0.5
          i2=dd2+0.5
          i3=dd3+0.5
          if (i1.eq.0) i1=n1sz
          if (i2.eq.0) i2=n2sz
          if (i3.eq.0) i3=n3sz
c          am=amass(j)
c rms shouldnt depend on particle mass
          am=1.
          masstot=masstot+amass(j)
          rho(i1,i2,i3,1)=rho(i1,i2,i3,1)+am
          rho(i1,i2,i3,2)=rho(i1,i2,i3,2)+xv(4,j)*am
          rho(i1,i2,i3,3)=rho(i1,i2,i3,3)+xv(5,j)*am
          rho(i1,i2,i3,4)=rho(i1,i2,i3,4)+xv(6,j)*am
         rho(i1,i2,i3,5)=rho(i1,i2,i3,5)+(xv(4,j)**2+xv(4,j)**2
     &             +xv(5,j)**2)*am
         rho(i1,i2,i3,6)=rho(i1,i2,i3,6)+amass(j)
       enddo
       do  i3=1,n3sz
          ni3=(i3-1)*nrhoim/n3sz+1
          do i2=1,n2sz
            ni2=(i2-1)*nrhoim/n2sz+1
            do i1=1,n1sz
              ni1=(i1-1)*nrhoim/n1sz+1
              if (rho(i1,i2,i3,1).gt.0.001) then 
              xpart=rho(i1,i2,i3,1)
              am=rho(i1,i2,i3,6)
              vmean1=rho(i1,i2,i3,2)/rho(i1,i2,i3,1)
              vmean2=rho(i1,i2,i3,3)/rho(i1,i2,i3,1)
              vmean3=rho(i1,i2,i3,4)/rho(i1,i2,i3,1)
              v2mean1=rho(i1,i2,i3,5)/rho(i1,i2,i3,1)
c make it unbiased
              if (xpart .gt. 1.1) then
                 sig2=xpart/(xpart-1.0)*(v2mean1-
     1                vmean1**2-vmean2**2-vmean3**2)
              else
                 sig2=(v2mean1-
     1                vmean1**2-vmean2**2-vmean3**2)/2
c just say mach number=1
              endif
              sig1=sqrt(abs(sig2))
              sig3=sig1*sig2
              rhorhobar=am/masstot*n1sz*n2sz*n3sz
              xy(i1,i2,ni3)=xy(i1,i2,ni3)+rhorhobar*sig2*consz
              yz(i2,i3,ni1)=yz(i2,i3,ni1)+rhorhobar*sig2*consz
              xz(i1,i3,ni2)=xz(i1,i3,ni2)+rhorhobar*sig2*consz
              xy1(i1,i2,ni3)=xy1(i1,i2,ni3)+rhorhobar*vmean3*conksz
              yz1(i2,i3,ni1)=yz1(i2,i3,ni1)+rhorhobar*vmean1*conksz
              xz1(i1,i3,ni2)=xz1(i1,i3,ni2)+rhorhobar*vmean2*conksz
              xy2(i1,i2,ni3)=xy2(i1,i2,ni3)+rhorhobar**2*sig1*conx
              yz2(i2,i3,ni1)=yz2(i2,i3,ni1)+rhorhobar**2*sig1*conx
              xz2(i1,i3,ni2)=xz2(i1,i3,ni2)+rhorhobar**2*sig1*conx
              xy3(i1,i2,ni3)=xy3(i1,i2,ni3)+rhorhobar**2*sig3*conx*te
              yz3(i2,i3,ni1)=yz3(i2,i3,ni1)+rhorhobar**2*sig3*conx*te
              xz3(i1,i3,ni2)=xz3(i1,i3,ni2)+rhorhobar**2*sig3*conx*te
             endif
            enddo
          enddo
        enddo
c write output
        write(21)aout
        write(21)(((xy(k,j,i), k=1,n2sz), j=1,n1sz), i=1,nrhoim)
        write(21)(((yz(k,j,i), k=1,n3sz), j=1,n2sz), i=1,nrhoim)
        write(21)(((xz(k,j,i), k=1,n3sz), j=1,n1sz), i=1,nrhoim)
        write(22)aout
        write(22)(((xy1(k,j,i), k=1,n2sz), j=1,n1sz), i=1,nrhoim)
        write(22)(((yz1(k,j,i), k=1,n3sz), j=1,n2sz), i=1,nrhoim)
        write(22)(((xz1(k,j,i), k=1,n3sz), j=1,n1sz), i=1,nrhoim)
        write(23)aout
        write(23)(((xy2(k,j,i), k=1,n2sz), j=1,n1sz), i=1,nrhoim)
        write(23)(((yz2(k,j,i), k=1,n3sz), j=1,n2sz), i=1,nrhoim)
        write(23)(((xz2(k,j,i), k=1,n3sz), j=1,n1sz), i=1,nrhoim)
        write(24)aout
        write(24)(((xy3(k,j,i), k=1,n2sz), j=1,n1sz), i=1,nrhoim)
        write(24)(((yz3(k,j,i), k=1,n3sz), j=1,n2sz), i=1,nrhoim)
        write(24)(((xz3(k,j,i), k=1,n3sz), j=1,n1sz), i=1,nrhoim)
        call wimage(xy,n2sz,n1sz,'out')
      return
      end
      

c now for gas:
      subroutine sz_gas(u,def,phi,a)
!      implicit none
      include 'globalpa.fi'
      include 'relaxgl.fi'
      include 'proj.fi'
      real u(5,ng1,ng2,ng3),def(ng1,ng2,ng3),phi(ng1,ng2,ng3)
      real*8 a

! locals
      integer n,nrhoim
      parameter (n=ndimproj,nrhoim=1)
c      common /sz/ rho
      real*4 dx,omegam,omegav,h0,epsilon,aout,sigma8l
      integer i,j,k,n1sz,n2sz,n3sz,i1,i2,i3,ni1,ni2,ni3
      parameter (n1sz=n,n2sz=n,n3sz=n)
      real*4 xy(n1sz,n2sz,nrhoim),yz(n2sz,n3sz,nrhoim)
     1     ,xz(n1sz,n3sz,nrhoim)
      real*4 xy1(n1sz,n2sz,nrhoim),yz1(n2sz,n3sz,nrhoim)
     1     ,xz1(n1sz,n3sz,nrhoim)
      real*4 xy2(n1sz,n2sz,nrhoim),yz2(n2sz,n3sz,nrhoim)
     1     ,xz2(n1sz,n3sz,nrhoim)
      real*4 xy3(n1sz,n2sz,nrhoim),yz3(n2sz,n3sz,nrhoim)
     1     ,xz3(n1sz,n3sz,nrhoim)
      real*4 xy4(n1sz,n2sz,nrhoim),yz4(n2sz,n3sz,nrhoim)
     1     ,xz4(n1sz,n3sz,nrhoim)
      real*4 xy5(n1sz,n2sz,nrhoim),yz5(n2sz,n3sz,nrhoim)
     1     ,xz5(n1sz,n3sz,nrhoim)
      equivalence (xy,yz),(xy,xz)
      equivalence (xy1,yz1),(xy1,xz1)
      equivalence (xy2,yz2),(xy2,xz2)
      equivalence (xy3,yz3),(xy3,xz3)
      equivalence (xy4,yz4),(xy4,xz4)
      equivalence (xy5,yz5),(xy5,xz5)
      common /phidot/ xy5,rcut
      real me,conksz,consz,yhe,ne,te,conx,dxscale,sig1,sig2,sig3
      real rhorhobar,xpart,vmean1,vmean2,vmean3,v2mean1,v2mean2
      real v2mean3,masstot,dd1,dd2,dd3,am
      real*8 checksum, checksum1, checksum2, checksum3
      integer icount
      logical firsttime
      save firsttime,icount
      data firsttime /.true./
      data icount /0/
c     coord stuff:
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13,a22
     &     , a23, a33, det,x,y,z,therm,alx,pi,alambda,tau,rhoc,amu
     &     ,engyu, velfact
      integer ip,jp,im,jm,kp,km
      
      nsub=n/ng1
c the size of each elliptical spheroid
      rcut=sqrt(2.)
      if (nsub*ng1 .ne. n) then
         write(*,*) 'sz_gas: n not muliple of ng1:',n,ng1
         stop
      endif
      icount=icount+1
      aout=a
      dx=boxsize/hubbleh0/n
      if (firsttime) then
         firsttime=.false.
         h0=hubbleh0
         omegam=omega0
         omegav=omegal
         sigma8l=sigma8
         open(25,file='szgas.dat',status='new',form='unformatted'
     &        ,err=20)
         rewind 25
         write(25)n,n,n,nrhoim
         write(25)dx,omegam,omegav,h0,sigma8l
         open(26,file='kinszgas.dat',status='new',form='unformatted'
     &        )
         rewind 26
         write(26)n,n,n,nrhoim
         write(26)dx,omegam,omegav,h0,sigma8l
         open(27,file='xraygas.dat',status='new',form='unformatted')
         rewind 27
         write(27)n,n,n,nrhoim
         write(27)dx,omegam,omegav,h0,sigma8l
         open(28,file='tempgas.dat',status='new',form='unformatted')
         rewind 28
         write(28)n,n,n,nrhoim
         write(28)dx,omegam,omegav,h0,sigma8l
        open(29,file='deltagas.dat',status='new',form='unformatted')
         rewind 29
         write(29)n,n,n,nrhoim
         write(29)dx,omegam,omegav,h0,sigma8l
        open(30,file='phidot.dat',status='new',form='unformatted')
         rewind 30
         write(30)n,n,n,nrhoim
         write(30)dx,omegam,omegav,h0,sigma8l
        open(107,file='dm.dat',status='new',form='unformatted')
         rewind 107
         write(107)n,n,n,nrhoim
         write(107)dx,omegam,omegav,h0,sigma8l
         goto 30
 20      continue
         open(25,file='szgas.dat',status='old',form='unformatted'
     &        ,access='append')
!     &        ,position='append')
         open(26,file='kinszgas.dat',status='old',form='unformatted'
     &        ,access='append')
         open(27,file='xraygas.dat',status='old',form='unformatted'
     &        ,access='append')
         open(28,file='tempgas.dat',status='old',form='unformatted'
     &        ,access='append')
         open(29,file='deltagas.dat',status='old',form='unformatted'
     &        ,access='append')
         open(30,file='phidot.dat',status='old',form='unformatted'
     &        ,access='append')
         open(107,file='dm.dat',status='old',form='unformatted'
     &        ,access='append')
 30      continue
      endif
      
c conversion factors
      yhe=0.24
c     optical depth assuming He ionized
      ne=(1-yhe/2.0)*2.3e-5
c     velocity in v/c: conversion factor for kinetic sz: v/c*tau
      conksz=1./2.998e5*ne
c     electron T 0.6mp(v/c)**2 or (400km/s)**2=1/3keV: multiply with 3-d 
c     velocity dispersion
      te=1./400.0**2/3.0*1000.0
c     Xray conversion factor in units of ergs/cm^2/s
      conx=5.9e-13*sqrt(te)
      
      dxscale=n*1./ng1
c     SZ conversion factor: 2y=2taukT/me multiplied with dx
c     (1+z)**2omega_bh^2 to be put in later
      me=511000.0 
      consz=2*te/me*ne
      pi=4*atan(1.)
      alambda=boxsize/hubbleh0/ng1*3.0856E22
      rhoc=hubbleh0**2 * 1.88E-26
      amu=omegab*rhoc*alambda**3
      tau=sqrt(omegab*alambda**3/(omega0*amu*6*pi*6.6726E-11 ))
      engyu=amu*alambda**2/tau**2
! to get into (proper) km/sec:
      velfact=alambda/(tau*1000*a)
      
      
      ni1=1
      ni2=1
      ni3=1
c     initialize
c$omp parallel do private(i,j) 
c$omp& shared(xy,xy1,xy2,xy3,xy4,xy5,yz,yz1,yz2,yz3,yz4,yz5,xz,xz1,xz2
c$omp&     ,xz3,xz4,xz5,ni1,ni2,ni3)
        do j=1,n
           do i=1,n
              xy(i,j,ni3)=0.0
              xy1(i,j,ni3)=0.0
              xy2(i,j,ni3)=0.0
              xy3(i,j,ni3)=0.0  
              xy4(i,j,ni3)=0.0
              xy5(i,j,ni3)=0.0
              yz(i,j,ni1)=0.0
              yz1(i,j,ni1)=0.0
              yz2(i,j,ni1)=0.0
              yz3(i,j,ni1)=0.0
              yz4(i,j,ni1)=0.0
              yz5(i,j,ni1)=0.0
              xz(i,j,ni2)=0.0
              xz1(i,j,ni2)=0.0
              xz2(i,j,ni2)=0.0
              xz3(i,j,ni2)=0.0
              xz4(i,j,ni2)=0.0
              xz5(i,j,ni2)=0.0
           enddo
        enddo
c$omp end parallel do
      masstot=ng1*ng2*ng3
      checksum=0
c the parallelization depends on the projection direction
      if (mod(icount,3) .ne. 1) then
c if mod(icount,3) = 0 project down z
c                  = 1              x
c                  = 2              y
c*$* assert do (concurrent)
c$omp parallel do default(private) 
c$omp&         shared(def,dxscale,u,velfact,conx,conksz,phi,a
c$omp&        ,icount,nsub,rcut,xy,xy1,xy2,xy3,xy4,xy5,yz,yz1,yz2,yz3
c$omp&        ,yz4,yz5,xz,xz1,xz2,xz3,xz4,xz5,masstot,consz,te)
c$omp& reduction(+:checksum) 
      do k=1,ng3
c parallelization relies on the fact that the tiling is
c very local:  cells only overlap with their nearest neighbors,
c so on any given tier, it is safe to do physically separated ones
c simultaneously.
c
         do j=1,ng2
            do i=1,ng1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
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
     &              -a12**2*a33
               x=i+(def(ip,j,k)-def(im,j,k))/2
               y=j+(def(i,jp,k)-def(i,jm,k))/2
               z=k+(def(i,j,kp)-def(i,j,km))/2
               x=mod(x+ng1,1.*ng1)
               y=mod(y+ng2,1.*ng2)
               z=mod(z+ng3,1.*ng3)
               
               dd1=x*dxscale
               dd2=y*dxscale
               dd3=z*dxscale
c     Nearest grid point, different than CIC or TSC.
               i1=nint(dd1)
               i2=nint(dd2)
               i3=nint(dd3)
               if (i1.le.0) i1=i1+n1sz
               if (i2.le.0) i2=i2+n2sz
               if (i3.le.0) i3=i3+n3sz
               am=u(1,i,j,k)
               therm=(u(5,i,j,k)-(u(2,i,j,k)**2+u(3,i,j,k)**2+
     &             u(4,i,j,k)**2)/u(1,i,j,k)/2)*
     &              velfact**2
               vmean1=u(2,i,j,k)*velfact/am
               vmean2=u(3,i,j,k)*velfact/am
               vmean3=u(4,i,j,k)*velfact/am
               sig2=therm/am
               sig1=sqrt(abs(sig2))
               sig3=sig1*sig2
               rhorhobar=am/masstot*n1sz*n2sz*n3sz
               ni1=(i1-1)*nrhoim/n1sz+1
               ni2=(i2-1)*nrhoim/n2sz+1
               ni3=(i3-1)*nrhoim/n3sz+1
               alx=conx*am**2*sqrt(abs(therm/am
     &              ))/det/n1sz
               tsz=rhorhobar*sig2*consz/n1sz
c reminder:
c icount=0,2 ->               
               if (mod(icount,3) .eq. 0) then
                  tksz=rhorhobar*vmean3*conksz/n1sz
               else
                  tksz=rhorhobar*vmean2*conksz/n1sz
               endif
	       ttalx=alx*sig1*te/n1sz
               tm=am*n1sz*n2sz/(ng1*ng2*ng3)
               phi2=phi(i,j,k)*det/n1sz
               checksum=checksum+tm

               b11=(a22*a33-a23**2)/det
               b12=(a13*a23-a12*a33)/det
               b13=(a12*a23-a13*a22)/det
               b22=(a11*a33-a13**2)/det
               b23=(a12*a13-a11*a23)/det
               b33=(a11*a22-a12**2)/det
               s11=(b11**2+b12**2+b13**2)
               s12=(b11*b12+b12*b22+b13*b23)
               s13=(b11*b13+b12*b23+b13*b33)
               s22=(b12**2+b22**2+b23**2)        
               s23=(b12*b13+b22*b23+b23*b33)
               s33=(b13**2+b23**2+b33**2)
               
               dx=(def(ip,j,k)-def(im,j,k))/2
               dy=(def(i,jp,k)-def(i,jm,k))/2
               dz=(def(i,j,kp)-def(i,j,km))/2
               if (mod(icount,3) .eq. 0) then
                  w11=s11 - s13**2/s33
                  w12=s12 - s13*s23/s33
                  w22=s22 - s23**2/s33

! we apply the same procedure again to find the bounding boxes

                  xmax=1/sqrt(w11-w12**2/w22)
                  ymax=1/sqrt(w22-w12**2/w11)
                  imax=xmax*nsub*rcut+1
                  jmax=ymax*nsub*rcut+1
                  narea=0
                  do j1=-jmax,jmax
                     do i1=-imax,imax
                        r2=(i1*w11*i1+2*i1*w12*j1+j1*w22*j1)/nsub**2
                        if (r2 .lt. rcut**2) then
                           narea=narea+1
                        endif
                     enddo
                  enddo
                  do j1=-jmax,jmax
                     do i1=-imax,imax
                        r2=(i1*w11*i1+2*i1*w12*j1+j1*w22*j1)/nsub**2
                        if (r2 .lt. rcut**2) then
                           ix=i1+nsub*(i+dx)+n
                           iy=j1+nsub*(j+dy)+n
                           ix=mod(ix+n,n)+1
                           iy=mod(iy+n,n)+1
                           xy(ix,iy,ni3)=xy(ix,iy,ni3)+tsz/narea
                           xy1(ix,iy,ni3)=xy1(ix,iy,ni3)+tksz/narea
                           xy2(ix,iy,ni3)=xy2(ix,iy,ni3)+alx/narea
                           xy3(ix,iy,ni3)=xy3(ix,iy,ni3)+ttalx/narea
                           xy4(ix,iy,ni3)=xy4(ix,iy,ni3)+tm/narea
                           xy5(ix,iy,ni3)=xy5(ix,iy,ni3)+phi2/narea/a
                        endif
                     enddo
                  enddo

               else if (mod(icount,3) .eq. 1) then
! now repeat for yz:
                  w11=s22 - s12**2/s11
                  w12=s23 - s12*s13/s11
                  w22=s33 - s13**2/s11
                  
! we apply the same procedure again to find the bounding boxes
                  
                  zmax=1/sqrt(w22-w12**2/w11)
                  
                  ymax=1/sqrt(w11-w12**2/w22)
                  jmax=ymax*nsub*rcut+1
                  kmax=zmax*nsub*rcut+1
                  narea=0
                  do k1=-kmax,kmax
                     do j1=-jmax,jmax
                        r2=(j1*w11*j1+2*j1*w12*k1+k1*w22*k1)/nsub**2
                        if (r2 .lt. rcut**2) then
                           narea=narea+1
                        endif
                     enddo
                  enddo
                  do k1=-kmax,kmax
                     do j1=-jmax,jmax
                        r2=(j1*w11*j1+2*j1*w12*k1+k1*w22*k1)/nsub**2
                        if (r2 .lt. rcut**2) then
                           iz=k1+nsub*(k+dz)+n
                           iy=j1+nsub*(j+dy)+n
                           iz=mod(iz+n,n)+1
                           iy=mod(iy+n,n)+1
                           yz(iy,iz,ni3)=yz(iy,iz,ni3)+tsz/narea
                           yz1(iy,iz,ni3)=yz1(iy,iz,ni3)+tksz/narea
                           yz2(iy,iz,ni3)=yz2(iy,iz,ni3)+alx/narea
                           yz3(iy,iz,ni3)=yz3(iy,iz,ni3)+ttalx/narea
                           yz4(iy,iz,ni3)=yz4(iy,iz,ni3)+tm/narea
                           yz5(iy,iz,ni3)=yz5(iy,iz,ni3)+phi2/narea/a
                        endif
                     enddo
                  enddo
!           if (mod(icount,3) .eq. 0)
!           else if (mod(icount,3) .eq. 1)
               else
! now repeat for xz:
                  w11=s11 - s12**2/s22
                  w12=s13 - s12*s23/s22
                  w22=s33 - s23**2/s22
                  
                  xmax=1/sqrt(w11-w12**2/w22)
                  zmax=1/sqrt(w22-w12**2/w11)
                  imax=xmax*nsub*rcut+1
                  kmax=zmax*nsub*rcut+1
                  narea=0
                  do k1=-kmax,kmax
                     do i1=-imax,imax
                        r2=(i1*w11*i1+2*i1*w12*k1+k1*w22*k1)/nsub**2
                        if (r2 .lt. rcut**2) then
                           narea=narea+1
                        endif
                     enddo
                  enddo
                  do k1=-kmax,kmax
                     do i1=-imax,imax
                        r2=(i1*w11*i1+2*i1*w12*k1+k1*w22*k1)/nsub**2
                        if (r2 .lt. rcut**2) then
                           iz=k1+nsub*(k+dz)+n
                           ix=i1+nsub*(i+dx)+n
                           iz=mod(iz+n,n)+1
                           ix=mod(ix+n,n)+1
                           xz(ix,iz,ni3)=xz(ix,iz,ni3)+tsz/narea
                           xz1(ix,iz,ni3)=xz1(ix,iz,ni3)+tksz/narea
                           xz2(ix,iz,ni3)=xz2(ix,iz,ni3)+alx/narea
                           xz3(ix,iz,ni3)=xz3(ix,iz,ni3)+ttalx/narea
                           xz4(ix,iz,ni3)=xz4(ix,iz,ni3)+tm/narea
                           xz5(ix,iz,ni3)=xz5(ix,iz,ni3)+phi2/narea/a
                        endif
                     enddo
                  enddo
!           if (mod(icount,3) .eq. 0)
!           else if (mod(icount,3) .eq. 1)
               endif
            enddo
         enddo
      enddo
!      if (mod(icount,3) .ne. 1)
      else
c we parallelize over the outer loop if we project over j:
c*$* assert do (concurrent)
c$omp parallel do default(private)
c$omp&         shared(def,dxscale,u,velfact,conx,conksz,phi,a,consz,te
c$omp&        ,icount,nsub,rcut,yz,yz1,yz2,yz3,yz4,yz5,masstot)
c$omp&    reduction(+:checksum) 
      do k=1,ng3
         do j=1,ng2
            do i=1,ng1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
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
     &              -a12**2*a33
               x=i+(def(ip,j,k)-def(im,j,k))/2
               y=j+(def(i,jp,k)-def(i,jm,k))/2
               z=k+(def(i,j,kp)-def(i,j,km))/2
               x=mod(x+ng1,1.*ng1)
               y=mod(y+ng2,1.*ng2)
               z=mod(z+ng3,1.*ng3)
               
               dd1=x*dxscale
               dd2=y*dxscale
               dd3=z*dxscale
c     Nearest grid point, different than CIC or TSC.
               i1=dd1+0.5
               i2=dd2+0.5
               i3=dd3+0.5
               if (i1.eq.0) i1=n1sz
               if (i2.eq.0) i2=n2sz
               if (i3.eq.0) i3=n3sz
               am=u(1,i,j,k)
               therm=(u(5,i,j,k)-(u(2,i,j,k)**2+u(3,i,j,k)**2+
     &             u(4,i,j,k)**2)/u(1,i,j,k)/2)*
     &              velfact**2*2
               vmean1=u(2,i,j,k)*velfact/am
               vmean2=u(3,i,j,k)*velfact/am
               vmean3=u(4,i,j,k)*velfact/am
               sig2=therm/am
               sig1=sqrt(abs(sig2))
               sig3=sig1*sig2
               rhorhobar=am/masstot*n1sz*n2sz*n3sz
               ni1=(i1-1)*nrhoim/n1sz+1
               ni2=(i2-1)*nrhoim/n2sz+1
               ni3=(i3-1)*nrhoim/n3sz+1
               alx=conx*am**2*sqrt(abs(therm/am
     &              ))/det/n1sz
               tsz=rhorhobar*sig2*consz/n1sz
c we know that we are projecting down X axis.
 	       tksz=rhorhobar*vmean1*conksz/n1sz
	       ttalx=alx*sig1*te/n1sz
               tm=am/masstot/n1sz
               phi2=phi(i,j,k)*det/n1sz
               checksum=checksum+tm

               b11=(a22*a33-a23**2)/det
               b12=(a13*a23-a12*a33)/det
               b13=(a12*a23-a13*a22)/det
               b22=(a11*a33-a13**2)/det
               b23=(a12*a13-a11*a23)/det
               b33=(a11*a22-a12**2)/det
               s11=(b11**2+b12**2+b13**2)
               s12=(b11*b12+b12*b22+b13*b23)
               s13=(b11*b13+b12*b23+b13*b33)
               s22=(b12**2+b22**2+b23**2)        
               s23=(b12*b13+b22*b23+b23*b33)
               s33=(b13**2+b23**2+b33**2)
               
               dx=(def(ip,j,k)-def(im,j,k))/2
               dy=(def(i,jp,k)-def(i,jm,k))/2
               dz=(def(i,j,kp)-def(i,j,km))/2
! we know that mod(icount,3)==1
               w11=s22 - s12**2/s11
               w12=s23 - s12*s13/s11
               w22=s33 - s13**2/s11

! we apply the same procedure again to find the bounding boxes

               zmax=1/sqrt(w22-w12**2/w11)
               
               ymax=1/sqrt(w11-w12**2/w22)
               jmax=ymax*nsub*rcut+1
               kmax=zmax*nsub*rcut+1
               narea=0
               do k1=-kmax,kmax
                  do j1=-jmax,jmax
                     r2=(j1*w11*j1+2*j1*w12*k1+k1*w22*k1)/nsub**2
                     if (r2 .lt. rcut**2) then
                        narea=narea+1
                     endif
                  enddo
               enddo
               do k1=-kmax,kmax
                  do j1=-jmax,jmax
                     r2=(j1*w11*j1+2*j1*w12*k1+k1*w22*k1)/nsub**2
                     if (r2 .lt. rcut**2) then
                        iz=k1+nsub*(k+dz)+n
                        iy=j1+nsub*(j+dy)+n
                        iz=mod(iz+n,n)+1
                        iy=mod(iy+n,n)+1
                        yz(iy,iz,ni3)=yz(iy,iz,ni3)+tsz/narea
                        yz1(iy,iz,ni3)=yz1(iy,iz,ni3)+tksz/narea
                        yz2(iy,iz,ni3)=yz2(iy,iz,ni3)+alx/narea
                        yz3(iy,iz,ni3)=yz3(iy,iz,ni3)+ttalx/narea
                        yz4(iy,iz,ni3)=yz4(iy,iz,ni3)+tm/narea
                        yz5(iy,iz,ni3)=yz5(iy,iz,ni3)+phi2/narea/a
                     endif
                  enddo
               enddo
               
            enddo
          enddo
        enddo
!      if (mod(icount,3) .ne. 1)
!      else
      endif
      checksum1=0
      checksum2=0
      checksum3=0
      ni1=1
      ni2=1
      ni3=1
c$omp parallel  do reduction(+:checksum1,checksum2,checksum3)
c$omp& shared(xy4,xz4,yz4,ni1,ni2,ni3)
      do j=1,n
         do i=1,n
            checksum1=checksum1+xy4(i,j,ni1)
            checksum2=checksum2+xz4(i,j,ni2)
            checksum3=checksum3+yz4(i,j,ni3)
         enddo
      enddo
      err=(1-checksum1/checksum)**2+(1-checksum2/checksum)**2
     &     +(1-checksum3/checksum)**2
      if (err .gt. 0.0001) then
c note that the three checksums are really the same because of
c equivalencing.
         write(*,*) 'checksum error: c,c1,c2,c3=',
     &        checksum,checksum1,checksum2,checksum3
      endif

c write output
      write(25)aout
      write(26)aout
      write(27)aout
      write(28)aout
      write(29)aout
      write(30)aout
      write(107)aout
      if (.false.) then
         call warray(xy,n1sz*n2sz,25)
         call warray(xy1,n1sz*n2sz,26)
         call warray(xy2,n1sz*n2sz,27)
         call warray(xy3,n1sz*n2sz,28)
         call warray(xy4,n1sz*n2sz,29)
      else
         if (mod(icount,3) .eq. 0) then
            write(25)(((xy(k,j,i), k=1,n2sz), j=1,n1sz), i=1,nrhoim)
            write(26)(((xy1(k,j,i), k=1,n2sz), j=1,n1sz), i=1,nrhoim)
            write(27)(((xy2(k,j,i), k=1,n2sz), j=1,n1sz), i=1,nrhoim)
            write(28)(((xy3(k,j,i), k=1,n2sz), j=1,n1sz), i=1,nrhoim)
            write(29)(((xy4(k,j,i), k=1,n2sz), j=1,n1sz), i=1,nrhoim)
         else if (mod(icount,3) .eq. 1) then
            write(25)(((yz(k,j,i), k=1,n3sz), j=1,n2sz), i=1,nrhoim)
            write(26)(((yz1(k,j,i), k=1,n3sz), j=1,n2sz), i=1,nrhoim)
            write(27)(((yz2(k,j,i), k=1,n3sz), j=1,n2sz), i=1,nrhoim)
            write(28)(((yz3(k,j,i), k=1,n3sz), j=1,n2sz), i=1,nrhoim)
            write(29)(((yz4(k,j,i), k=1,n3sz), j=1,n2sz), i=1,nrhoim)
         else
            write(25)(((xz(k,j,i), k=1,n3sz), j=1,n1sz), i=1,nrhoim)
            write(26)(((xz1(k,j,i), k=1,n3sz), j=1,n1sz), i=1,nrhoim)
            write(27)(((xz2(k,j,i), k=1,n3sz), j=1,n1sz), i=1,nrhoim)
            write(28)(((xz3(k,j,i), k=1,n3sz), j=1,n1sz), i=1,nrhoim)
            write(29)(((xz4(k,j,i), k=1,n3sz), j=1,n1sz), i=1,nrhoim)
         endif
      endif
      if (mod(icount,3) .eq. 0) then
         call wimage(xy,n2sz,n1sz,'gas')
         if (mod(icount,15) .eq. 0) then
            call wimage(xy1,n2sz,n1sz,'ksz')
            call wimage(xy2,n2sz,n1sz,'glx')
            call wimage(xy4,n2sz,n1sz,'dgs')
            call wimage(xy5,n2sz,n1sz,'phi')
         endif
      endif
      return
      end




