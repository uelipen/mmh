      subroutine projphidot(phi,def,a,dt)
c this routine is called one time step after the general
c projection.  Take the time difference wrt to the previous one.
!      implicit none
      include 'globalpa.fi'
      include 'relaxgl.fi'
      include 'proj.fi'
      real phi(ng1,ng2,ng3),phidot(ng1,ng2,ng3),def(ng1,ng2,ng3)
      real a

! locals
      integer n
      parameter (n=ndimproj)
      real*4 dx,omegam,omegav,h0,epsilon,aout
      integer i,j,k,n1sz,n2sz,n3sz,i1,i2,i3,ni1,ni2,ni3
      parameter (n1sz=n,n2sz=n,n3sz=n)
      real*4 phiproj(n1sz,n2sz)
      common /phidot/ phiproj,rcut
      real me,conksz,consz,yhe,ne,te,conx,dxscale,sig1,sig2,sig3
      real rhorhobar,xpart,vmean1,vmean2,vmean3,v2mean1,v2mean2
      real v2mean3,masstot,dd1,dd2,dd3,am
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
      if (nsub*ng1 .ne. n) then
         write(*,*) 'sz_gas: n not muliple of ng1:',n,ng1
         stop
      endif
      icount=icount+1
      dx=boxsize/hubbleh0/n
      if (firsttime) then
         firsttime=.false.
      endif
      
c in meters:     
      alambda=boxsize/hubbleh0/ng1*3.0856E22

c our grid phi satisfies \Nabla^2 \phi = (2/3) \delta
c while the proper physical phi_p satisfies
c \Nabla_p^2 \phi_p = 4 \pi G \bar{\rho} \delta
c or in comoving coord \Nabla_c^2 \phi_p = 4\pi G \rho_c^0 \Omega_0 \delta/a
c so \phi_c = \delta_x^2 6 \pi G \rho_c^0 \Omega_0 \phi
c and remember \rho_c = 3 H_0^2/4 \pi G
      conphi=(boxsize/ng1)**2*9*omega0/(4*3000)

      masstot=ng1*ng2*ng3
      dxscale=n*1./ng1
      nrhoim=1
      prms=0
      phitot=0
c$omp parallel default(private) 
c$omp& shared(def,conphi,phi,a,icount,rcut,phiproj,nsub,nrhoim,dxscale
c$omp& ,dt)
c$omp& reduction(+:prms,phitot)
      do k=1,ng3
c parallelization relies on the fact that the tiling is
c very local:  cells only overlap with their nearest neighbors,
c so on any given tier, it is safe to do physically separated ones
c simultaneously.
c
c*$* assert do (concurrent)
c$omp do
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
               prms=prms+(conphi*phi(i,j,k)/a)**2*det
               phitot=phitot+phi(i,j,k)*det
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
               ni1=(i1-1)*nrhoim/n1sz+1
               ni2=(i2-1)*nrhoim/n2sz+1
               ni3=(i3-1)*nrhoim/n3sz+1
               phi2=phi(i,j,k)*det/n1sz

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
                  ix=i1+nsub*(i+dx)
                  iy=j1+nsub*(j+dy)
                  ix=mod(ix+n,n)+1
                  iy=mod(iy+n,n)+1
                  phiproj(ix,iy)=phiproj(ix,iy)-phi2/narea/a
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
               if (r2 .lt. 1) then
                  iz=k1+nsub*(k+dz)
                  iy=j1+nsub*(j+dy)
                  iz=mod(iz+n,n)+1
                  iy=mod(iy+n,n)+1
                  phiproj(iy,iz)=phiproj(iy,iz)-phi2/narea/a
               endif
            enddo
         enddo
         else
! now repeat for xz:
         w11=s11 - s12**2/s22
         w12=s13 - s12*s23/s22
         w22=s33 - s23**2/s22

         xmax=1/sqrt(w11-w12**2/w22)
         zmax=1/sqrt(w22-w12**2/w11)
         imax=xmax*nsub+1
         kmax=zmax*nsub+1
         narea=0
         do k1=-kmax,kmax
            do i1=-imax,imax
               r2=(i1*w11*i1+2*i1*w12*k1+k1*w22*k1)/nsub**2
               if (r2 .lt. 1) then
                  narea=narea+1
               endif
            enddo
         enddo
         do k1=-kmax,kmax
            do i1=-imax,imax
               r2=(i1*w11*i1+2*i1*w12*k1+k1*w22*k1)/nsub**2
               if (r2 .lt. 1) then
                  iz=k1+nsub*(k+dz)
                  ix=i1+nsub*(i+dx)
                  iz=mod(iz+n,n)+1
                  ix=mod(ix+n,n)+1
               phiproj(ix,iz)=phiproj(ix,iz)-phi2/narea/a
               endif
            enddo
         enddo

         endif

            enddo
          enddo
        enddo
c$omp do
      do j=1,n
        do i=1,n
           phiproj(i,j)=conphi*phiproj(i,j)/dt
        enddo
      enddo
c$omp end parallel
      write(*,*) 'projphi: rms phi,tot=',sqrt(prms/ng1/ng2/ng3)
     &	,phitot/ng1/ng2/ng3
c write output
        write(30) phiproj
c        call warray(phiproj,n**2,30)
      return
      end
      





