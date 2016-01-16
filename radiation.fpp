c-*- Fortran -*-
      subroutine bremsstrahlung(u,def,a,t,tau,dt)
      implicit none
      include 'relaxgl.fi'
      include 'globalpa.fi'
      real u(5,ng1,ng2,ng3), def(ng1,ng2,ng3), dt,a,t,tau

c locals
      real rho,rhocgs,engy,engyunit,elx,temp,api,alambda,dx,rhocrit
     &     ,amass,umean,tempunit, tunit, edot, rne, gff, yfrac, acount
     &     ,rnp,rnhe
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &      , a23, a33, det, b11,b12,b13,b22,b23,b33
      integer i,j,k,ip,im,jp,jm,kp,km
c zpj: add definition of u1
      real u1



      alambda=boxsize/hubbleh0*3.0856E24/ng1
      api=4*atan(1.)
c in cgs:
      dx=a*alambda
      rhocrit=1.88E-29*hubbleh0**2  ! 3H^2/8\pi G
c mass unit for baryons: does not change with time    
      amass=rhocrit*alambda**3*omegab
c time unit: does not change, but the code uses \tau units which
c do change. for example, for energy, we take grid units and divide by a^2.
      tunit=sqrt(omegab*alambda**3/(omega0*amass*6*api*6.6726E-8))
c dx already contains two powers of a
      engyunit=amass*dx**2/tunit**2/a**4
c helium fraction by mass:      
      yfrac=0.24
c to calculate mean particle mass, fully ionized, is two particles per
c hydrogen, and three per helium.
c n = (1-y)*2+y*3/4
c mean molecular weight is
      acount=(1-yfrac)*2+yfrac*.75
      umean=1.67E-24/acount
c note that (1-y+y/2)/n=(1-y/2)/n of the particles are electrons,
c (1-y)/n are protons, and y/(4n) are helium
      tempunit=engyunit*(umean/amass)*(2./3.)/1.38e-16
      
c$omp parallel default(private) 
c$omp& firstprivate(ng1,ng2,ng3,yfrac,dx,amass,umean,acount,tempunit
c$omp&     ,tunit,engyunit)
c$omp& shared(def,u,a,dt)
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
     &              -a12**2*a33
               rho=u(1,i,j,k)/det
               engy=(u(5,i,j,k)-(u(2,i,j,k)**2+u(3,i,j,k)**2
     &              +u(4,i,j,k)**2)/2/u(1,i,j,k))/det
               rhocgs=rho*amass/dx**3
               rne=rhocgs/umean*(1-yfrac/2)/acount
               rnp=rhocgs/umean*(1-yfrac)/acount
               rnhe=rhocgs/umean*(yfrac*.25)/acount
c rne is approximately the number density  of electrons  (cm^-3)
c for now assume rni=rne, equal number of ions and electrons.               
               temp=engy/rho*tempunit
c gaunt factor.               
               gff=1.2
c the emission rate per cell in units of erg/sec:               
               elx=gff*rne*(rnp+4*rnhe)*sqrt(temp)*dx**3*det*1.426E-27
c convert to grid units               
               edot=elx*tunit*a**2/engyunit
c second order runge kutta:
               u1=engy-edot*dt/2/det
               temp=u1/rho*tempunit
               elx=gff*rne*(rnp+4*rnhe)*sqrt(temp)*dx**3*det*1.426E-27
               edot=elx*tunit*a**2/engyunit
               if (engy*det-edot*dt .gt. 0) then
                  u(5,i,j,k)=u(5,i,j,k)-edot*dt
               endif
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end











