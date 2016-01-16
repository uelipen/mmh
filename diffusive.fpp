c-*- Fortran -*-
c written August 21, 1997 by Ue-Li Pen, upen@cfa.harvard.edu
c to implement diffusive effects in cluster simulations
c
      subroutine diffuse(u,def,a,t,tau,dt)
      implicit none
      include 'relaxgl.fi'
      include 'globalpa.fi'
      real u(5,ng1,ng2,ng3), def(ng1,ng2,ng3), dt,a,t,tau

c locals
      real rho,rhocgs,engy,engyunit,temp,api,alambda,dx,rhocrit
     &     ,amass,umean,tempunit, tunit, edot, rne, gff, yfrac, acount
     &     ,rnp,rnhe
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &      , a23, a33, det, b11,b12,b13,b22,b23,b33
c zpj: add the definition of private variables in parallel
      real tempmx,temppyz,temppymx,temppymz,temppzmx,gmxp1,gmyp1,gmzp1,fmxp1
     &     ,fmyp1,fmzp1,fep1,gep1,temppzmy,hmxp1,hmup1,hmzp1,hep1,tk,gep,hep
     &     ,fephmyp,fmzp,gmzp,hmzp,vxph,vyph,vzph,fep,fmxp1,hmyp1,vzdzp1
     &     ,divvp1,tempcgs,mup,fmxp,gmxp,hmxp,fmyp,gmyp,vydxp1,vydyp1
     &     ,vydzp1,vzdzp,b11p,b12p,b13p,b22p,b23p,b33p,detp,tempdxp1
     &     ,vxdxp,vxdyp,vxdzp,vydxp,vydyp,tempdyp1,tempdzp1,vxdxp1
     &     ,vxdyp1,vxdzp1,vydzp,vzdxp,vzdyp,b11p,tempmz,vzdxp1,vzdyp1
     &     ,temppxy,temppxz,temppxmy,temppxmz,tempdyp,tempdzp,vxp,vyp
     &     ,vydyp,vzp,rgu,vx,vy,vz,rgup,temppx,tempdxp,temppy,temppz
     &     ,tempmy,tempmz,tempdzp,alpha,eta,hmyp,hep,dept
      real gbaj(7,ng1,ng2,ng3),fmps1(4,ng1,ng2,2),gmps1(4,ng1,ng2)
     &     ,hmps1(4,ng1,ng2,2), uflux(4,ng1,ng2)   
      integer ii,k1,ki
c zpj: end

      integer i,j,k,ip,im,jp,jm,kp,km
      include 'gmetric.fi'

#ifndef GMETRIC
      write(*,*) 'viscosity only implemented for global metric'
      stop
#endif


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
      
      alpha=0.3

      eta=5500
      
c$omp parallel default(private) shared(gbaj,fmps1,gmps1,hmps1,uflux)
c$omp& firstprivate(k,amass,dx,tempunit,eta,tunit)
      do k1=0,ng3+1
         ki=mod(k1,2)
c j loop should be fully parallel

c$omp do 
        do j=1,ng2
            do i=1,ng1
               k=mod(k1,ng3)+1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               b11=gbaj(1,i,j,k)
               b12=gbaj(2,i,j,k)
               b13=gbaj(3,i,j,k)
               b22=gbaj(4,i,j,k)
               b23=gbaj(5,i,j,k)
               b33=gbaj(6,i,j,k)
               det=gbaj(7,i,j,k)
               rgu=u(1,i,j,k)
               rho=rgu/det
               engy=(u(5,i,j,k)-(u(2,i,j,k)**2+u(3,i,j,k)**2
     &              +u(4,i,j,k)**2)/2/u(1,i,j,k))/det
               rhocgs=rho*amass/dx**3
               temp=engy/rho
               vx=u(2,i,j,k)/rgu
               vy=u(3,i,j,k)/rgu
               vz=u(4,i,j,k)/rgu



c compute all velocity gradients at (i+1/2)
               rgup=u(1,ip,j,k)
               temppx=(u(5,ip,j,k)-(u(2,ip,j,k)**2+u(3,ip,j,k)**2
     &              +u(4,ip,j,k)**2)/2/u(1,ip,j,k))/rgup
               tempdxp=temppx-temp
               temppy=(u(5,i,jp,k)-(u(2,i,jp,k)**2+u(3,i,jp,k)**2
     &              +u(4,i,jp,k)**2)/2/u(1,i,jp,k))/u(1,i,jp,k)
               temppz=(u(5,i,j,kp)-(u(2,i,j,kp)**2+u(3,i,j,kp)**2
     &              +u(4,i,j,kp)**2)/2/u(1,i,j,kp))/u(1,i,j,kp)
               tempmy=(u(5,i,jm,k)-(u(2,i,jm,k)**2+u(3,i,jm,k)**2
     &              +u(4,i,jm,k)**2)/2/u(1,i,jm,k))/u(1,i,jm,k)
               tempmz=(u(5,i,j,km)-(u(2,i,j,km)**2+u(3,i,j,km)**2
     &              +u(4,i,j,km)**2)/2/u(1,i,j,km))/u(1,i,j,km)
               temppxy=(u(5,ip,jp,k)-(u(2,ip,jp,k)**2+u(3,ip,jp,k)**2
     &              +u(4,ip,jp,k)**2)/2/u(1,ip,jp,k))/u(1,ip,jp,k)
               temppxz=(u(5,ip,j,kp)-(u(2,ip,j,kp)**2+u(3,ip,j,kp)**2
     &              +u(4,ip,j,kp)**2)/2/u(1,ip,j,kp))/u(1,ip,j,kp)             
               temppxmy=(u(5,ip,jm,k)-(u(2,ip,jm,k)**2+u(3,ip,jm,k)**2
     &              +u(4,ip,jm,k)**2)/2/u(1,ip,jm,k))/u(1,ip,jm,k)
               temppxmz=(u(5,ip,j,km)-(u(2,ip,j,km)**2+u(3,ip,j,km)**2
     &              +u(4,ip,j,km)**2)/2/u(1,ip,j,km))/u(1,ip,j,km)             

               tempdyp=(temppy-tempmy+temppxy-temppxmy)/4
               tempdzp=(temppz-tempmz+temppxz-temppxmz)/4
               vxp=u(2,ip,j,k)/rgup
               vyp=u(3,ip,j,k)/rgup
               vzp=u(4,ip,j,k)/rgup
               vxdxp=vxp-vx
               vxdyp=(u(2,ip,jp,k)-u(2,ip,jm,k)+u(2,i,jp,k)-u(2,i,jm,k))/4
               vxdzp=(u(2,ip,j,kp)-u(2,ip,j,km)+u(2,i,j,kp)-u(2,i,j,km))/4
               vydxp=vyp-vy
               vydyp=(u(3,ip,jp,k)-u(3,ip,jm,k)+u(3,i,jp,k)-u(3,i,jm,k))/4
               vydzp=(u(3,ip,j,kp)-u(3,ip,j,km)+u(3,i,j,kp)-u(3,i,j,km))/4
               vzdxp=vzp-vz
               vzdyp=(u(4,ip,jp,k)-u(4,ip,jm,k)+u(4,i,jp,k)-u(4,i,jm,k))/4
               vzdzp=(u(4,ip,j,kp)-u(4,ip,j,km)+u(4,i,j,kp)-u(4,i,j,km))/4
               b11p=(b11+gbaj(1,ip,j,k))/2
               b12p=(b12+gbaj(2,ip,j,k))/2
               b13p=(b13+gbaj(3,ip,j,k))/2
               b22p=(b22+gbaj(4,ip,j,k))/2
               b23p=(b23+gbaj(5,ip,j,k))/2
               b33p=(b33+gbaj(6,ip,j,k))/2
               detp=(det+gbaj(7,ip,j,k))/2
               
               tempdxp1=tempdxp*b11p+tempdyp*b12p+tempdzp*b13p
               tempdyp1=tempdxp*b12p+tempdyp*b22p+tempdzp*b23p
               tempdzp1=tempdxp*b13p+tempdyp*b23p+tempdzp*b33p
               vxdxp1=vxdxp*b11p+vxdyp*b12p+vxdzp*b13p
               vxdyp1=vxdxp*b12p+vxdyp*b22p+vxdzp*b23p
               vxdzp1=vxdxp*b13p+vxdyp*b23p+vxdzp*b23p
               vydxp1=vydxp*b11p+vydyp*b12p+vydzp*b13p
               vydyp1=vydxp*b12p+vydyp*b22p+vydzp*b23p
               vydzp1=vydxp*b13p+vydyp*b23p+vydzp*b23p
               vzdxp1=vzdxp*b11p+vzdyp*b12p+vzdzp*b13p
               vzdyp1=vzdxp*b12p+vzdyp*b22p+vzdzp*b23p
               vzdzp1=vzdxp*b13p+vzdyp*b23p+vzdzp*b23p

               divvp1=vxdxp1+vydyp1+vzdzp1

               tempcgs=(temp+temppx)/2*tempunit
c assumes a ln\Lambda=40, from Sarazin, eqn 5.45
c zpj: I change ^ to ** before 2.5
               mup=eta*(tempcgs/1.e8)**2.5/(amass*dx/tunit**2)

c and all fluxes
               fmxp=2*mup*vxdxp1-(2./3.)*divvp1
               gmxp=mup*(vxdyp1+vydxp1)
               hmxp=mup*(vxdzp1+vzdxp1)
               fmyp=mup*(vxdyp1+vydxp1)
               gmyp=2*mup*vydyp1-(2./3.)*divvp1
               hmyp=mup*(vzdxp1+vxdzp1)
               fmzp=mup*(vxdzp1+vzdxp1)
               gmzp=mup*(vydzp1+vzdyp1)
               hmzp=2*mup*vzdzp1-(2./3.)*divvp1
c and the energy flux
               vxph=(vx+vxp)/2
               vyph=(vy+vyp)/2
               vzph=(vz+vzp)/2
               fep=2*mup*vxdxp1*vxph+mup*(vydxp1+vxdyp1)*vyph
     &              +mup*(vzdxp1+vxdzp1)*vzph-(2./3.)*mup*divvp1*vxph
     &              +tk*tempdxp1
               gep=2*mup*vydyp1*vyph+mup*(vydxp1+vxdyp1)*vxph
     &              +mup*(vzdyp1+vydzp1)*vzph-(2./3.)*mup*divvp1*vyph
     &              +tk*tempdyp1
               hep=2*mup*vzdzp1*vyph+mup*(vzdxp1+vxdzp1)*vxph
     &              +mup*(vydzp1+vzdyp1)*vyph-(2./3.)*mup*divvp1*vzph
     &              +tk*tempdzp1
c now rotate all fluxes into the curvilinear frame
               fmxp1=(fmxp*b11p+gmxp*b12p+hmxp*b13p)*detp
               fmyp1=(fmyp*b11p+gmyp*b12p+hmyp*b13p)*detp
               fmzp1=(fmzp*b11p+gmzp*b12p+hmzp*b13p)*detp
               fep1=(fep*b11p+gep*b12p+hep*b13p)*detp
               fmps1(1,i,j,ki)=fmxp1
               fmps1(2,i,j,ki)=fmyp1
               fmps1(3,i,j,ki)=fmzp1
               fmps1(4,i,j,ki)=fep1

c now do the whole thing again for the y flux:
c compute all velocity gradients at (j+1/2)
c all trailing 'p's now indicate the j offset
               rgup=u(1,i,jp,k)
               tempdyp=temppy-temp
               tempmx=(u(5,im,j,k)-(u(2,im,j,k)**2+u(3,im,j,k)**2
     &              +u(4,im,j,k)**2)/2/u(1,im,j,k))/u(1,im,j,k)
               temppyz=(u(5,i,jp,kp)-(u(2,i,jp,kp)**2+u(3,i,jp,kp)**2
     &              +u(4,i,jp,kp)**2)/2/u(1,i,jp,kp))/u(1,i,jp,kp)
               temppymx=(u(5,im,jp,k)-(u(2,im,jp,k)**2+u(3,im,jp,k)**2
     &              +u(4,im,jp,k)**2)/2/u(1,im,jp,k))/u(1,im,jp,k)
               temppymz=(u(5,i,jp,km)-(u(2,i,jp,km)**2+u(3,i,jp,km)**2
     &              +u(4,i,jp,km)**2)/2/u(1,i,jp,km))/u(1,i,jp,km)             

               tempdxp=(temppx-tempmx+temppxy-temppymx)/4
               tempdzp=(temppz-tempmz+temppyz-temppymz)/4
               vxp=u(2,i,jp,k)/rgup
               vyp=u(3,i,jp,k)/rgup
               vzp=u(4,i,jp,k)/rgup
               vxdyp=vxp-vx
               vxdxp=(u(2,ip,jp,k)-u(2,im,jp,k)+u(2,ip,j,k)-u(2,im,j,k))/4
               vxdzp=(u(2,i,jp,kp)-u(2,i,jp,km)+u(2,i,j,kp)-u(2,i,j,km))/4
               vydyp=vyp-vy
               vydxp=(u(3,ip,jp,k)-u(3,im,jp,k)+u(3,ip,j,k)-u(3,im,j,k))/4
               vydzp=(u(3,i,jp,kp)-u(3,i,jp,km)+u(3,i,j,kp)-u(3,i,j,km))/4
               vzdyp=vzp-vz
               vzdxp=(u(4,ip,jp,k)-u(4,im,jp,k)+u(4,ip,j,k)-u(4,im,j,k))/4
               vzdzp=(u(4,i,jp,kp)-u(4,i,jp,km)+u(4,i,j,kp)-u(4,i,j,km))/4
               b11p=(b11+gbaj(1,i,jp,k))/2
               b12p=(b12+gbaj(2,i,jp,k))/2
               b13p=(b13+gbaj(3,i,jp,k))/2
               b22p=(b22+gbaj(4,i,jp,k))/2
               b23p=(b23+gbaj(5,i,jp,k))/2
               b33p=(b33+gbaj(6,i,jp,k))/2
               detp=(det+gbaj(7,i,jp,k))/2
               
               tempdxp1=tempdxp*b11p+tempdyp*b12p+tempdzp*b13p
               tempdyp1=tempdxp*b12p+tempdyp*b22p+tempdzp*b23p
               tempdzp1=tempdxp*b13p+tempdyp*b23p+tempdzp*b33p
               vxdxp1=vxdxp*b11p+vxdyp*b12p+vxdzp*b13p
               vxdyp1=vxdxp*b12p+vxdyp*b22p+vxdzp*b23p
               vxdzp1=vxdxp*b13p+vxdyp*b23p+vxdzp*b23p
               vydxp1=vydxp*b11p+vydyp*b12p+vydzp*b13p
               vydyp1=vydxp*b12p+vydyp*b22p+vydzp*b23p
               vydzp1=vydxp*b13p+vydyp*b23p+vydzp*b23p
               vzdxp1=vzdxp*b11p+vzdyp*b12p+vzdzp*b13p
               vzdyp1=vzdxp*b12p+vzdyp*b22p+vzdzp*b23p
               vzdzp1=vzdxp*b13p+vzdyp*b23p+vzdzp*b23p

               divvp1=vxdxp1+vydyp1+vzdzp1

               tempcgs=(temp+temppy)/2*tempunit
c assumes a ln\Lambda=40, from Sarazin, eqn 5.45
c zpj: I change ^ to ** before 2.5
               mup=eta*(tempcgs/1.e8)**2.5/(amass*dx/tunit**2)


c and all fluxes
               fmxp=2*mup*vxdxp1-(2./3.)*divvp1
               gmxp=mup*(vxdyp1+vydxp1)
               hmxp=mup*(vxdzp1+vzdxp1)
               fmyp=mup*(vxdyp1+vydxp1)
               gmyp=2*mup*vydyp1-(2./3.)*divvp1
               hmyp=mup*(vzdxp1+vxdzp1)
               fmzp=mup*(vxdzp1+vzdxp1)
               gmzp=mup*(vydzp1+vzdyp1)
               hmzp=2*mup*vzdzp1-(2./3.)*divvp1
c and the energy flux
               vxph=(vx+vxp)/2
               vyph=(vy+vyp)/2
               vzph=(vz+vzp)/2
               fep=2*mup*vxdxp1*vxph+mup*(vydxp1+vxdyp1)*vyph
     &              +mup*(vzdxp1+vxdzp1)*vzph-(2./3.)*mup*divvp1*vxph
     &              +tk*tempdxp1
               gep=2*mup*vydyp1*vyph+mup*(vydxp1+vxdyp1)*vxph
     &              +mup*(vzdyp1+vydzp1)*vzph-(2./3.)*mup*divvp1*vyph
     &              +tk*tempdyp1
               hep=2*mup*vzdzp1*vyph+mup*(vzdxp1+vxdzp1)*vxph
     &              +mup*(vydzp1+vzdyp1)*vyph-(2./3.)*mup*divvp1*vzph
     &              +tk*tempdzp1
c now rotate all fluxes into the curvilinear frame
               gmxp1=(fmxp*b12p+gmxp*b22p+hmxp*b23p)*detp
               gmyp1=(fmyp*b12p+gmyp*b22p+hmyp*b23p)*detp
               gmzp1=(fmzp*b12p+gmzp*b22p+hmzp*b23p)*detp
               gep1=(fep*b12p+gep*b22p+hep*b23p)*detp

c to keep the j-direction parallelizable, we will store the whole tier,
c and difference at the end.
               gmps1(1,i,j)=gmxp1
               gmps1(2,i,j)=gmyp1
               gmps1(3,i,j)=gmzp1
               gmps1(4,i,j)=gep1


c and finally for the z flux:
c compute all velocity gradients at (k+1/2)
c all trailing 'p's now indicate the k offset
               rgup=u(1,i,j,kp)
               tempdzp=temppz-temp
               temppzmx=(u(5,im,j,kp)-(u(2,im,j,kp)**2+u(3,im,j,kp)**2
     &              +u(4,im,j,kp)**2)/2/u(1,im,j,kp))/u(1,im,j,kp)             
               temppzmy=(u(5,i,jm,kp)-(u(2,i,jm,kp)**2+u(3,i,jm,kp)**2
     &              +u(4,i,jm,kp)**2)/2/u(1,i,jm,kp))/u(1,i,jm,kp)             

               tempdxp=(temppx-tempmx+temppxz-temppzmx)/4
               tempdyp=(temppy-tempmy+temppyz-temppzmy)/4
               vxp=u(2,i,j,kp)/rgup
               vyp=u(3,i,j,kp)/rgup
               vzp=u(4,i,j,kp)/rgup
               vxdzp=vxp-vx
               vxdxp=(u(2,ip,j,kp)-u(2,im,j,kp)+u(2,ip,j,k)-u(2,im,j,k))/4
               vxdyp=(u(2,i,jp,kp)-u(2,i,jm,kp)+u(2,i,jp,k)-u(2,i,jm,k))/4
               vydzp=vyp-vy
               vydxp=(u(3,ip,j,kp)-u(3,im,j,kp)+u(3,ip,j,k)-u(3,im,j,k))/4
               vydyp=(u(3,i,jp,kp)-u(3,i,jm,kp)+u(3,i,jp,k)-u(3,i,j,km))/4
               vzdzp=vzp-vz
               vzdxp=(u(4,ip,j,kp)-u(4,im,j,kp)+u(4,ip,j,k)-u(4,im,j,k))/4
               vzdyp=(u(4,i,jp,kp)-u(4,i,jm,kp)+u(4,i,jp,k)-u(4,i,j,km))/4
               b11p=(b11+gbaj(1,i,j,kp))/2
               b12p=(b12+gbaj(2,i,j,kp))/2
               b13p=(b13+gbaj(3,i,j,kp))/2
               b22p=(b22+gbaj(4,i,j,kp))/2
               b23p=(b23+gbaj(5,i,j,kp))/2
               b33p=(b33+gbaj(6,i,j,kp))/2
               detp=(det+gbaj(7,i,j,kp))/2
               
               tempdxp1=tempdxp*b11p+tempdyp*b12p+tempdzp*b13p
               tempdyp1=tempdxp*b12p+tempdyp*b22p+tempdzp*b23p
               tempdzp1=tempdxp*b13p+tempdyp*b23p+tempdzp*b33p
               vxdxp1=vxdxp*b11p+vxdyp*b12p+vxdzp*b13p
               vxdyp1=vxdxp*b12p+vxdyp*b22p+vxdzp*b23p
               vxdzp1=vxdxp*b13p+vxdyp*b23p+vxdzp*b23p
               vydxp1=vydxp*b11p+vydyp*b12p+vydzp*b13p
               vydyp1=vydxp*b12p+vydyp*b22p+vydzp*b23p
               vydzp1=vydxp*b13p+vydyp*b23p+vydzp*b23p
               vzdxp1=vzdxp*b11p+vzdyp*b12p+vzdzp*b13p
               vzdyp1=vzdxp*b12p+vzdyp*b22p+vzdzp*b23p
               vzdzp1=vzdxp*b13p+vzdyp*b23p+vzdzp*b23p

               divvp1=vxdxp1+vydyp1+vzdzp1

               tempcgs=(temp+temppz)/2*tempunit
c zpj: I change ^ to ** in the next two lines and add a * after 4.6e13
c assumes a ln\Lambda=40, from Sarazin, eqn 5.45
               mup=eta*(tempcgs/1.e8)**2.5/(amass*dx/tunit**2)
c and heat conduction from (5.37)
               tk=4.6e13*(tempcgs/1.e8)**2.5/(amass*dx**2/tunit**3)

c and all fluxes
               fmxp=2*mup*vxdxp1-(2./3.)*divvp1
               gmxp=mup*(vxdyp1+vydxp1)
               hmxp=mup*(vxdzp1+vzdxp1)
               fmyp=mup*(vxdyp1+vydxp1)
               gmyp=2*mup*vydyp1-(2./3.)*divvp1
               hmyp=mup*(vzdxp1+vxdzp1)
               fmzp=mup*(vxdzp1+vzdxp1)
               gmzp=mup*(vydzp1+vzdyp1)
               hmzp=2*mup*vzdzp1-(2./3.)*divvp1
c and the energy flux
               vxph=(vx+vxp)/2
               vyph=(vy+vyp)/2
               vzph=(vz+vzp)/2
               fep=2*mup*vxdxp1*vxph+mup*(vydxp1+vxdyp1)*vyph
     &              +mup*(vzdxp1+vxdzp1)*vzph-(2./3.)*mup*divvp1*vxph
     &              +tk*tempdxp1
               gep=2*mup*vydyp1*vyph+mup*(vydxp1+vxdyp1)*vxph
     &              +mup*(vzdyp1+vydzp1)*vzph-(2./3.)*mup*divvp1*vyph
     &              +tk*tempdyp1
               hep=2*mup*vzdzp1*vyph+mup*(vzdxp1+vxdzp1)*vxph
     &              +mup*(vydzp1+vzdyp1)*vyph-(2./3.)*mup*divvp1*vzph
     &              +tk*tempdzp1
c now rotate all fluxes into the curvilinear frame
               hmxp1=(fmxp*b13p+gmxp*b23p+hmxp*b33p)*detp
               hmyp1=(fmyp*b13p+gmyp*b23p+hmyp*b33p)*detp
               hmzp1=(fmzp*b13p+gmzp*b23p+hmzp*b33p)*detp
               hep1=(fep*b13p+gep*b23p+hep*b33p)*detp

c to keep the j-direction parallelizable, we will store the whole tier,
c and difference at the end.
               hmps1(1,i,j,ki)=hmxp1
               hmps1(2,i,j,ki)=hmyp1
               hmps1(3,i,j,ki)=hmzp1
               hmps1(4,i,j,ki)=hep1
            enddo
         enddo
c update the fluxes at the very end, two tiers back:
c$omp do
         do j=1,ng2
            do i=1,ng1
               im=mod(i+ng1-2,ng1)+1
               jm=mod(j+ng2-2,ng2)+1
               km=mod(k+ng3-2,ng3)+1
               if (k .ge. 2) then
                  u(2,i,j,km)=uflux(1,i,j)
                  u(3,i,j,km)=uflux(2,i,j)
                  u(4,i,j,km)=uflux(3,i,j)
                  u(5,i,j,km)=uflux(4,i,j)
               endif
               if (k .ge. 1) then
                  do ii=1,4
                     uflux(ii,i,j)=fmps1(ii,i,j)-fmps1(ii,im,j)
     &+gmps1(1,ii,j)-gmps1(1,ii,jm)+hmps1(1,ii,j,ki)-hmps1(1,ii,j,1-ki)

                  enddo
               endif
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end




