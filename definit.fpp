c-*- fortran -*-

c#define TEST
      subroutine definit(def,r0,compress,compress2)
c compress the grid to dx=1/compress within r<n*r0
      include 'relaxgl.fi'
      real def(ng1,ng2,ng3),r0,compress

c locals
      integer ntable
      parameter(ntable=ng3*100)
      real rhotable(0:ntable), deftable(0:2*ntable)
      real r1,r2,dr1,v1,v2,v3,term1,term2,term3,r,coef1,coef2,r3
     &            ,vcum,sol1,rint,v, rr,ir,dr
      integer i,j,k
      real*4 densinit(ng1,ng2,ng3)
cmydist densinit(*,block,*)
      common /densi/densinit
     
      r1=.8*r0+.2
      r2=.3*r0+.7
      r3=.15*r0+.85
      v1=1./compress**3
      v2=1./compress2**3
      vcum=r0**3*v1+(v1+r0*(v1-v2)/(r1-r0))* (r1**3-r0**3)
     &     + (3./4.)*(r1**4-r0**4)*(v2-v1)/(r1-r0)
     &    +v2*(r2**3-r1**3)
#ifdef TEST
      write(*,*) 'vcum=',vcum
#endif
c now we need to solve
c      (V2+r2*(V2-V3)/(r3-r2)) (r3**3-r2**3)
c     + (3./4.)*(r3**4-r2**4)*(v3-v2)/(r3-r2) + next term = 1-vcum
c
c next term= (v3+r3*(v3-1))*(1**3-r3**3)/(r3-r2)
c     + (3./4.)*(1**4-r3**4)*(1-v3)/(1-r3)
      term1= (v2+r2*v2/(r3-r2))*(r3**3-r2**3)
     &      + (3./4.)*(r3**4-r2**4)*(-v2)/(r3-r2)
c coefficient of v3-v2
      term2=(r3*(-1))*(1**3-r3**3)/(r3-r2) 
     &    + (3./4.)*(1**4-r3**4)/(1-r3)
      coef1= -(r2)*(r3**3-r2**3)/(r3-r2) 
     &         + (3./4.)*(r3**4-r2**4)*(1)/(r3-r2)
      coef2=(1+r3/(r3-r2))*(1**3-r3**3)
     &       + (3./4.)*(1**4-r3**4)*(-1)/(1-r3)
      sol1=(1-vcum-term1-term2)/(coef1+coef2)
      v3=sol1
#ifdef TEST
      write(*,*) 'term1,2=',term1,term2,term1+term2
      write(*,*) 'coef1,2=',coef1,coef2,coef1+coef2
      write(*,*) 'vcum,v3=',vcum,v3
#endif
      do i=0,ntable
         r=(i+0.5)/ntable
         if ( r .lt. r0 ) then
            rhotable(i)=v1
c cover volume fraction r0^3 v1
         else if ( r .lt. r1 ) then
            dr1=(r-r0)/(r1-r0)
            rhotable(i)=(1-dr1)*v1+(dr1)*v2
c let V1=1/compress^3, V2=1/compress2**3
c shell volume is (V1+r0(V1-V2)) (r1^3-r0^3)/(r1-r0)
c                + (3/4)(r1^4-r0^4)(V2-V1)/(r1-r0)
         else if ( r .lt. r2 ) then
            rhotable(i)=v2
         else if ( r .lt. r3 ) then
            dr1=(r-r2)/(r3-r2)
            rhotable(i)=(1-dr1)*v2+dr1*v3
         else
            dr1=(r-r3)/(1-r3)
c we want it to go smoothly onto the background flow
            rhotable(i)=(1-dr1)*v3+dr1
         endif
      enddo
c one may wish to smooth rhotable out a bit to make flow more smooth
      
      
      testsum=0
      deftable(0)=0
      dr=(ng3/2.)/ntable
      rint=0
      do i=1,ntable
         r=(i-0.5)*0.5*ng3/ntable
         rint=rint+rhotable(i-1)*dr*r**2*3
         v=r**3*rint
         
c         deftable(i)=deftable(i-1)+r*(-1+(rhotable(i-1))**(1./3.))*dr
         deftable(i)=deftable(i-1)+(rint**(1./3.)-r)*dr
         testsum=testsum+3*(r/(0.5*ng3))**2*rhotable(i-1)/ntable
#ifdef TEST
         write(10,*) r,deftable(i),rhotable(i)
#endif
      enddo
      do i=ntable+1,2*ntable
         deftable(i)=deftable(ntable)
      enddo

#ifdef TEST
      write(*,*) 'test integral=',testsum,' (should be 1)'
#endif

c zpj move "real rr,ir,dr" to the begining
c$omp parallel default(private) shared(def,deftable)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               r=sqrt((i-ng1/2-0.5)**2+(j-ng2/2-0.5)**2
     & +(k-ng3/2-0.5)**2)

               rr=r*ntable*2./ng3
               ir=rr
               dr=rr-ir
               if (ir .ge. 2*ntable) then
                  write(*,*) 'definit: ir=',ir
                  ir=ntable
               endif
               def(i,j,k)=(1-dr)*deftable(ir)+dr*deftable(ir+1)             
            enddo
         enddo
      enddo
c$omp end  parallel
      return   
      end

#ifdef TEST
      include 'relaxgl.fi'
      real def(ng1,ng2,ng3)
      call definit(def,0.5,4.,2.)

      j=ng2/2
      do i=1,ng1
         do k=1,ng3
            ip=mod(i,ng1)+1
            im=mod(i+ng1-2,ng1)+1
            kp=mod(k,ng3)+1
            km=mod(k+ng3-2,ng3)+1
            x=i+(def(ip,j,k)-def(im,j,k))/2
            z=k+(def(i,j,kp)-def(i,j,km))/2
            write(15,*) x,z,j
         enddo
      enddo
      i=ng1/2
      do k=1,ng3
         ip=mod(i,ng1)+1
         im=mod(i+ng1-2,ng1)+1
         jp=mod(j,ng2)+1
         jm=mod(j+ng2-2,ng2)+1
         kp=mod(k,ng3)+1
         km=mod(k+ng3-2,ng3)+1
         phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
         phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
         phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
c 3*( 2 A, 2 M )            
         phixy=(def(ip,jp,k)-def(im,jp,k)
     &        -def(ip,jm,k)+def(im,jm,k))/4
         phiyz=(def(i,jp,kp)-def(i,jp,km)
     &        -def(i,jm,kp)+def(i,jm,km))/4
         phixz=(def(ip,j,kp)-def(im,j,kp)
     &        -def(ip,j,km)+def(im,j,km))/4
c 3*( 3 A, 1 M )         
         a11=1+phixx
         a12=phixy
         a13=phixz
         a22=1+phiyy
         a23=phiyz
         a33=1+phizz
c det is \sqrt{g}            
         det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &           -a12**2*a33

         write(19,*) k,det,1/abs(det)**(1/3.)
      enddo
      end

      

#endif







