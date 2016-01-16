! this routine is so wasteful of memory that it shouldn't be used
subroutine xvmap(xv,xveucl,def,np,n,a)
  implicit none
  include 'globalpa.fi90'
  integer n,np
  real def(n,n,n)
  real*4 xv(6,np),xveucl(6,np)
  real a
! locals
  real*4 xpos(3,n,n,n), dx(3)
  real velfact,x,y,z
  integer i, idim,ix,iy,iz,ixp,iyp,izp,ir
  real rho000,rho001,rho010,rho011,rho100,rho101,rho110 &
       ,rho111,wx,wy,wz,w000,w001,w010,w011,w100,w101,w110 &
       ,w111,p1,rhoc,tau,alambda,pi,amu,engyu
  
  ir(i)=mod(i-1+n,n)+1
  
  pi=4*atan(1.)
! the fun part: unit conversions into SI
  alambda=boxsize/hubbleh0/n*3.0856E22
  rhoc=hubbleh0**2 * 1.88E-26
  amu=omegab*rhoc*alambda**3
  tau=sqrt(omegab*alambda**3/(omega0*amu*6*pi*6.6726E-11 ))
  engyu=amu*alambda**2/tau**2
! to get into (proper) km/sec:
  velfact=alambda/(tau*1000*a)


  xpos(1,:,:,:)=cshift(def,1)-def  ! xpos(1,1,1,1)=def(2,1,1)-def(1,1,1)
  xpos(2,:,:,:)=cshift(def,1,2)-def
  xpos(3,:,:,:)=cshift(def,1,3)-def
  xpos(1,:,:,:)=(cshift(xpos(1,:,:,:),1,2)+xpos(1,:,:,:))/2
  xpos(1,:,:,:)=(cshift(xpos(1,:,:,:),1,3)+xpos(1,:,:,:))/2
  xpos(2,:,:,:)=(cshift(xpos(2,:,:,:),1,1)+xpos(2,:,:,:))/2
  xpos(2,:,:,:)=(cshift(xpos(2,:,:,:),1,3)+xpos(2,:,:,:))/2
  xpos(3,:,:,:)=(cshift(xpos(3,:,:,:),1,1)+xpos(3,:,:,:))/2
  xpos(3,:,:,:)=(cshift(xpos(3,:,:,:),1,2)+xpos(3,:,:,:))/2

  
  do i=1,np
     x=xv(1,i)-0.5
     y=xv(2,i)-0.5
     z=xv(3,i)-0.5
! xpos(1) -> rho(1.5)
     ix=x
     iy=y
     iz=z
     wx=x-ix
     wy=y-iy
     wz=z-iz
     w111=wx*wy*wz
     w110=wx*wy*(1-wz)
     w101=wx*(1-wy)*wz
     w100=wx*(1-wy)*(1-wz)
     w011=(1-wx)*wy*wz
     w010=(1-wx)*wy*(1-wz)
     w001=(1-wx)*(1-wy)*wz
     w000=(1-wx)*(1-wy)*(1-wz)
     ixp=modulo(ix,n)+1
     iyp=modulo(iy,n)+1
     izp=modulo(iz,n)+1
     ix=modulo(ix-1,n)+1
     iy=modulo(iy-1,n)+1
     iz=modulo(iz-1,n)+1
     do idim=1,3
        rho000=xpos(idim,ix,iy,iz)
        rho001=xpos(idim,ix,iy,izp)
        rho010=xpos(idim,ix,iyp,iz)
        rho011=xpos(idim,ix,iyp,izp)
        rho100=xpos(idim,ixp,iy,iz)
        rho101=xpos(idim,ixp,iy,izp)
        rho110=xpos(idim,ixp,iyp,iz)
        rho111=xpos(idim,ixp,iyp,izp)
        p1=rho000*w000+rho001*w001+rho010*w010	&
             +rho011*w011+rho100*w100+rho101*w101	&
             +rho110*w110+rho111*w111
        dx(idim)=p1
     enddo
     xveucl(:3,i)=xv(:3,i)+dx
! velocity conversions:
! the resulting velocities are in proper km/sec
!
     xveucl(4:,i)=xv(4:,i)*velfact
  enddo
  call limitx2(xveucl,np,n)
return
end

subroutine limitx2(xv,npart,n)
  implicit none
  integer npart,n
  real*4 xv(6,npart)
! local
  integer i
  
  do i=1,npart
     if (xv(1,i) .lt. 0) xv(1,i)=xv(1,i)+n
     if (xv(2,i) .lt. 0) xv(2,i)=xv(2,i)+n
     if (xv(3,i) .lt. 0) xv(3,i)=xv(3,i)+n
     if (xv(1,i) .ge. n) xv(1,i)=xv(1,i)-n
     if (xv(2,i) .ge. n) xv(2,i)=xv(2,i)-n
     if (xv(3,i) .ge. n) xv(3,i)=xv(3,i)-n
  enddo
  return
end subroutine limitx2




subroutine xvproj(xv,def,np,n,a)
  implicit none
  include 'globalpa.fi90'
  integer n,np
  real def(n,n,n)
  real*4 xv(6,np)
  real a
! locals
  real*4 xpos(3,n,n,n), dx(3)
  real velfact,x,y,z
  integer i, idim,ix,iy,iz,ixp,iyp,izp,ir
  real rho000,rho001,rho010,rho011,rho100,rho101,rho110 &
       ,rho111,wx,wy,wz,w000,w001,w010,w011,w100,w101,w110 &
       ,w111,p1,rhoc,tau,alambda,pi,amu,engyu
  
  ir(i)=mod(i-1+n,n)+1
  
  pi=4*atan(1.)
! the fun part: unit conversions into SI
  alambda=boxsize/hubbleh0/n*3.0856E22
  rhoc=hubbleh0**2 * 1.88E-26
  amu=omegab*rhoc*alambda**3
  tau=sqrt(omegab*alambda**3/(omega0*amu*6*pi*6.6726E-11 ))
  engyu=amu*alambda**2/tau**2
! to get into (proper) km/sec:
  velfact=alambda/(tau*1000*a)


  xpos(1,:,:,:)=cshift(def,1)-def  ! xpos(1,1,1,1)=def(2,1,1)-def(1,1,1)
  xpos(2,:,:,:)=cshift(def,1,2)-def
  xpos(3,:,:,:)=cshift(def,1,3)-def
  xpos(1,:,:,:)=(cshift(xpos(1,:,:,:),1,2)+xpos(1,:,:,:))/2
  xpos(1,:,:,:)=(cshift(xpos(1,:,:,:),1,3)+xpos(1,:,:,:))/2
  xpos(2,:,:,:)=(cshift(xpos(2,:,:,:),1,1)+xpos(2,:,:,:))/2
  xpos(2,:,:,:)=(cshift(xpos(2,:,:,:),1,3)+xpos(2,:,:,:))/2
  xpos(3,:,:,:)=(cshift(xpos(3,:,:,:),1,1)+xpos(3,:,:,:))/2
  xpos(3,:,:,:)=(cshift(xpos(3,:,:,:),1,2)+xpos(3,:,:,:))/2

  
  do i=1,np
     x=xv(1,i)-0.5
     y=xv(2,i)-0.5
     z=xv(3,i)-0.5
! xpos(1) -> rho(1.5)
     ix=x
     iy=y
     iz=z
     wx=x-ix
     wy=y-iy
     wz=z-iz
     w111=wx*wy*wz
     w110=wx*wy*(1-wz)
     w101=wx*(1-wy)*wz
     w100=wx*(1-wy)*(1-wz)
     w011=(1-wx)*wy*wz
     w010=(1-wx)*wy*(1-wz)
     w001=(1-wx)*(1-wy)*wz
     w000=(1-wx)*(1-wy)*(1-wz)
     ixp=modulo(ix,n)+1
     iyp=modulo(iy,n)+1
     izp=modulo(iz,n)+1
     ix=modulo(ix-1,n)+1
     iy=modulo(iy-1,n)+1
     iz=modulo(iz-1,n)+1
     do idim=1,3
        rho000=xpos(idim,ix,iy,iz)
        rho001=xpos(idim,ix,iy,izp)
        rho010=xpos(idim,ix,iyp,iz)
        rho011=xpos(idim,ix,iyp,izp)
        rho100=xpos(idim,ixp,iy,iz)
        rho101=xpos(idim,ixp,iy,izp)
        rho110=xpos(idim,ixp,iyp,iz)
        rho111=xpos(idim,ixp,iyp,izp)
        p1=rho000*w000+rho001*w001+rho010*w010	&
             +rho011*w011+rho100*w100+rho101*w101	&
             +rho110*w110+rho111*w111
        dx(idim)=p1
     enddo
!     =xv(:3,i)+dx
  enddo
return
end
