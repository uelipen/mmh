! note: compile with -r8
! generate power spectrum and correlation function based on particle
! data.  looks for rawpart_chk.dat and pmass0.dat and def_chk.dat.
! Should rawpart_chk.dat not be present, it uses prho.dat instead.
! if that is also not present, uses pmass0 to reconstruct initial conditions.
! all correlations are linearly extrapolated to z=0.
!
! to read SGI data on an alpha, compile with -convert big_endian
!
integer, parameter :: n=32,nred=1,ng=nred*n,np=n**3
real*4, dimension(n,n,n) :: def,prho
real*4, dimension(np) :: prho1
equivalence(prho,prho1)
real*4, dimension(ng+2,ng,ng) :: rho
real*4, dimension(6,np) :: xv
complex*8, dimension(ng/2+1,ng,ng) :: crho
equivalence(crho,rho)
real, dimension(ng)::pk,weight
real*8 a
logical fxv
integer i, idim,ix,iy,iz,ixp,iyp,izp,ir
real wx,wy,wz,w000,w001,w010,w011,w100,w101,w110 &
       ,w111,p1,rhoc,tau,alambda,pi,amu,engyu
real, dimension(3) :: dx3,rho000,rho001,rho010,rho011,rho100,rho101 &
     ,rho110,rho111
  
pi=4*atan(1.)

open(10,file='rawpart_chk.dat',form='unformatted',status='old', err=100)
read(10) xv
close(10)
open(10,file='pmass0.dat',form='unformatted',status='old')
read(10) prho1
close(10)
fxv=.true.
goto 200
100 continue
fxv=.false.
if (.true.) then
   open(10,file='prho.dat',form='unformatted',status='old', err=300)
   read(10) prho
else
   ! read the gas density field
   open(10,file='u_chk.dat',form='unformatted',status='old', err=300)
   read(10) (prho1(i),u2,u3,u4,u5,i=1,n**3)
endif
close(10)
goto 200
300 continue
open(10,file='pmass0.dat',form='unformatted',status='old')
read(10) prho1
read(10) a
close(10)
a=a*5./3.
def=0
goto 400
200 continue
open(10,file='def_chk.dat',form='unformatted',status='old')
read(10) def
close(10)
open(10,file='state_chk.dat',form='unformatted',status='old')
read(10) loopni,a
close(10)
400 continue
write(*,*) 'a,z=',a,1/a-1
rho=0
if (fxv) then
   do i=1,np
      x=xv(1,i)
      y=xv(2,i)
      z=xv(3,i)
      ix=x
      iy=y
      iz=z
      wx=x-ix
      wy=y-iy
      wz=z-iz
      ix=modulo(ix-1,n)+1
      iy=modulo(iy-1,n)+1
      iz=modulo(iz-1,n)+1
      w111=wx*wy*wz
      w110=wx*wy*(1-wz)
      w101=wx*(1-wy)*wz
      w100=wx*(1-wy)*(1-wz)
      w011=(1-wx)*wy*wz
      w010=(1-wx)*wy*(1-wz)
      w001=(1-wx)*(1-wy)*wz
      w000=(1-wx)*(1-wy)*(1-wz)
      ixp=mod(ix,n)+1
      iyp=mod(iy,n)+1
      izp=mod(iz,n)+1   
      rho000=xpos(ix,iy,iz)
      rho001=xpos(ix,iy,izp)
      rho010=xpos(ix,iyp,iz)
      rho011=xpos(ix,iyp,izp)
      rho100=xpos(ixp,iy,iz)
      rho101=xpos(ixp,iy,izp)
      rho110=xpos(ixp,iyp,iz)
      rho111=xpos(ixp,iyp,izp)
      dx3=rho000*w000+rho001*w001+rho010*w010  &
           +rho011*w011+rho100*w100+rho101*w101       &
           +rho110*w110+rho111*w111
      x=modulo(x+dx3(1),n*1.)*nred+1
      y=modulo(y+dx3(2),n*1.)*nred+1
      z=modulo(z+dx3(3),n*1.)*nred+1
      p1=prho1(i)
      call pmap
   enddo
else
   do k=1,n
      do j=1,n
         do i=1,n
            ip=mod(i,n)+1
            im=mod(i+n-2,n)+1
            jp=mod(j,n)+1
            jm=mod(j+n-2,n)+1
            kp=mod(k,n)+1
            km=mod(k+n-2,n)+1
            dx=(def(ip,j,k)-def(im,j,k))/2
            dy=(def(i,jp,k)-def(i,jm,k))/2
            dz=(def(i,j,kp)-def(i,j,km))/2
            x=modulo(i+dx,n*1.)*nred+1
            y=modulo(j+dy,n*1.)*nred+1
            z=modulo(k+dz,n*1.)*nred+1
            p1=prho(i,j,k)
            call pmap
         enddo
      end do
   end do
endif
call sfft3(rho,ng,1)
crho=crho*conjg(crho)/ng**3
if (fxv) then
   anoise=sum(prho**2)/2/ng**3
   crho=crho-anoise
   write(*,*) 'using rawpart_chk.dat with poisson noise removal',anoise
end if
pk=0
weight=0
do k=1,ng
   kk=k-1
   if (k>ng/2) kk=kk-ng
   do j=1,ng
      jj=j-1
      if (j>ng/2) jj=jj-ng
      do i=1,ng/2+1
         ii=i-1
         r=sqrt(ii**2+jj**2+kk**2*1.)
         ir=r+1
         w=r+1-ir
         irp=ir+1
         weight(ir)=weight(ir)+(1-w)
         weight(irp)=weight(irp)+w
         pk(ir)=pk(ir)+crho(i,j,k)*(1-w)
         pk(irp)=pk(irp)+crho(i,j,k)*w
      end do
   end do
end do
where (weight>0) pk=pk/weight
open(10,file='pspect.dat')
do i=1,ng
   write(10,*) i-1,(i-1)**3*pk(i)*4*pi/ng**3/a**2,weight(i)
end do
close(10)
crho(1,1,1)=0
call sfft3(rho,ng,-1)
pk=0
weight=0
do k=1,ng
   kk=k-1
   if (k>ng/2) kk=kk-ng
   do j=1,ng
      jj=j-1
      if (j>ng/2) jj=jj-ng
      do i=1,ng
         ii=i-1
         if (i>ng/2) ii=ii-ng
         r=sqrt(ii**2+jj**2+kk**2*1.)
         ir=r+1
         w=r+1-ir
         irp=ir+1
         weight(ir)=weight(ir)+(1-w)
         weight(irp)=weight(irp)+w
         pk(ir)=pk(ir)+rho(i,j,k)*(1-w)
         pk(irp)=pk(irp)+rho(i,j,k)*w
      end do
   end do
end do
where (weight>0) pk=pk/weight
open(10,file='corr.dat')
do i=1,ng
   write(10,*) i-1,pk(i)/a**2,weight(i)
end do
close(10)

contains
  function xpos(i,j,k)
    real, dimension(3)::xpos
    
    ip=mod(i,n)+1
    im=mod(i-2+n,n)+1
    jp=mod(j,n)+1
    jm=mod(j-2+n,n)+1
    kp=mod(k,n)+1
    km=mod(k-2+n,n)+1
    xpos(1)=(def(ip,j,k)-def(im,j,k))/2
    xpos(2)=(def(i,jp,k)-def(i,jm,k))/2
    xpos(3)=(def(i,j,kp)-def(i,j,km))/2
    return
  end function xpos

  subroutine pmap
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
         ixp=mod(ix,ng)+1
         iyp=mod(iy,ng)+1
         izp=mod(iz,ng)+1
         p1=p1*nred**3
         rho(ix,iy,iz)=rho(ix,iy,iz)+w000*p1
         rho(ix,iy,izp)=rho(ix,iy,izp)+w001*p1
         rho(ix,iyp,iz)=rho(ix,iyp,iz)+w010*p1
         rho(ix,iyp,izp)=rho(ix,iyp,izp)+w011*p1
         rho(ixp,iy,iz)=rho(ixp,iy,iz)+w100*p1
         rho(ixp,iy,izp)=rho(ixp,iy,izp)+w101*p1
         rho(ixp,iyp,iz)=rho(ixp,iyp,iz)+w110*p1
         rho(ixp,iyp,izp)=rho(ixp,iyp,izp)+w111*p1

  end subroutine pmap
end
