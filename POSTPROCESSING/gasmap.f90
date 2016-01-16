! note: compile with -r8
! to avoid float overflows in unit conversions.
!
integer, parameter :: n=128,np=n**3
real*4, dimension(n,n,n) :: def
real*4, dimension(5,n,n,n) :: u
real*4, dimension(3,n,n,n) :: pos

open(10,file='u_chk.dat',form='unformatted',status='old')
read(10) u
close(10)

open(10,file='def_chk.dat',form='unformatted',status='old')
read(10) def
close(10)

call xvmap(pos,u,def,n)
pos=modulo(pos-1,n*1.)+1.
open(10,file='gascartesian.dat')
write(10,'(5f12.4)') (pos(:,i,1,1),u(1,i,1,1),u(5,i,1,1),i=1,np)
! outputs x,y,z,mass,pressure
close(10)

end

subroutine xvmap(pos,u,def,n)
  implicit none
  integer n
  real*4 def(n,n,n)
  real*4 pos(3,n,n,n),u(5,n,n,n)
  real a
! locals

  pos(1,:,:,:)=cshift(def,1)-cshift(def,-1)
  pos(2,:,:,:)=cshift(def,1,2)-cshift(def,-1,2)
  pos(3,:,:,:)=cshift(def,1,3)-cshift(def,-1,3)
  pos=pos/2
  u(5,:,:,:)=u(5,:,:,:)-sum(u(2:3,:,:,:)**2,1)/2/u(1,:,:,:)
return
end
