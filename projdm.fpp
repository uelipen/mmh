      subroutine projdm(def)
      include 'relaxgl.fi'
      include 'nbody.fi'
      integer*4 idum
      real*4 part(6)
      include 'proj.fi'
      parameter(no=ndimproj)
      real def(ng1,ng2,ng3)
      real xpos(3,ng1,ng2,ng3), dx(3)
      real*4 rho(no,no),rholocal(no,no)
      integer icount
      logical firsttime
      save firsttime,icount
      data firsttime /.true./
      data icount /0/
      integer i,j,k
#ifdef _ALPHA
       character*24   fdate
       external       fdate
#endif

      nsub=no/ng1
      if (nsub*ng1 .ne. no) then
          write(*,*) 'projdm: nsub*n <> no',nsub,ng1,no
       endif
      icount=icount+1
      if (firsttime) then
      firsttime = .false.
      endif
c this first loop determines the memory placement.

#ifdef _ALPHA
      write(*,*) 'starting projdm ',fdate()
#endif
c$omp parallel do default(none) shared(xpos,def)
      do j=1,ng2
c xpos(1,1,1,1)=def(2,1,1)-def(1,1,1)
      xpos(1,:,j,:)=cshift(def(:,j,:),1)-def(:,j,:) 
      xpos(3,:,j,:)=cshift(def(:,j,:),1,2)-def(:,j,:)
      enddo
c$omp parallel do default(none) shared(xpos,def)
      do k=1,ng3
      xpos(2,:,:,k)=cshift(def(:,:,k),1,2)-def(:,:,k)
      enddo
c$omp parallel do default(none) shared(xpos)
      do i=1,ng1
      xpos(1,i,:,:)=(cshift(xpos(1,i,:,:),1,1)+xpos(1,i,:,:))/2
      xpos(1,i,:,:)=(cshift(xpos(1,i,:,:),1,2)+xpos(1,i,:,:))/2
      enddo
c$omp parallel do default(none) shared(xpos)
      do k=1,ng3
      xpos(3,:,:,k)=(cshift(xpos(3,:,:,k),1,1)+xpos(3,:,:,k))/2
      xpos(3,:,:,k)=(cshift(xpos(3,:,:,k),1,2)+xpos(3,:,:,k))/2
      enddo
c$omp parallel do default(none) shared(xpos)
      do j=1,ng2
      xpos(2,:,j,:)=(cshift(xpos(2,:,j,:),1,1)+xpos(2,:,j,:))/2
      xpos(2,:,j,:)=(cshift(xpos(2,:,j,:),1,2)+xpos(2,:,j,:))/2
      enddo
c make sure rho is distributed in memory
c$omp parallel do default(none) shared(rho)
      do i=1,no
      rho(:,i)=0
      enddo
c$omp parallel default(private) shared(xv,rho,xpos,nsub,icount)
      rholocal=0
c$omp do
      do i=1,npart
! note that the projections are computed at half grid cells:
         x=xv(1,i)-0.5
         y=xv(2,i)-0.5
         z=xv(3,i)-0.5
         ix=floor(x)
         iy=floor(y)
         iz=floor(z)
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
         ixp=modulo(ix,ng1)+1
         iyp=modulo(iy,ng2)+1
         izp=modulo(iz,ng3)+1
         ix=modulo(ix-1,ng1)+1
         iy=modulo(iy-1,ng2)+1
         iz=modulo(iz-1,ng3)+1
         
         do idim=1,3
            rho000=xpos(idim,ix,iy,iz)
            rho001=xpos(idim,ix,iy,izp)
            rho010=xpos(idim,ix,iyp,iz)
            rho011=xpos(idim,ix,iyp,izp)
            rho100=xpos(idim,ixp,iy,iz)
            rho101=xpos(idim,ixp,iy,izp)
            rho110=xpos(idim,ixp,iyp,iz)
            rho111=xpos(idim,ixp,iyp,izp)
            p1=rho000*w000+rho001*w001+rho010*w010	
     &           +rho011*w011+rho100*w100+rho101*w101	
     &           +rho110*w110+rho111*w111
            dx(idim)=p1
         enddo
         x=x+dx(1)
         y=y+dx(2)
         z=z+dx(3)
         ixo=modulo(int(x*nsub),no)+1
         iyo=modulo(int(y*nsub),no)+1
         izo=modulo(int(z*nsub),no)+1
         if (mod(icount,3) .eq. 0) then
            rholocal(ixo,iyo)=rholocal(ixo,iyo)+1
         else if (mod(icount,3) .eq. 1) then
            rholocal(ixo,izo)=rholocal(ixo,izo)+1
         else
            rholocal(iyo,izo)=rholocal(iyo,izo)+1
         endif
      enddo
c       write(*,*) 'projdm: enter critical'
c$omp critical (projdmc)
      rho=rho+rholocal
c$omp end critical (projdmc)
c       write(*,*) 'projdm: done critical'
c$omp end parallel
      rho=rho*no**2/npart
c      call warray(rho,no**2,107)
#ifdef _ALPHA
c      write(*,*) 'done projdm ',fdate()
#endif
      write(107) rho
      if (mod(icount,15) .eq. 0) then
         call wimage(rho,no,no,'pdm')
      endif
      end
      
