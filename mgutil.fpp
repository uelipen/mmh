
      subroutine prolong(u,uhalf,nx,ny,nz)
c WORKER routine
      implicit none
      integer nx,ny,nz
      real u(nx,ny,nz),uhalf((nx+1)/2,(ny+1)/2,(nz+1)/2)
cdir$ shared *u(:block,:block,:),*uhalf(:block,:block,:)
! locals
      integer i,j,k,ibl,ibh,jbl,jbh,kbl,kbh
c$omp parallel default(private) firstprivate(nx,ny,nz) shared(u,uhalf)
      do k=1,nz
cdir$ doshared(j,i) on u(i,j,1)      
c$omp do
         do j=1,ny
            do i=1,nx
               ibl=(i+1)/2
               ibh=(i/2)+1
               if (ibh .gt. nx/2) ibh=1
               jbl=(j+1)/2
               jbh=(j/2)+1
               if (jbh .gt. ny/2) jbh=1
               kbl=(k+1)/2
               kbh=(k/2)+1
               if (kbh .gt. nz/2) kbh=1
               u(i,j,k)=(uhalf(ibl,jbl,kbl)+uhalf(ibl,jbl,kbh)
     &              +uhalf(ibl,jbh,kbl)+uhalf(ibl,jbh,kbh)
     &              +uhalf(ibh,jbh,kbl)+uhalf(ibh,jbh,kbh)
     &              +uhalf(ibh,jbl,kbl)+uhalf(ibh,jbl,kbh))/8
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end



      subroutine inject(u,uhalf,nx,ny,nz)
c WORKER routine
      implicit none
      integer nx,ny,nz
      real u(nx,ny,nz),uhalf((nx+1)/2,(ny+1)/2,(nz+1)/2)
cdir$ shared *u(:block,:block,:),*uhalf(:block,:block,:)
! locals
      integer i,j,k,ibl,ibh,jbl,jbh,kbl,kbh

c$omp parallel default(private) firstprivate(nx,ny,nz) shared(u,uhalf)
      do k=1,nz
cdir$ doshared(j,i) on u(i,j,1)   
c$omp do   
         do j=1,ny
            do i=1,nx
               ibl=(i+1)/2
               ibh=(i/2)+1
               if (ibh .gt. nx/2) ibh=1
               jbl=(j+1)/2
               jbh=(j/2)+1
               if (jbh .gt. ny/2) jbh=1
               kbl=(k+1)/2
               kbh=(k/2)+1
               if (kbh .gt. nz/2) kbh=1
               u(i,j,k)=u(i,j,k)+(uhalf(ibl,jbl,kbl)+uhalf(ibl,jbl,kbh)
     &              +uhalf(ibl,jbh,kbl)+uhalf(ibl,jbh,kbh)
     &              +uhalf(ibh,jbh,kbl)+uhalf(ibh,jbh,kbh)
     &              +uhalf(ibh,jbl,kbl)+uhalf(ibh,jbl,kbh))/8
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end


      subroutine restrict(uhalf,u,nx,ny,nz)
c WORKER routine
      implicit none
      integer nx,ny,nz
      real u(nx,ny,nz),uhalf((nx+1)/2,(ny+1)/2,(nz+1)/2)
cdir$ shared *u(:block,:block,:),*uhalf(:block,:block,:)
! locals
      integer i,j,k,n1,n2,n3,ip,im,jp,jm,kp,km,i2,j2,k2
      n1=(nx+1)/2
      n2=(ny+1)/2
      n3=(nz+1)/2
c      if (n1 .gt. 128) write(*,*)'remove maxtrip directive'
c$dir scalar

c$omp parallel default(private) firstprivate(n1,n2,n3,nx,ny,nz)
c$omp& shared(u,uhalf)      
      do k=1,n3
c$dir force_parallel
cdir$ doshared(j,i) on u(i,j,1)  
c$omp do    
         do j=1,n2
            do i=1,n1
            j2=j*2-1
            k2=k*2-1
            jp=mod(j2+ny,ny)+1
            jm=mod(j2+ny-2,ny)+1
            kp=mod(k2,nz)+1
            km=mod(k2+nz-2,nz)+1
c$dir force_vector
               i2=i*2-1
               ip=mod(i2+nx,nx)+1
               im=mod(i2+nx-2,nx)+1
               uhalf(i,j,k)=u(i2,j2,k2)/8+(u(ip,j2,k2)+u(im,j2,k2)
     &  +u(i2,jp,k2)+u(i2,jm,k2)+u(i2,j2,kp)+u(i2,j2,km))/16
     &  +(u(ip,jp,k2)+u(ip,jm,k2)+u(im,jp,k2)+u(im,jm,k2)+u(ip,j2,kp)
     &   +u(ip,j2,km)+u(im,j2,kp)+u(im,j2,km)+u(i2,jp,kp)+u(i2,jp,km)
     &   +u(i2,jm,kp)+u(i2,jm,km))/32
     &  +(u(ip,jp,kp)+u(ip,jp,km)+u(ip,jm,kp)+u(ip,jm,km)
     &   +u(im,jp,kp)+u(im,jp,km)+u(im,jm,kp)+u(im,jm,km))/64
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end


      subroutine restrictu(uhalf,u,nx,ny,nz,nu)
c WORKER routine
      implicit none
      integer nx,ny,nz,nu
      real u(nu,nx,ny,nz),uhalf(1,(nx+1)/2,(ny+1)/2,(nz+1)/2)
cdir$ shared *u(:,:block,:block,:),*uhalf(:,:block,:block,:)
! locals
      integer i,j,k,n1,n2,n3,ip,im,jp,jm,kp,km,i2,j2,k2
      n1=(nx+1)/2
      n2=(ny+1)/2
      n3=(nz+1)/2
c      if (n1 .gt. 128) write(*,*)'remove maxtrip directive'
c$dir scalar   

c$omp parallel default(private) firstprivate(n1,n2,n3,nx,ny,nz)
c$omp& shared(u,uhalf)  
      do k=1,n3
c$dir force_parallel
cdir$ doshared(j,i) on u(1,i,j,1)      
c$omp do
         do j=1,n2
            do i=1,n1
            j2=j*2-1
            k2=k*2-1
            jp=mod(j2+ny,ny)+1
            jm=mod(j2+ny-2,ny)+1
            kp=mod(k2,nz)+1
            km=mod(k2+nz-2,nz)+1
c$dir force_vector
               i2=i*2-1
               ip=mod(i2+nx,nx)+1
               im=mod(i2+nx-2,nx)+1
        uhalf(1,i,j,k)=u(1,i2,j2,k2)/8+(u(1,ip,j2,k2)+u(1,im,j2,k2)
     &  +u(1,i2,jp,k2)+u(1,i2,jm,k2)+u(1,i2,j2,kp)+u(1,i2,j2,km))/16
     &  +(u(1,ip,jp,k2)+u(1,ip,jm,k2)+u(1,im,jp,k2)
     &  +u(1,im,jm,k2)+u(1,ip,j2,kp)
     &   +u(1,ip,j2,km)+u(1,im,j2,kp)+u(1,im,j2,km)
     &   +u(1,i2,jp,kp)+u(1,i2,jp,km)
     &   +u(1,i2,jm,kp)+u(1,i2,jm,km))/32
     &  +(u(1,ip,jp,kp)+u(1,ip,jp,km)+u(1,ip,jm,kp)+u(1,ip,jm,km)
     &   +u(1,im,jp,kp)+u(1,im,jp,km)+u(1,im,jm,kp)+u(1,im,jm,km))/64
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end

#if 0
      subroutine truncerr(trnerr,uhalf,residh,u,rho,defph
     &      ,dx,nx,ny,nz,iopt)
      implicit none
      integer nx,ny,nz,iopt
      real trnerr,u(nx,ny,nz),rho(nx,ny,nz)
     &      ,residh((nx+1)/2,(ny+1)/2,(nz+1)/2)
     &      ,uhalf((nx+1)/2,(ny+1)/2,(nz+1)/2)
     &      ,defph((nx+1)/2,(ny+1)/2,(nz+1)/2)
     &     ,dx
cdir$ shared *rho(:block,:block,:),*uhalf(:block,:block,:)
cdir$ shared *u(:block,:block,:),*residh(:block,:block,:)
cdir$ shared *defph(:block,:block,:)
c estimate truncation error.
c the idea is to restrict and prolong, and see how much one looses in
c the infinity norm.
c  The quantity of interest is || L[R[u]]-R[L[u]] ||_\infty
c the array RES contains R[\rho-L[u]]
c so all we need is L[R[u]]-R[\rho]
c   we build each of R[u], R[\rho] and run the relaxor with zero iterations.
c
      integer nxh,nyh,nzh,irelaxopt
      real rtmp

      nxh=(nx+1)/2
      nyh=(ny+1)/2
      nzh=(nz+1)/2
c just borrow uhalf for a second...      
      call restrict(uhalf,rho,nx,ny,nz)
      call matsub(residh,uhalf,residh,nxh,nyh,nzh)
c now residh contains R[L[u]].  It changed sign.      
      call restrict(uhalf,u,nx,ny,nz)
c infinity norm residual      
      irelaxopt=4+iopt*8
c returns |R[L[u]]-L[R[u]]|      
      call relax(uhalf,rtmp,residh,defph,trnerr,dx*2,nxh,nyh,nzh
     &    ,0,irelaxopt)
      return
      end
#endif




      subroutine matsub(a,b,c,nx,ny,nz)
c WORKER routine
      implicit none
      integer nx,ny,nz
      real a(nx,ny,nz), b(nx,ny,nz), c(nx,ny,nz)
cdir$ shared *a(:block,:block,:),*c(:block,:block,:)
cdir$ shared *b(:block,:block,:)
c potential violation of ANSI conventions: a,b,c may not be unique.
c locals
      integer i,j,k

c$omp parallel private(i,j,k) firstprivate(nx,ny,nz) shared(a,b,c)
      do k=1,nz
cdir$ doshared(j,i) on a(i,j,1)
c$omp do      
         do j=1,ny
            do i=1,nx
               a(i,j,k)=b(i,j,k)-c(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end


      subroutine usum(sum,arr,nx,ny,nz)
c WORKER routine
      implicit none
      integer nx,ny,nz
      real sum, arr(nx,ny,nz)
cdir$ shared *arr(:block,:block,:)
c locals
      integer i,j,k

      sum=0

c$omp parallel private(i,j,k) shared(nx,ny,nz,arr) reduction(+:sum)
      do k=1,nz
cdir$ doshared(j,i) on a(i,j,1)      
c$omp do
         do j=1,ny
            do i=1,nx
               sum=sum+arr(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end



      subroutine rgmult(arr,def,na1,na2,na3)
c WORKER routine.  expensive.
c multiple arr by \sqrt{g}.
c
      implicit none
      integer na1,na2,na3,k
      real arr(na1,na2,na3),dx ,def(na1,na2,na3)
cdir$ shared *def(:block,:block,:), *arr(:block,:block,:)
c locals
      integer kp,km,j,jp,jm,i,ip,im,ns1,ns2,ns3,ndim,ib,jb,kb
     &     ,io,ii,jo,ji,ko,ki,kbm,jbm,ibm
      real phixx,phiyy,phizz,phixy,phixz,phiyz,a11,a12,a13,a22,a23,a33
     &    ,det,dsum,asum,dfact

c$omp parallel default(private) firstprivate(na1,na2,na3)
c$omp& shared(def,arr)
      do k=1,na3
         kp=mod(k,na3)+1
         km=mod(k+na3-2,na3)+1
c$dir prefer_vector
cdir$ doshared(j,i) on arr(i,j,1)
c$omp do
         do j=1,na2
            jp=mod(j,na2)+1
            jm=mod(j+na2-2,na2)+1
            do i=1,na1
               ip=mod(i,na1)+1
               im=mod(i+na1-2,na1)+1
               phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
               phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
               phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
               phixy=(def(ip,jp,k)-def(im,jp,k)
     &              -def(ip,jm,k)+def(im,jm,k))/4
               phiyz=(def(i,jp,kp)-def(i,jp,km)
     &              -def(i,jm,kp)+def(i,jm,km))/4
               phixz=(def(ip,j,kp)-def(im,j,kp)
     &              -def(ip,j,km)+def(im,j,km))/4
               phixz=(def(ip,j,kp)-def(im,j,kp)
     &              -def(ip,j,km)+def(im,j,km))/4
               a11=1+phixx
               a12=phixy
               a13=phixz
               a22=1+phiyy
               a23=phiyz
               a33=1+phizz
               det=(a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &              -a12**2*a33)
               arr(i,j,k)=arr(i,j,k)*det
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end
            
 

      subroutine rgdiv(arr,def,na1,na2,na3)
c WORKER routine.  expensive.
c divide arr by \sqrt{g}.
c
      implicit none
      integer na1,na2,na3,k
      real arr(na1,na2,na3),dx ,def(na1,na2,na3)
cdir$ shared *def(:block,:block,:), *arr(:block,:block,:)
c locals
      integer kp,km,j,jp,jm,i,ip,im,ns1,ns2,ns3,ndim,ib,jb,kb
     &     ,io,ii,jo,ji,ko,ki,kbm,jbm,ibm
      real phixx,phiyy,phizz,phixy,phixz,phiyz,a11,a12,a13,a22,a23,a33
     &    ,det,dsum,asum,dfact

c$omp parallel default(private) firstprivate(na1,na2,na3)
c$omp& shared(def,arr)
      do k=1,na3
         kp=mod(k,na3)+1
         km=mod(k+na3-2,na3)+1
c$dir prefer_vector
cdir$ doshared(j,i) on arr(i,j,1)
c$omp do
         do j=1,na2
            jp=mod(j,na2)+1
            jm=mod(j+na2-2,na2)+1
            do i=1,na1
               ip=mod(i,na1)+1
               im=mod(i+na1-2,na1)+1
               phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
               phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
               phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
               phixy=(def(ip,jp,k)-def(im,jp,k)
     &              -def(ip,jm,k)+def(im,jm,k))/4
               phiyz=(def(i,jp,kp)-def(i,jp,km)
     &              -def(i,jm,kp)+def(i,jm,km))/4
               phixz=(def(ip,j,kp)-def(im,j,kp)
     &              -def(ip,j,km)+def(im,j,km))/4
               phixz=(def(ip,j,kp)-def(im,j,kp)
     &              -def(ip,j,km)+def(im,j,km))/4
               a11=1+phixx
               a12=phixy
               a13=phixz
               a22=1+phiyy
               a23=phiyz
               a33=1+phizz
               det=(a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &              -a12**2*a33)
               arr(i,j,k)=arr(i,j,k)/det
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end
            
 

c some diagnostic routines:
      subroutine diagarray(arr,nx,ny,nz)
      implicit none
      integer nx,ny,nz
      real arr(nx,ny,nz)
c locals
      integer i,j,k,nj,nji

      return
      
      nj=ny
      nji=1
      if (nj .gt. 4) nji=nj/2
      do i=1,nx,nji
        if (nx .gt. 1)  write(*,*) 'i=',i
        do j=1,nj,nji
           write(*,'(1x,6g13.6)') (arr(i,j,k),k=1,nz)
        enddo
      enddo
      return
      end











