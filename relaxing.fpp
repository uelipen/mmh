c -*- Fortran -*- 77 file relaxing.fpp
c**************************************************************************
c File: relaxing.fpp
c started Sept 3, 1994 by Ue-Li Pen, upen@astro.princeton.edu
c Purpose: to implement a single step of the relaxing TVD scheme in
c curvilinear coordinates.
c
c
c Feb 6, 1995: everything seems to work well, but we will now add the
c exact energy conservation expression.  This will require three extra
c arrays to implement to second order accuracy.
c     
c Feb 11: energy conservation seems to work.  We underestimate \sigma
c         and my current suspicion is that it stems from the TVD limiter.
c     
c     
c**************************************************************************
c Note that ng3 must be 2**n >= 16 for the relaxation to work.
c     
c     
#ifdef EXACTENERGY
      subroutine relaxing
      write(*,*) 'wrong subroutine relaxing for exact grav engy'
      stop
      end
      
      subroutine rhorelaxing(u,rhoflux,def,defp,cmax,dt, nfluidcmp)
#elif defined(_ALPHA)
      subroutine relaxing(u,def,defp,cmax,dt, nfluidcmp)
      include 'relaxgl.fi'

      real flat(ng1*ng2,
     &     (20*maxfluidcmp+7)*nstore+ncstore+1+1024/(ng1*ng2))

c locals
      logical firsttime
      save firsttime,ioff
      data firsttime /.true./

      integer iwordsize,ioff,ichk
cdec$ alias get_page_size, "getpagesize"
      integer get_page_size,ipagesize
      external get_page_size

      if (firsttime) then
         firsttime=.false.
      ipagesize=get_page_size()
      iwordsize=%loc(flat(2,1))-%loc(flat(1,1))
      write(*,*) 'wordsize=',iwordsize
      write(*,*) 'flat at ', mod(%loc(flat(1,1)),ipagesize),
     &     mod(%loc(flat(2,1)),ipagesize)
      ioff=(ipagesize-mod(%loc(flat(1,1)),ipagesize))/iwordsize+1
c      ioff=ioff+512
      write(*,*) 'ioff=',ioff
      write(*,*) 'flat ioff ', mod(%loc(flat(ioff,1)),ipagesize),
     &     mod(%loc(flat(ioff+1,1)),ipagesize)
      endif
      ichk= mod(%loc(flat(ioff,1)),ipagesize)
      if (ichk .ne. 0) write(*,*) 'relaxing ioff ',ichk
      call relaxing1(u,def,defp, flat(ioff,1), flat(ioff,7*nstore+1), 
     &     flat(ioff,(3*maxfluidcmp+7)*nstore+1), ! res
     &     flat(ioff,(4*maxfluidcmp+7)*nstore+1), ! u1
     &     flat(ioff,(5*maxfluidcmp+7)*nstore+1), ! vtf
     &     flat(ioff,(8*maxfluidcmp+7)*nstore+1), ! vtfstore
     &     flat(ioff,(20*maxfluidcmp+7)*nstore+1), ! c
     &     flat(ioff,(20*maxfluidcmp+7)*nstore+ncstore+1), ! u1store
     &     cmax,dt, nfluidcmp)
      return
      end
      subroutine relaxing1(u,def,defp, baj, vt, res, u1, vtf, 
     &     vtfstore, c, u1store, cmax,dt, nfluidcmp)


#else
      subroutine relaxing(u,def,defp,cmax,dt, nfluidcmp)
#endif
c
c we need to pass nfluidcmp as an argument, since I want to be able
c to simulate several phases if necessary.  E.g., one anisotropic dark
c matter plus one ideal gas.
c 'c' is the freezing soundspeed.  Since we use units where dx=1,
c  It must satisfy dt < 1/c.
c To calculate c, one needs to construct the whole metric, which
c is a sacrifice of computational efficiency over memory.  It should
c be a small effect compared to the time integrator.
c
      implicit none
      include 'relaxgl.fi'
      integer nfluidcmp
c
c both def and defp are defined at the half time step between the
c current u and the next one that this routines computes.
c
      real u(nfluidcmp,ng1,ng2,ng3), def(ng1,ng2,ng3), defp(ng1,ng2,ng3)
cmydist u(*,*,block,*), def(*,block,*)
cmydist  defp(*,block,*)
#ifdef EXACTENERGY
      real rhoflux(ng1,ng2,ng3,3)
#endif
      real cmax,dt
c
c Notes:
c    defp must be defined for the middle of the time step, so that the
c runge kutta has grid velocities known to second order.      
c
c locals:
c
c the parameter nstore is the number of tiers of metric and various
c other temporary stores required to achieve the second order
c runge-kutta time integrator.
c
c
c arrays:
c we assume everything is statically allocated.
c
      real baj(7,ng1,ng2,nstore), vt(ng1,ng2,nstore,3,maxfluidcmp)
     &    ,res(ng1,ng2,maxfluidcmp), u1(maxfluidcmp,ng1,ng2,nstore)
     &    ,vtf(ng1,ng2,nstore,3,maxfluidcmp)
c we need to retain a copy of the top layer:
      real vtfstore(ng1,ng2,4,3,maxfluidcmp),c(ng1,ng2,ncstore)
     &	  ,u1store(maxfluidcmp,ng1,ng2,4)
cmydist baj(*,*,block,*),vt(*,block,*,*,*),res(*,block,*)
cmydist u1(*,*,block,*),vtf(*,block,*,*,*)
cmydist vtfstore(*,block,*,*,*)
cmydist c(*,block,*),u1store(*,*,block,*)
      integer indx(ng3), indxu(ng3), indxc(ng3)
c
c array indx has image range (1:nstore), and maps the desired index
c onto the shorter temporary arrays.
c
c procedure scope:
      integer idxptr, idxcptr
c
c short range locals:
      integer k,kg,kh,kf,j,i,nf,idim
c

c----------------------------------------------------------------------
c begin executables: initialization
c      
#ifdef _ALPHANNN
cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(baj,vt,res,u1,vtf)
cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(vtfstore,c,u1store)
      write(*,*) '%c1=',%loc(c(1,1,1)),mod(%loc(c(1,1,1)),8192)
#endif


      do k=1,ng3
         indxu(k)=k
      enddo
c
c begin main body
c
      idxptr=0

c build the whole c-index pointers.  One could have done the same with
c the u and vt arrays.      
      idxcptr=0
c we will need the first six again at the end, and 7 running
      do k=1,ncstore
         indxc(k)=k
      enddo
      idxcptr=6
      do k=ncstore+1,ng3
         idxcptr=mod(idxcptr-6,ncstore-6)+7
         indxc(k)=idxcptr
      enddo

c bootstrap the whole process:
c we need the metric over a long range:
      do kg=1,4
         idxptr=mod(idxptr,nstore)+1
         indx(kg)=idxptr
         call calcbaj(baj,def,indx,kg)
         call cflkeuler(c,u,baj,defp,indxc,indx,kg)
         call flux(vt,defp,u,baj,indx,indxu,nfluidcmp,ng3,kg)
      enddo

cfpp$ skip
      do k=3,6
c dpv calculates the upwind TVD limited flux 
c the domain of dependence of dpv is 2 on each side.
         kg=k+2
         kh=k
         idxptr=mod(idxptr,nstore)+1
         indx(kg)=idxptr
         call calcbaj(baj,def,indx,kg)
         call cflkeuler(c,u,baj,defp,indxc,indx,kg)
         call flux(vt,defp,u,baj,indx,indxu,nfluidcmp,ng3,kg)
         call dpv(res,u,vt,c,indx,indxu,indxc,nfluidcmp,ng3,k)
c Runge-kutta half step
         call stepu(u1,u,res,dt/2,indx,nstore,nfluidcmp,kh)
         call flux(vtf,defp,u1,baj,indx,indx,nfluidcmp,nstore,kh)
      enddo
c we will need these 4 vtf layers again at the very end:

c$omp parallel default(private) 
c$omp& firstprivate(nfluidcmp) 
c$omp& shared(vtfstore,vtf,indx,u1store,u1)
cfpp$ skip
cc*$* assert do(serial)
      do nf=1,nfluidcmp
cfpp$ skip
         do idim=1,3
cfpp$ skip

            do k=1,4
c$dir prefer_parallel_ext
c*$* assert do (concurrent)
c$omp do
               do j=1,ng2
                  do i=1,ng1
                     vtfstore(i,j,k,idim,nf)=vtf(i,j,indx(k+2),idim,nf)
                  enddo
               enddo
c$omp enddo nowait
            enddo
         enddo
      enddo
cfpp$ skip
cc*$* assert do(serial)
!cdir novector
      do nf=1,nfluidcmp
         do k=1,4
c$dir prefer_parallel_ext               
c*$* assert do (concurrent)
c$omp do
            do j=1,ng2
!cdir vector
               do i=1,ng1
                  u1store(nf,i,j,k)=u1(nf,i,j,indx(k+2))
               enddo
            enddo
c$omp enddo nowait
         enddo
      enddo
c$omp end parallel

c main loop down the layers.  k=5 will be the first finished layer.
      do k=7,ng3+2
c the main iteration variable kg tracks the level at which we update
c the gravity
         kh=mod(k+ng3-1,ng3)+1
         kg=mod(kh+1,ng3)+1
         idxptr=mod(idxptr,nstore)+1
         indx(kg)=idxptr
c kh tracks the half step runge kutta update
         kf=mod(kh-3+ng3,ng3)+1
c kf tracks the full step runge kutta index.  At the end of the loop,
c layer kf is completely updated.
         call calcbaj(baj,def,indx,kg)
c
c flux() uses only the metric and u() at the current level, and returns
c the 3 x nfluid flux matrix at the current level.
         call cflkeuler(c,u,baj,defp,indxc,indx,kg)
         call flux(vt,defp,u,baj,indx,indxu,nfluidcmp,ng3,kg)
c
c dpv() uses u and vt at two levels each way.  It solves the system of
c constant coefficient equations
c
c  \dot{u}   + c (v^1,x+v^2,y) = 0
c  \dot{v^1} + c  u,x          = 0   ; -> split (F1-v^1)/\epsilon
c  \dot{v^2} + c        u,y    = 0   ;          (F2-v^2)/\epsilon
c
         call dpv(res,u,vt,c,indx,indxu,indxc,nfluidcmp,ng3,kh)
c Runge-Kutta half step
c
c update first order half step estimate:  u^(1/2)=u-dpv*dt/2
         call stepu(u1,u,res,dt/2,indx,nstore,nfluidcmp,kh)
c now redo to second order.
         call flux(vtf,defp,u1,baj,indx,indx,nfluidcmp,nstore,kh)
c         write(*,'(5f15.9)')(vtf(1,1,indx(kf),3,j),j=1,5)
#ifdef EXACTENERGY
         call rhodpv(res,rhoflux,u1,vtf,c,indx,indx,indxc,nfluidcmp
     &                  ,nstore,kf)
#else
         call dpv(res,u1,vtf,c,indx,indx,indxc,nfluidcmp,nstore,kf)
#endif
c FULL step         
c u^+ = u-dpv(u1)*dt
c         call stepu(u,u,res,dt,indxu,ng3,nfluidcmp,kf)
         call stepuu(u,res,dt,nfluidcmp,kf)
      enddo

c the last four layers need to use the buffer zone:
cfpp$ skip
      do kf=1,4
         kg=mod(kf+3,ng3)+1
         idxptr=mod(idxptr,nstore)+1
         indx(kg)=idxptr
c$omp parallel default(private) firstprivate(kf,nfluidcmp) 
c$omp& shared(indx,vtf,vtfstore,u1,u1store)
cfpp$ skip
cc*$* assert do prefer(serial)
         do nf=1,nfluidcmp
cfpp$ skip
            do idim=1,3

c$dir prefer_parallel_ext
c*$* assert do (concurrent)
c$omp do
               do j=1,ng2
                  do i=1,ng1
                   vtf(i,j,indx(kf+2),idim,nf)=vtfstore(i,j,kf,idim,nf)
                  enddo
               enddo
c$omp enddo nowait
            enddo
         enddo
cfpp$ skip
cc*$* assert do(serial)
!cdir novector
         do nf=1,nfluidcmp
c$dir prefer_parallel_ext               
c*$* assert do (concurrent)
c$omp do
            do j=1,ng2
!cdir vector
               do i=1,ng1
                  u1(nf,i,j,indx(kf+2))=u1store(nf,i,j,kf)
               enddo
            enddo
c$omp enddo nowait
         enddo
c$omp end parallel

#ifdef EXACTENERGY
         call rhodpv(res,rhoflux,u1,vtf,c,indx,indx,indxc,nfluidcmp
     &                  ,nstore,kf)
#else
         call dpv(res,u1,vtf,c,indx,indx,indxc,nfluidcmp,nstore,kf)
#endif
c FULL step         
         call stepuu(u,res,dt,nfluidcmp,kf)
      enddo
      
c that was simple!  we are done.
c
      return
      end


c**********************************************************************
#ifdef _ALPHA
      subroutine dpv(res,u,vt,c,indxr, indxu, indxc, nfluid, ngu3, k)
      include 'relaxgl.fi'

      real flat(ng1*ng2,
     &     15*maxfluidcmp+1+1024/(ng1*ng2))

c locals
      real*4 dtime, delta,systime,stime,tarray(2)
      integer*4 time, itarray(9), itime
      external dtime,time,ltime
      integer iwordsize,ioff,ichk
      logical firsttime
      save firsttime,ioff,itcum
      data firsttime /.true./
cdec$ alias get_page_size, "getpagesize"
      integer get_page_size,ipagesize
      external get_page_size

      if (firsttime) then
         firsttime=.false.
      ipagesize=get_page_size()
      iwordsize=%loc(flat(2,1))-%loc(flat(1,1))
      write(*,*) 'wordsize=',iwordsize
      write(*,*) 'flat at ', mod(%loc(flat(1,1)),ipagesize),
     &     mod(%loc(flat(2,1)),ipagesize)
      ioff=(ipagesize-mod(%loc(flat(1,1)),ipagesize))/iwordsize+1
      write(*,*) 'ioff=',ioff
      write(*,*) 'flat ioff ', mod(%loc(flat(ioff,1)),ipagesize),
     &     mod(%loc(flat(ioff+1,1)),ipagesize)
      endif
c      ichk= mod(%loc(flat(ioff,1)),ipagesize)
c      if (ichk .ne. 0) write(*,*) 'flat ioff ',ichk,k



c      call ltime(itime,itarray)
c      if (k .eq. 3) write(*,*) itime,itarray
c      delta=dtime(tarray)
c      write(*,*) delta
c      systime=tarray(2)
      call  dpva(res,u,vt,c,
     &     flat(ioff,1), ! vk
     &     flat(ioff,(5*maxfluidcmp)+1), ! vij
     &     flat(ioff,(6*maxfluidcmp)+1), ! u1
     &     flat(ioff,(11*maxfluidcmp)+1), ! u2
     &     flat(ioff,(12*maxfluidcmp)+1), ! flz
     &     flat(ioff,(14*maxfluidcmp)+1), ! flx
     &     indxr, indxu, indxc, nfluid, ngu3, k)
c      delta=dtime(tarray)

!      write(*,*) k
      return
      end
      subroutine dpva(res,u,vt,c,
     &     vk,vij,u1,u2,flz,flx,
     &     indxr, indxu, indxc, nfluid, ngu3, k)

#else
      subroutine dpv(res,u,vt,c,indxr, indxu, indxc, nfluid, ngu3, k)
#endif
      implicit none
      include 'relaxgl.fi'
      integer k, indxu(ng3), indxr(ng3), indxc(ng3), ngu3,nfluid
      real res(ng1,ng2,nfluid), u(nfluid,ng1,ng2,ngu3)
     &     , vt(ng1,ng2,nstore,3,nfluid), c(ng1,ng2,ncstore)
cmydist u(*,*,block,*),vt(*,block,*,*,*),res(*,block,*)
cmydist c(*,block,*)

c
c locals
      integer i,j,nf,l,lk,im,jm,kp,ip,jp,km, nf5, lkp, lkm
      parameter(nf5=5)
c 5 is the number of buffer zones for a domain of dependence of 2.
      real  vk(5,maxfluidcmp,ng1,ng2)
     &    ,vij(ng1,ng2,maxfluidcmp), u1(5,maxfluidcmp,ng1,ng2)
     &    ,u2(ng1,ng2,maxfluidcmp), flz(2,ng1,ng2,maxfluidcmp)
     &    ,flx(ng1,ng2,maxfluidcmp),fp,fm
cmydist vk(*,*,*,block),vij(*,block,*),u1(*,*,*,block)
cmydist u2(*,block,*),flz(*,*,block,*),flx(*,block,*)
c to conserve memory, one could equivalence some of the arrays.  But
c I am worried that equivalences break many optimizations and
c parallelizations.
      integer lastk
      save lastk
      data lastk /0/
#ifdef COLD    
      include 'cold.fi'
#endif      
#ifdef _ALPHA


cdec$ alias addr_to_rad, "_OtsAddrToRad"
      integer  addr_to_rad,irad,irad1,irad2
      external addr_to_rad
      integer omp_get_thread_num,cpuid,cpu_get_rad
      external omp_get_thread_num,cpuid,cpu_get_rad

ccdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(vk,vij,u1,u2,flz,flx)

#endif

           
      if (nfluid .ne. nf5) then
         write(*,*) 'dpv: only nfluid=5 currently support, nf=',nfluid
         stop
      endif
c first sweep the k dimension:
c zpj: I make c$doacross comments because it interfere with c$omp
cc$doacross local(j,nf,l,i,lk)
c*KAP*parallel region local(j,nf,l,i,lk)
c*KAP*parallel do

c$omp parallel do default(private) firstprivate(k) 
c$omp& shared(vk,vt,c,indxr,indxc,u,u1,indxu)
      do j=1,ng2

         do nf=1,nf5
c l covers the domain of dependence
            do l=1,5
               lk = mod(k+l-4+ng3,ng3)+1
               do i=1,ng1
                  vk(l,nf,i,j)=vt(i,j,indxr(lk),3,nf)/c(i,j,indxc(lk))
               enddo
            enddo
         enddo
         do l=1,5
         lk = mod(k+l-4+ng3,ng3)+1
            do i=1,ng1
                  u1(l,1,i,j)=u(1,i,j,indxu(lk))
                  u1(l,2,i,j)=u(2,i,j,indxu(lk))
                  u1(l,3,i,j)=u(3,i,j,indxu(lk))
                  u1(l,4,i,j)=u(4,i,j,indxu(lk))
                  u1(l,5,i,j)=u(5,i,j,indxu(lk))
            enddo
         enddo
#if 0
            if (k .eq. 1) then
               irad1=addr_to_rad(u1(1,1,1,j))
               irad2=addr_to_rad(u(1,1,j,1))
               irad=cpu_get_rad()
               if ((irad1 .ne. irad) .or. (irad2 .ne. irad)) then
                       write(*,*) j,irad1,irad2,irad
               endif
            endif
#endif
      enddo
c$omp end parallel do

c*KAP*end parallel region
#ifdef _ALPHA
c      if (k .eq. 1)write(*,*) '%vk=',mod(%loc(vk(1,1,1,1)),8192)
#endif
      call wphalfz(flz,u1,vk,nfluid)

      kp=mod(k,ng3)+1
      km=mod(k+ng3-2,ng3)+1
      lk = mod(k+3-4+ng3,ng3)+1
      lkm=mod(k+2-4+ng3,ng3)+1
      lkp=mod(k+4-4+ng3,ng3)+1
      
cc*$* assert do (concurrent)

c zpj: I make c$doacross comments
cc$doacross local(j,nf,i) shared(ng2,ng1,res,indxc,kp,k,c,km,flz
cc$&  ,vij,vt,indxr,u,u2,indxu)

c*KAP*parallel region local(j,nf,i)
c*KAP*parallel do
c zpj: end
#ifdef COLD
c$omp parallel do default(private) 
c$omp& shared(flz,vt,indxr,indxc,c,res,vij,u2,u,indxu,cold, 
c$omp& k,kp,km,lk,lkm,lkp)
#else
c$omp parallel do default(private) 
c$omp& shared(flz,vt,indxr,indxc,c,res,vij,u2,u,indxu,
c$omp& k,kp,km,lk,lkm,lkp)
#endif 
      do j=1,ng2
#ifdef COLD         
         do i=1,ng1
            if (cold(i,j,k) .and. cold(i,j,kp) ) then
               do nf=1,5
                  flz(2,i,j,nf)=(vt(i,j,indxr(lk),3,nf)+
     &      vt(i,j,indxr(lkp),3,nf))/(c(i,j,indxc(kp))+c(i,j,indxc(k)))
               enddo
            endif
            if (cold(i,j,km) .and. cold(i,j,k) ) then
               do nf=1,5
                  flz(1,i,j,nf)=(vt(i,j,indxr(lk),3,nf)+
     &      vt(i,j,indxr(lkm),3,nf))/(c(i,j,indxc(km))+c(i,j,indxc(k)))
               enddo
            endif
         enddo
#endif            
         do nf=1,5
            do i=1,ng1
               res(i,j,nf)=(c(i,j,indxc(kp))+c(i,j,indxc(k)))/2*
     &              flz(2,i,j,nf)-(c(i,j,indxc(k))+c(i,j,indxc(km)))/2
     &              *flz(1,i,j,nf)
               vij(i,j,nf)=vt(i,j,indxr(k),1,nf)/c(i,j,indxc(k))
               u2(i,j,nf)=u(nf,i,j,indxu(k))
            enddo
         enddo
      enddo
c$omp end parallel do
c*KAP*end parallel region
      call wphalfx(flx,u2,vij,nfluid)
cc*$* assert do (concurrent)s
c zpj: I make c$doacross comment
cc$doacross local(j,nf,i,ip,im,fp,fm)
#ifdef COLD
c$omp parallel do default(private) 
c$omp& shared(flx,vt,indxr,c,indxc,res,vij,cold,k,lk)
#else
c$omp parallel do default(private) 
c$omp& shared(flx,vt,indxr,c,indxc,res,vij,k,lk)
#endif
      do j=1,ng2
#ifdef COLD         
         do i=1,ng1
            ip=mod(i,ng1)+1
            im=mod(i-2+ng1,ng1)+1
            if (cold(i,j,k) .and. cold(ip,j,k) ) then
               do nf=1,5
                  flx(i,j,nf)=(vt(i,j,indxr(lk),1,nf)+
     &      vt(ip,j,indxr(lk),1,nf))/(c(ip,j,indxc(k))+c(i,j,indxc(k)))
               enddo
            endif
         enddo
#endif            
         do nf=1,nf5
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i-2+ng1,ng1)+1
               fp=(c(ip,j,indxc(k))+c(i,j,indxc(k)))/2*flx(i,j,nf)
               fm=(c(i,j,indxc(k))+c(im,j,indxc(k)))/2*flx(im,j,nf)
               res(i,j,nf)=res(i,j,nf)+fp-fm
            enddo
         enddo

c lastly the y dimension flux:
         do nf=1,nf5
            do i=1,ng1
               vij(i,j,nf)=vt(i,j,indxr(k),2,nf)/c(i,j,indxc(k))
            enddo
         enddo
      enddo
c$omp end parallel do
      call wphalfy(flx,u2,vij,nfluid)
#ifdef COLD
c$omp parallel do default(private) 
c$omp& shared(flx,vt,indxr,c,indxc,res,cold,k,lk)
#else
c$omp parallel do default(private) 
c$omp& shared(flx,vt,indxr,c,indxc,res,k,lk)
#endif
#ifdef COLD         
c*$* assert do (concurrent)
      do j=1,ng2
         jp=mod(j,ng2)+1
         jm=mod(j-2+ng2,ng2)+1
         do i=1,ng1
            if (cold(i,j,k) .and. cold(i,jp,k) ) then
               do nf=1,5
                  flx(i,j,nf)=(vt(i,j,indxr(lk),2,nf)+
     &      vt(i,jp,indxr(lk),2,nf))/(c(i,jp,indxc(k))+c(i,j,indxc(k)))
               enddo
            endif
         enddo
      enddo
c$omp end parallel do
c$omp parallel do default(private) 
c$omp& shared(flx,vt,indxr,c,indxc,res,cold,k,lk)
#endif            
c*$* assert do (concurrent)
      do j=1,ng2
         do nf=1,nf5
            do i=1,ng1
               jp=mod(j,ng2)+1
               jm=mod(j-2+ng2,ng2)+1
               res(i,j,nf)=res(i,j,nf)
     &              +(c(i,jp,indxc(k))+c(i,j,indxc(k)))/2*flx(i,j,nf)
     &              -(c(i,j,indxc(k))+c(i,jm,indxc(k)))/2*flx(i,jm,nf)
            enddo
         enddo
      enddo
c$omp end parallel do

      if (k .gt. lastk) then
c         write(*,'(5f15.9)') (flz(1,1,2,j),j=1,5)
         lastk=k
      else
c          write(*,'(5f15.9)') (vk(1,1,5,j),j=1,5)
      endif
      return
      end


#ifdef EXACTENERGY



c**********************************************************************

      subroutine rhodpv(res,rhoflux,u,vt,c,indxr, indxu, indxc, nfluid,
     &     ngu3, k)
      implicit none
      include 'relaxgl.fi'
      integer k, indxu(ng3), indxr(ng3), indxc(ng3), ngu3,nfluid
      real res(ng1,ng2,nfluid), u(nfluid,ng1,ng2,ngu3)
     &     , vt(ng1,ng2,nstore,3,nfluid), c(ng1,ng2,ncstore)
     &      , rhoflux(ng1,ng2,ng3,3)
cmydist u(*,*,block,*),vt(*,block,*,*,*),res(*,block,*)
cmydist c(*,block,*)      
c
c locals
      integer i,j,nf,l,lk,im,jm,kp,ip,jp,km, nf5, lkp, lkm
      parameter(nf5=5)
c 5 is the number of buffer zones for a domain of dependence of 2.
      real  vk(5,maxfluidcmp,ng1,ng2)
     &    ,vij(ng1,ng2,maxfluidcmp), u1(5,maxfluidcmp,ng1,ng2)
     &    ,u2(ng1,ng2,maxfluidcmp), flz(2,ng1,ng2,maxfluidcmp)
     &    ,flx(ng1,ng2,maxfluidcmp),fp,fm
cmydist vk(*,*,*,block),vij(*,block,*),u1(*,*,*,block)
cmydist u2(*,block,*),flz(*,*,block,*),flx(*,block,*)
c to conserve memory, one could equivalence some of the arrays.  But
c I am worried that equivalences break many optimizations and
c parallelizations.
      integer lastk
      save lastk
      data lastk /0/
#ifdef COLD
      include 'cold.fi'
#endif      
           
      if (nfluid .ne. nf5) then
         write(*,*) 'dpv: only nfluid=5 currently support, nf=',nfluid
         stop
      endif
c first sweep the k dimension:
c$doacross local(j,nf,l,i,lk)

c$omp parallel do default(private)
c$omp& shared(k,vk,vt,indxr,c,indxc,u1,indxu))
      do j=1,ng2
         do nf=1,nf5
c l covers the domain of dependence
            do l=1,5
               lk = mod(k+l-4+ng3,ng3)+1
               do i=1,ng1
                  vk(l,nf,i,j)=vt(i,j,indxr(lk),3,nf)/c(i,j,indxc(lk))
               enddo
            enddo
         enddo
         do l=1,5
         lk = mod(k+l-4+ng3,ng3)+1
            do i=1,ng1
                  u1(l,1,i,j)=u(1,i,j,indxu(lk))
                  u1(l,2,i,j)=u(2,i,j,indxu(lk))
                  u1(l,3,i,j)=u(3,i,j,indxu(lk))
                  u1(l,4,i,j)=u(4,i,j,indxu(lk))
                  u1(l,5,i,j)=u(5,i,j,indxu(lk))
            enddo
         enddo
      enddo
      call wphalfz(flz,u1,vk,nfluid)
      kp=mod(k,ng3)+1
      km=mod(k+ng3-2,ng3)+1
      lk = mod(k+3-4+ng3,ng3)+1
      lkm=mod(k+2-4+ng3,ng3)+1
      lkp=mod(k+4-4+ng3,ng3)+1
      
cc*$* assert do (concurrent)
c$doacross local(j,nf,i) shared(ng2,ng1,nf,res,x,indxc,kp,k,c,km,flz
c$&  ,vij,vt,indxr,u,u2,indxu)
#ifdef COLD
c$omp parallel default(private) firstprivate(k,kp,km,lk,lkm,lkp)
c$omp& shared(cold,flz,vt,indxr,c,indxc,res,vij,u2,indxu,rhoflux)
#else
c$omp parallel default(private) firstprivate(k,kp,km,lk,lkm,lkp)
c$omp& shared(flz,vt,indxr,c,indxc,res,vij,u2,indxu,rhoflux)
#endif
c$omp do
      do j=1,ng2
#ifdef COLD         
         do i=1,ng1
            if (cold(i,j,k) .and. cold(i,j,kp) ) then
               do nf=1,5
                  flz(2,i,j,nf)=(vt(i,j,indxr(lk),3,nf)+
     &      vt(i,j,indxr(lkp),3,nf))/(c(i,j,indxc(kp))+c(i,j,indxc(k)))
               enddo
            endif
            if (cold(i,j,km) .and. cold(i,j,k) ) then
               do nf=1,5
                  flz(1,i,j,nf)=(vt(i,j,indxr(lk),3,nf)+
     &      vt(i,j,indxr(lkm),3,nf))/(c(i,j,indxc(km))+c(i,j,indxc(k)))
               enddo
            endif
         enddo
#endif            
         do nf=1,5
            do i=1,ng1
               res(i,j,nf)=(c(i,j,indxc(kp))+c(i,j,indxc(k)))/2*
     &              flz(2,i,j,nf)-(c(i,j,indxc(k))+c(i,j,indxc(km)))/2
     &              *flz(1,i,j,nf)
            enddo
         enddo
         do nf=1,5
            do i=1,ng1
               vij(i,j,nf)=vt(i,j,indxr(k),1,nf)/c(i,j,indxc(k))
            enddo
         enddo
         do i=1,ng1
               u2(i,j,1)=u(1,i,j,indxu(k))
               u2(i,j,2)=u(2,i,j,indxu(k))
               u2(i,j,3)=u(3,i,j,indxu(k))
               u2(i,j,4)=u(4,i,j,indxu(k))
               u2(i,j,5)=u(5,i,j,indxu(k))
         enddo
      enddo
c$omp enddo nowait


      nf=1
c$dir prefer_parallel_ext            
cc$doacross local(i,j)

c$omp do
      do j=1,ng2
         do i=1,ng1
            rhoflux(i,j,k,3)=rhoflux(i,j,k,3)
     &          +(c(i,j,indxc(kp))+c(i,j,indxc(k)))/2*  flz(2,i,j,nf)
         enddo
      enddo
c$omp end parallel
      
      call wphalfx(flx,u2,vij,nfluid)
cc*$* assert do (concurrent)
c$doacross local(j,nf,i,ip,im,fp,fm)
#ifdef COLD
c$omp parallel default(private)
c$omp& shared(cold,k,flx,vt,indxr,lk,c,indxc,res,rhoflux)
#else
c$omp parallel default(private)
c$omp& shared(k,flx,vt,indxr,lk,c,indxc,res,rhoflux)
#endif
c$omp do
      do j=1,ng2
#ifdef COLD         
         do i=1,ng1
            ip=mod(i,ng1)+1
            im=mod(i-2+ng1,ng1)+1
            if (cold(i,j,k) .and. cold(ip,j,k) ) then
               do nf=1,5
                  flx(i,j,nf)=(vt(i,j,indxr(lk),1,nf)+
     &      vt(ip,j,indxr(lk),1,nf))/(c(ip,j,indxc(k))+c(i,j,indxc(k)))
               enddo
            endif
         enddo
#endif            
         do nf=1,nf5
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i-2+ng1,ng1)+1
               fp=(c(ip,j,indxc(k))+c(i,j,indxc(k)))/2*flx(i,j,nf)
               fm=(c(i,j,indxc(k))+c(im,j,indxc(k)))/2*flx(im,j,nf)
               res(i,j,nf)=res(i,j,nf)+fp-fm
            enddo
         enddo

c lastly the y dimension flux:
         do nf=1,nf5
            do i=1,ng1
               vij(i,j,nf)=vt(i,j,indxr(k),2,nf)/c(i,j,indxc(k))
            enddo
         enddo
      enddo
c$omp enddo nowait


      nf=1
c$dir prefer_parallel_ext            
cc$doacross local(i,j)

c$omp do
      do j=1,ng2
         do i=1,ng1
            ip=mod(i,ng1)+1
            rhoflux(i,j,k,1)=rhoflux(i,j,k,1)
     &            +(c(ip,j,indxc(k))+c(i,j,indxc(k)))/2  *flx(i,j,nf)
         enddo
      enddo
c$omp end parallel
      
      call wphalfy(flx,u2,vij,nfluid)

#ifdef COLD         
c*$* assert do (concurrent)
c$omp parallel do default(private) 
c$omp& shared(cold,flx,vt,indxr,indxc,res,k,lk)
      do j=1,ng2
         jp=mod(j,ng2)+1
         jm=mod(j-2+ng2,ng2)+1
         do i=1,ng1
            if (cold(i,j,k) .and. cold(i,jp,k) ) then
               do nf=1,5
                  flx(i,j,nf)=(vt(i,j,indxr(lk),2,nf)+
     &      vt(i,jp,indxr(lk),2,nf))/(c(i,jp,indxc(k))+c(i,j,indxc(k)))
               enddo
            endif
         enddo
      enddo
#endif            
c*$* assert do (concurrent)

c$omp parallel default(private) 
c$omp& shared(flx,vt,indxr,indxc,res,k,lk)
c$omp do
      do j=1,ng2
         do nf=1,nf5
            do i=1,ng1
               jp=mod(j,ng2)+1
               jm=mod(j-2+ng2,ng2)+1
               res(i,j,nf)=res(i,j,nf)
     &              +(c(i,jp,indxc(k))+c(i,j,indxc(k)))/2*flx(i,j,nf)
     &              -(c(i,j,indxc(k))+c(i,jm,indxc(k)))/2*flx(i,jm,nf)
            enddo
         enddo
      enddo
c$omp enddo nowait


      nf=1
c$dir prefer_parallel_ext            
cc$doacross local(i,j)

c$omp do
      do j=1,ng2
         do i=1,ng1
            jp=mod(j,ng2)+1
            rhoflux(i,j,k,2)=rhoflux(i,j,k,2)
     &             +(c(i,jp,indxc(k))+c(i,j,indxc(k)))/2*flx(i,j,nf)
         enddo
      enddo
c$omp end parallel
      return
      end
c endif exactenergy
#endif


c**********************************************************************
      subroutine wphalfx(flx,u,v,nfluid)
c the result flx is one half cell offset to the right.  So
c flx(1)<->u(1.5)
      implicit none
      include 'relaxgl.fi'
      integer nfluid
      real u(ng1,ng2,nfluid),v(ng1,ng2,nfluid),flx(ng1,ng2,nfluid)
cmydist u(*,block,*),v(*,block,*),flx(*,block,*)
c locals
      integer i,j,ip,im,ipp,nf
      real a,b,ap,bp,da,dam,db,dbm, minmod,mm1,mm2


c minmod is a real misnomer, since it is only minmod if wslopelim=1 or wbee=0.
c
      mm1(a,b)=min(wslopelim*a,wbee*b+(1-wbee)*a)
      mm2(a,b)=mm1(min(a,b),max(a,b))
      minmod(a,b)=(sign(1.,a)+sign(1.,b))*mm2(abs(a),abs(b))/2


c$omp parallel default(private) shared(u,v,flx,nfluid)
cfpp$ skip
cc*$* assert do(serial)      
      do nf=1,nfluid
c$dir prefer_parallel_ext            
cc$doacross local(i,j)
c*$* assert do (concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               ip=mod(i,ng1)+1
               ipp=mod(ip,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               a=(u(i,j,nf)+v(i,j,nf))/2
               b=(u(ip,j,nf)-v(ip,j,nf))/2
               da=(u(ip,j,nf)+v(ip,j,nf))/2-a
               dam=a-(u(im,j,nf)+v(im,j,nf))/2
               ap=a+minmod(dam,da)/2
               db=(u(ipp,j,nf)-v(ipp,j,nf))/2-b
               dbm=b-(u(i,j,nf)-v(i,j,nf))/2
               bp=b-minmod(db,dbm)/2
               flx(i,j,nf)=ap-bp
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end

c**********************************************************************
      subroutine wphalfy(flx,u,v,nfluid)
c the result flx is one half cell offset to the right.  So
c flx(1)<->u(1.5)
      implicit none
      include 'relaxgl.fi'
      integer nfluid
      real u(ng1,ng2,nfluid),v(ng1,ng2,nfluid),flx(ng1,ng2,nfluid)
cmydist u(*,block,*),v(*,block,*),flx(*,block,*)

c locals
      integer i,j,jp,jm,jpp,nf
      real a,b,ap,bp,da,dam,db,dbm,minmod,mm1,mm2

      mm1(a,b)=min(wslopelim*a,wbee*b+(1-wbee)*a)
      mm2(a,b)=mm1(min(a,b),max(a,b))
      minmod(a,b)=(sign(1.,a)+sign(1.,b))*mm2(abs(a),abs(b))/2


c$omp parallel default(private) shared(nfluid,flx,u,v)
cfpp$ skip
cc*$* assert do(serial)      
      do nf=1,nfluid
c$dir prefer_parallel_ext            
cc$doacross local(i,j)
c*$* assert do (concurrent)
c$omp do
         do j=1,ng2
            jp=mod(j,ng2)+1
            jpp=mod(jp,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            do i=1,ng1
               a=(u(i,j,nf)+v(i,j,nf))/2
               b=(u(i,jp,nf)-v(i,jp,nf))/2
               da=(u(i,jp,nf)+v(i,jp,nf))/2-a
               dam=a-(u(i,jm,nf)+v(i,jm,nf))/2
               ap=a+minmod(dam,da)/2
               db=(u(i,jpp,nf)-v(i,jpp,nf))/2-b
               dbm=b-(u(i,j,nf)-v(i,j,nf))/2
               bp=b-minmod(db,dbm)/2
               flx(i,j,nf)=ap-bp
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end



c**********************************************************************
      subroutine wphalfz(flz,u,v,nfluid)
c the result flx is one half cell offset to the right.  So
c flx(1)<->u(1.5)
      implicit none
      include 'relaxgl.fi'
      integer nfluid
      real u(5,nfluid,ng1,ng2),v(5,nfluid,ng1,ng2),flz(2,ng1,ng2,nfluid)
cmydist u(*,*,*,block),v(*,*,*,block),flz(*,*,block,*)

c locals
      integer i,j,nf,k1,k,km,kp,kpp
      real a,b,ap,bp,da,dam,db,dbm,minmod,mm1,mm2

      mm1(a,b)=min(wslopelim*a,wbee*b+(1-wbee)*a)
      mm2(a,b)=mm1(min(a,b),max(a,b))
      minmod(a,b)=(sign(1.,a)+sign(1.,b))*mm2(abs(a),abs(b))/2

c$omp parallel do default(private) shared(u,v,flz)
c*$* assert do (concurrent)
      do j=1,ng2
         do i=1,ng1
            nf=1
            k1=1
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
            k1=2
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
            nf=2
            k1=1
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
            k1=2
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
            nf=3
            k1=1
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
            k1=2
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
            nf=4
            k1=1
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
            k1=2
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
            nf=5
            k1=1
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
            k1=2
            k=k1+1
            km=k1
            kp=k1+2
            kpp=k1+3
                  a=(u(k,nf,i,j)+v(k,nf,i,j))/2
                  b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
                  da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
                  dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
                  ap=a+minmod(dam,da)/2
                  db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
                  dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
                  bp=b-minmod(db,dbm)/2
                  flz(k1,i,j,nf)=ap-bp
         enddo
      enddo
      return
      end



c**********************************************************************
      subroutine stepu(unew, uold, res, dt, indxunew,  ngul3,
     &      nfluidcmp, k)
c
      implicit none
      include 'relaxgl.fi'
      integer k, indxunew(ng3), ngul3, nfluidcmp
      real dt, unew(nfluidcmp,ng1,ng2,ngul3), 
     &      uold(nfluidcmp,ng1,ng2,ng3), res(ng1,ng2,nfluidcmp)
cmydist unew(*,*,block,*),uold(*,*,block,*)
cmydist res(*,block,*)
c
c locals
      integer i,j,nf

cfpp$ skip
cc*$* assert do(serial)  
!cdir loopcnt=5
c$omp parallel private(i,j) firstprivate(nfluidcmp)
c$omp& shared(unew,indxunew,uold,k,res,dt)   
      do nf=1,nfluidcmp
c$dir prefer_parallel_ext            
cc$doacross local(i,j)
c*$* assert do (concurrent)
c$omp do
         do j=1,ng2
!cdir vector
            do i=1,ng1
               unew(nf,i,j,indxunew(k))=uold(nf,i,j,k)
     &                -res(i,j,nf)*dt
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end



c**********************************************************************
      subroutine stepuu( u, res, dt, nfluidcmp, k)
c
      implicit none
      include 'relaxgl.fi'
      integer k, nfluidcmp
      real dt, u(nfluidcmp,ng1,ng2,ng3), res(ng1,ng2,nfluidcmp)
cmydist u(*,*,block,*), res(*,block,*)
c locals
      integer i,j,nf

cfpp$ skip
cc*$* assert do(serial)
c$omp parallel default(private) shared(nfluidcmp,u,res,k,dt)
!cdir loopcnt=5
      do nf=1,nfluidcmp
c$dir prefer_parallel_ext
c*$* assert do (concurrent)
c$omp do
          do j=1,ng2
!cdir vector
            do i=1,ng1
               u(nf,i,j,k)=u(nf,i,j,k)-res(i,j,nf)*dt
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end   
      


c**********************************************************************
      subroutine calcbaj(baj,def,indx,k)
      implicit none
      include 'relaxgl.fi'
      include 'gmetric.fi'
      integer k,indx(ng3)
      real baj(7,ng1,ng2,nstore),def(ng1,ng2,ng3)
cmydist baj(*,*,block,*),def(*,block,*)
c
c locals
      integer i,j,ip,jp,im,jm,kp,km
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &      , a23, a33, det, b11,b12,b13,b22,b23,b33


c$omp parallel do default(private) shared(baj,k,indx,def)
c$dir prefer_parallel_ext            
cc$doacross local(i,j)
      do j=1,ng2
         do i=1,ng1
#ifdef GMETRIC
            baj(1,i,j,indx(k))=gbaj(1,i,j,k)
            baj(2,i,j,indx(k))=gbaj(2,i,j,k)
            baj(3,i,j,indx(k))=gbaj(3,i,j,k)
            baj(4,i,j,indx(k))=gbaj(4,i,j,k)
            baj(5,i,j,indx(k))=gbaj(5,i,j,k)
            baj(6,i,j,indx(k))=gbaj(6,i,j,k)
            baj(7,i,j,indx(k))=gbaj(7,i,j,k)
#else
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
     &           -def(ip,jm,k)+def(im,jm,k))/4
            phiyz=(def(i,jp,kp)-def(i,jp,km)
     &           -def(i,jm,kp)+def(i,jm,km))/4
            phixz=(def(ip,j,kp)-def(im,j,kp)
     &           -def(ip,j,km)+def(im,j,km))/4
            a11=1+phixx
            a12=phixy
            a13=phixz
            a22=1+phiyy
            a23=phiyz
            a33=1+phizz
            det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &              -a12**2*a33
            baj(7,i,j,indx(k))=det
            b11=(a22*a33-a23**2)/det
            b12=(a13*a23-a12*a33)/det
            b13=(a12*a23-a13*a22)/det
            b22=(a11*a33-a13**2)/det
            b23=(a12*a13-a11*a23)/det
            b33=(a11*a22-a12**2)/det
            baj(1,i,j,indx(k))=b11
            baj(2,i,j,indx(k))=b12
            baj(3,i,j,indx(k))=b13
            baj(4,i,j,indx(k))=b22
            baj(5,i,j,indx(k))=b23
            baj(6,i,j,indx(k))=b33
#endif
         enddo
      enddo
      return
      end


#ifdef GMETRIC
      subroutine gcalcbaj(def)
      implicit none
      include 'relaxgl.fi'
      include 'gmetric.fi'
      real def(ng1,ng2,ng3)
c
c locals
      integer i,j,ip,jp,im,jm,kp,km,k
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &      , a23, a33, det, b11,b12,b13,b22,b23,b33



c$dir prefer_parallel_ext
c$omp parallel default(private) shared(def,gbaj)
      do k=1,ng3
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
     &           -def(ip,jm,k)+def(im,jm,k))/4
            phiyz=(def(i,jp,kp)-def(i,jp,km)
     &           -def(i,jm,kp)+def(i,jm,km))/4
            phixz=(def(ip,j,kp)-def(im,j,kp)
     &           -def(ip,j,km)+def(im,j,km))/4
            a11=1+phixx
            a12=phixy
            a13=phixz
            a22=1+phiyy
            a23=phiyz
            a33=1+phizz
            det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &              -a12**2*a33
            gbaj(7,i,j,k)=det
            b11=(a22*a33-a23**2)/det
            b12=(a13*a23-a12*a33)/det
            b13=(a12*a23-a13*a22)/det
            b22=(a11*a33-a13**2)/det
            b23=(a12*a13-a11*a23)/det
            b33=(a11*a22-a12**2)/det
            b23=(a12*a13-a11*a23)/det
            b33=(a11*a22-a12**2)/det
            gbaj(1,i,j,k)=b11
            gbaj(2,i,j,k)=b12
            gbaj(3,i,j,k)=b13
            gbaj(4,i,j,k)=b22
            gbaj(5,i,j,k)=b23
            gbaj(6,i,j,k)=b33
         enddo
      enddo
      enddo
c$omp end parallel
      return
      end
#endif
c**********************************************************************
      subroutine flux(vt,defp,u,baj,indxb,indxu,nfluidcmp,ngu3,k)
      implicit none
      include 'relaxgl.fi'
      integer k, nfluidcmp, indxb(ng3),indxu(ng3),ngu3
      real vt(ng1,ng2,nstore,3,nfluidcmp),defp(ng1,ng2,ng3)
     &      ,baj(7,ng1,ng2,nstore),u(nfluidcmp,ng1,ng2,ngu3)
     &      , dascale
      external dascale
cmydist vt(*,block,*,*,*),defp(*,block,*)
cmydist baj(*,*,block,*),u(*,*,block,*)
c
c     locals
      real rootg(ng1,ng2),vt1(3,maxfluidcmp,ng1,ng2)
     &     ,ul(maxfluidcmp,ng1,ng2),rg,pdxg,pdyg,pdzg
     &     ,b11,b12,b13,b22,b23,b33,vtt1,vtt2,vtt3
     &     ,pdcache(29,ng1,ng2),t
cmydist rootg(*,block),vt1(*,*,*,block),ul(*,*,block)
cmydist pdcache(*,*,block)
      integer i,j,ndim,nf,ip,jp,kp,im,jm,km,l, ik
      integer lastk
      save lastk
      data lastk /0/
#ifdef _SX5
      if (nfluidcmp .ne. 5) then
         write(*,*) 'NEC/flux: nfluidcmp != 5',nfluidcmp
         stop
      endif
#endif
#ifdef _ALPHANNN
cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(rootg,vt1,ul,pdcache)
#endif


      ik=indxb(k)
c copy the array into a temporary for efficiency and portability
c across distributed memory machines and f77 compilers.
c 7.2 f77
c zpj: I made c$doacros comments
cc$doacross local(i,nf)

c$omp parallel do default(private) 
c$omp& shared(rootg,baj,ul,u,indxu,ik,k,nfluidcmp)
      do j=1,ng2
!cdir vector
         do i=1,ng1
            rootg(i,j)=baj(7,i,j,ik)
!cdir loopcnt=5 nounroll assert(nfluidcmp=5)
            do nf=1,nfluidcmp
               ul(nf,i,j)=u(nf,i,j,indxu(k))
            enddo
         enddo
      enddo 
      call fluxeuler(vt1,ul,rootg,nfluidcmp,k)
      
c 7.2 f77
c zpj: make  c$doacros comments
cc$doacross local(i,ip,im,jp,jm,kp,km,rg,b11,b12,b13,b22,b23,b33,pdxg
cc$&     ,pdyg,pdzg,vtt1,vtt2,vtt3,l)
c$omp parallel do default(private) 
c$omp& shared(rootg,defp,baj,vt1,u,indxu,vt,nfluidcmp,ik,k)
!cdir nodep
      do j=1,ng2
!cdir vector, nodep
         do i=1,ng1
            ip=mod(i,ng1)+1
            im=mod(i+ng1-2,ng1)+1
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            kp=mod(k,ng3)+1
            km=mod(k+ng3-2,ng3)+1
            rg=rootg(i,j)
            pdxg=(defp(ip,j,k)-defp(im,j,k))/2/rg
            pdyg=(defp(i,jp,k)-defp(i,jm,k))/2/rg
            pdzg=(defp(i,j,kp)-defp(i,j,km))/2/rg
            b11=baj(1,i,j,ik)*rg
            b12=baj(2,i,j,ik)*rg
            b13=baj(3,i,j,ik)*rg
            b22=baj(4,i,j,ik)*rg
            b23=baj(5,i,j,ik)*rg
            b33=baj(6,i,j,ik)*rg
!cdir loopcnt=5
            do l=1,nfluidcmp
               vtt1=vt1(1,l,i,j)-u(l,i,j,indxu(k))*pdxg
               vtt2=vt1(2,l,i,j)-u(l,i,j,indxu(k))*pdyg
               vtt3=vt1(3,l,i,j)-u(l,i,j,indxu(k))*pdzg
               vt(i,j,ik,1,l)=b11*vtt1+b12*vtt2+b13*vtt3
               vt(i,j,ik,2,l)=b12*vtt1+b22*vtt2+b23*vtt3
               vt(i,j,ik,3,l)=b13*vtt1+b23*vtt2+b33*vtt3
            enddo
         enddo
      enddo
      if (k .gt. lastk) then
c         write(*,'(3f20.9)') (vt1(1,1,j,3),j=1,3)
         lastk=k
      endif
      return
      end

c**********************************************************************
      subroutine fluxeuler(vt,u,rootg,nfluidcomp,k)
      implicit none
      include 'relaxgl.fi'
      integer nfluidcomp, k
      real vt(3,nfluidcomp,ng1,ng2), u(nfluidcomp,ng1,ng2)
     &        ,rootg(ng1,ng2)
cmydist vt(*,*,*,block),u(*,*,block),rootg(*,block)
c
c locals
      integer i,j
      real rho,mx,my,mz,vx,vy,vz,engy,kinetic,pressure,econv, gamma, rg
      parameter(gamma=5./3.)
#ifdef COLD
      include 'cold.fi'
#endif      


      if (nfluidcomp .ne. 5) then
         write(*,*) 'fluxeuler: invalid argument: nfluid=',nfluidcomp
         stop
      endif

c$dir prefer_parallel_ext            
cc$doacross local(i,j)
#ifdef COLD
c$omp parallel do default(private) shared(rootg,u,vt,k,cold)
#else
c$omp parallel do default(private) shared(rootg,u,vt,k)
#endif
      do j=1,ng2
         do i=1,ng1
            rg=rootg(i,j)
            rho=u(1,i,j)/rg
            mx=u(2,i,j)/rg
            my=u(3,i,j)/rg
            mz=u(4,i,j)/rg
            vx=mx/rho
            vy=my/rho
            vz=mz/rho
            engy=u(5,i,j)/rg
            kinetic=rho*(vx**2+vy**2+vz**2)/2
            pressure=(gamma-1)*(engy-kinetic)
            pressure=max(0.,pressure)
#ifdef  COLD
            if (cold(i,j,k)) pressure=0
#endif            
            econv=engy+pressure
            vt(1,1,i,j)=mx
            vt(2,1,i,j)=my
            vt(3,1,i,j)=mz
            vt(1,2,i,j)=mx*vx+pressure
            vt(2,2,i,j)=mx*vy
            vt(3,2,i,j)=mx*vz
            vt(1,3,i,j)=my*vx
            vt(2,3,i,j)=my*vy+pressure
            vt(3,3,i,j)=my*vz
            vt(1,4,i,j)=mz*vx
            vt(2,4,i,j)=mz*vy
            vt(3,4,i,j)=mz*vz+pressure
            vt(1,5,i,j)=econv*vx
            vt(2,5,i,j)=econv*vy
            vt(3,5,i,j)=econv*vz
         enddo
      enddo
      return
      end         

c**********************************************************************
      subroutine eulercfl(c,cfluid,u,def,defp)
c C: CFL condition (output)
c cfluid: CFL condition ignoring grid motion (output)
c u,def,defp: as usual      
c
c to save on memory, we calculate the metric from scratch.
c
      implicit none
      include 'relaxgl.fi'
      real c,cfluid,u(5,ng1,ng2,ng3),def(ng1,ng2,ng3),defp(ng1,ng2,ng3)
cmydist u(*,*,block,*),def(*,block,*),defp(*,block,*)
c locals
      integer i,j,k,ip,im,jp,jm,kp,km, negrhocount,imax,jmax,kmax
      real a11,a12,a13,a22,a23,a33,b11,b12,b13,b22,b23,b33,phixx,phiyy
     &     ,phizz,phixy,phiyz,phixz,det,pdx,pdy,pdz,vx,vy,vz,rho
     &     ,kin,ed,pressure,cs,bnorm,tcfl,timax,v11,v12,v13,v21
     &     ,v22,v23,v31,v32,v33,vv11,vv12,vv13,vv21
     &     ,vv22,vv23,vv31,vv32,vv33,rhosum,d11,d22,a,b,am,theta
     &     ,trace,r1,r2,r3,pi,blim1, blim2, blim3, v1, v2, v3
     &     ,vv1, vv2, vv3, rhomax, detmax, csmax, vvmax, vdmax

      real gamma
      parameter(gamma=5./3.)
      parameter(pi=3.14159265358979323846264338328)
#ifdef GMETRIC
      include 'gmetric.fi'
#endif

      c=0
      cfluid=0
      rhosum=0
      negrhocount=0
cATTENTION: this routine does not parallelize under 7.2 f77
cc$doacross local(i,j,k), reduction(c)
c 7.2 f77
#ifdef GMETRIC
c$omp parallel default(private) shared(def,defp,u,gbaj)
c$omp&  reduction(+:rhosum,negrhocount)
c$omp&  reduction(max:c,cfluid,vvmax,vdmax)  
#else
c$omp parallel default(private) shared(def,defp,u)
c$omp&  reduction(+:rhosum,negrhocount)
c$omp&  reduction(max:c,cfluid,vvmax,vdmax) 
#endif   
c*$* assert do(serial)
      do k=1,ng3
c$dir prefer_parallel_ext         
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
#ifdef GMETRIC
               det=gbaj(7,i,j,k)
               b11=gbaj(1,i,j,k)
               b12=gbaj(2,i,j,k)
               b13=gbaj(3,i,j,k)
               b22=gbaj(4,i,j,k)
               b23=gbaj(5,i,j,k)
               b33=gbaj(6,i,j,k)
#else

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
     &                 -a12**2*a33
               b11=(a22*a33-a23**2)/det
               b12=(a13*a23-a12*a33)/det
               b13=(a12*a23-a13*a22)/det
               b22=(a11*a33-a13**2)/det
               b23=(a12*a13-a11*a23)/det
               b33=(a11*a22-a12**2)/det
#endif
               pdx=(defp(ip,j,k)-defp(im,j,k))/2
               pdy=(defp(i,jp,k)-defp(i,jm,k))/2
               pdz=(defp(i,j,kp)-defp(i,j,km))/2
c apply a density limiter while we are at it to prevent CFL blowup
               if ( u(1,i,j,k) .lt. 0.00001) negrhocount=negrhocount+1
               u(1,i,j,k)=max(u(1,i,j,k),0.00001)
               rhosum=rhosum+u(1,i,j,k)
               vx=u(2,i,j,k)/u(1,i,j,k)
               vy=u(3,i,j,k)/u(1,i,j,k)
               vz=u(4,i,j,k)/u(1,i,j,k)
               rho=u(1,i,j,k)/det
               kin=rho*(vx**2+vy**2+vz**2)/2
               ed=u(5,i,j,k)/det-kin
               pressure=(gamma-1)*ed
               cs=sqrt(gamma*abs(pressure/rho))

               v1=vx-pdx
               v2=vy-pdy
               v3=vz-pdz
               vv1=abs(b11*v1+b12*v2+b13*v3)
               vv2=abs(b12*v1+b22*v2+b23*v3)
               vv3=abs(b13*v1+b23*v2+b33*v3)
               blim1=cs*sqrt(b11**2+b12**2+b13**2)
               blim2=cs*sqrt(b12**2+b22**2+b23**2)
               blim3=cs*sqrt(b13**2+b23**2+b33**2)
               tcfl=max(blim1+vv1,blim2+vv2,blim3+vv3)
               cfluid=max(cfluid,blim1,blim2,blim3)
#if (DEBUG>0)
#define EULERDEBUG
#endif
c zpj: define a _OMP to offset the EULERDEBUG
#define _OMP
#if defined(_SGI_SOURCE) || defined(_SX5) || defined(_OMP)
#undef EULERDEBUG
#endif
#ifndef EULERDEBUG
               c=max(c,tcfl)
#else
                if (c .lt. tcfl) then
                  c=tcfl
                  detmax=det
                  rhomax=u(1,i,j,k)
                  imax=i
                  jmax=j
                  kmax=k
                  csmax=cs
                  vvmax=max(abs(v1),abs(v2),abs(v3))
                  vdmax=max(abs(vx),abs(vy),abs(vz))
                endif
#endif
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
#ifdef EULERDEBUG
110   format('eulercfl: detmax, rhomax=',2f8.4,3e9.2,3i4)
      write(*,110)detmax,rhomax,csmax,vvmax,vdmax,imax,jmax,kmax
#endif
      if (negrhocount .gt. 0) then
         write(*,*) 'eulercfl: rho adjusted ',negrhocount,' times'
         rhosum=rhosum/(ng1*ng2*ng3)
cc$doacross local(i,j,k)
c$omp parallel default(private) shared(u,rhosum)
         do k=1,ng3
c$dir prefer_parallel_ext
c*$* assert do(concurrent)
c$omp do
            do j=1,ng2
               do i=1,ng1
                  u(1,i,j,k)=u(1,i,j,k)/rhosum
               enddo
            enddo
c$omp enddo nowait
         enddo
c$omp end parallel
      endif
      return
      end

c**********************************************************************
      subroutine cflkeuler(c,u,baj,defp,indxc,indx,k)
c
c to save on memory, we calculate the metric from scratch.
c
      implicit none
      include 'relaxgl.fi'
      integer k,indx(ng3),indxc(ng3)
      real c(ng1,ng2,ncstore),u(5,ng1,ng2,ng3),baj(7,ng1,ng2,nstore)
     &     ,defp(ng1,ng2,ng3)
cmydist c(*,block,*),u(*,*,block,*),baj(*,*,block,*)
cmydist defp(*,block,*)
c locals
      integer i,j,ip,im,jp,jm,kp,km
      real det,pdx,pdy,pdz,vx,vy,vz,rho,b11,b12,b13,b22,b23,b33
     &      ,kin,ed,pressure,cs,bnorm,tcfl, blim1, blim2, blim3
     &     ,vv1,vv2,vv3, v1, v2, v3
      real gamma
      parameter(gamma=5./3.)

c$omp parallel do default(private) 
c$omp& shared(defp,baj,u,c,indxc,indx,k)
c$dir prefer_parallel_ext         
cc$doacross local(i,j)
      do j=1,ng2
         do i=1,ng1
            ip=mod(i,ng1)+1
            im=mod(i+ng1-2,ng1)+1
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            kp=mod(k,ng3)+1
            km=mod(k+ng3-2,ng3)+1
            pdx=(defp(ip,j,k)-defp(im,j,k))/2
            pdy=(defp(i,jp,k)-defp(i,jm,k))/2
            pdz=(defp(i,j,kp)-defp(i,j,km))/2
            b11=baj(1,i,j,indx(k))
            b12=baj(2,i,j,indx(k))
            b13=baj(3,i,j,indx(k))
            b22=baj(4,i,j,indx(k))
            b23=baj(5,i,j,indx(k))
            b33=baj(6,i,j,indx(k))
            det=baj(7,i,j,indx(k))
            vx=u(2,i,j,k)/u(1,i,j,k)
            vy=u(3,i,j,k)/u(1,i,j,k)
            vz=u(4,i,j,k)/u(1,i,j,k)
            rho=u(1,i,j,k)/det
            kin=rho*(vx**2+vy**2+vz**2)/2
            ed=u(5,i,j,k)/det-kin
            pressure=(gamma-1)*ed
            cs=sqrt(gamma*abs(pressure/rho))
c one should really calculate the largest eigenvalue of the
c dreibein, but I''m too lazy, so will use the infinity norm
c as a bound.
c            bnorm=max(abs(b11)+abs(b12)+abs(b13),abs(b12)
c     &           +abs(b22)+abs(b23),abs(b13),abs(b23),abs(b33))
c            v2=sqrt((vx-pdx)**2+(vy-pdy)**2+(vz-pdz)**2)
c            tcfl=bnorm*(v2+cs)

c this interesting formula is apparently an exact result, as
c listed in Yee and confirmed by gnedin.
            v1=vx-pdx
            v2=vy-pdy
            v3=vz-pdz
            vv1=b11*v1+b12*v2+b13*v3
            vv2=b12*v1+b22*v2+b23*v3
            vv3=b13*v1+b23*v2+b33*v3
            blim1=cs*sqrt(b11**2+b12**2+b13**2)+abs(vv1)
            blim2=cs*sqrt(b12**2+b22**2+b23**2)+abs(vv2)
            blim3=cs*sqrt(b13**2+b23**2+b33**2)+abs(vv3)
            tcfl=max(blim1,blim2,blim3)
            c(i,j,indxc(k))=tcfl
         enddo
      enddo
      return
      end



c*********************************************************************
      subroutine matpp(a,b,s)
      implicit none
      include 'relaxgl.fi'
      real a(ng1,ng2,ng3),b(ng1,ng2,ng3),s
cmydist a(*,block,*),b(*,block,*)
c locals
      integer i,j,k

c$omp parallel private(i,j,k) shared(a,b,s)
cc$doacross local(i,j,k)
      do k=1,ng3
c$dir prefer_parallel_ext
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               a(i,j,k)=a(i,j,k)+b(i,j,k)*s
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end

#if 0

c*********************************************************************
      subroutine fakezerosum(arr,def)
c WORKER routine.  expensive.
c set the mean of arr() to zero, by subtracting a constant multiple
c of \sqrt{g}.
c
      implicit none
      include 'relaxgl.fi'
      real arr(ng1,ng2,ng3), def(ng1,ng2,ng3)
c locals
      integer kp,km,j,jp,jm,i,ip,im,ns1,ns2,ns3,ndim,ib,jb,kb
     &     ,io,ii,jo,ji,ko,ki,kbm,jbm,ibm,k
      real phixx,phiyy,phizz,phixy,phixz,phiyz,a11,a12,a13,a22,a23,a33
     &     ,b11,b12,b13,b22,b23,b33,det,dsum,asum,dfact

c$omp parallel default(private)
c$omp& shared(def,arr)
      do k=1,ng3
         kp=mod(k,ng3)+1
         km=mod(k+ng3-2,ng3)+1
c$dir prefer_vector   
c$omp do   
         do j=1,ng2
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            do i=1,ng1
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
               det=(a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &              -a12**2*a33)
               arr(i,j,k)=arr(i,j,k)-det
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end


c*********************************************************************
      subroutine sphericalgravity(phi,deltae,rho,u,def,nfluid,dts)
c fake gravity by using a spherically symmetric solver:  We simply
c put all the mass on a shell centered at the center of the grid
c and solve assuming it were spherically distributed.      
      implicit none
      include 'relaxgl.fi'
      integer nfluid
      real phi(ng1,ng2,ng3),rho(ng1,ng2,ng3),u(nfluid,ng1,ng2,ng3)
     &     ,def(ng1,ng2,ng3),deltae(ng1,ng2,ng3)
      real dts
c locals
      integer nrad
      parameter(nrad=100000)
      real fpig
      parameter(fpig=2./3.)
      integer i,j,k,ip,im,jp,jm,kp,km,kg,ir
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &      , a23, a33, det, b11,b12,b13,b22,b23,b33,phix,phiy,phiz
     &      ,phi1,phi2,phi3,du,de,ko,koo,p1p,p1m,p2p,p2m,p3p,p3m
      real massradial(0:nrad), phiradial(0:nrad),force
      real cumass,r,x,y,z,newtong,cx,cy,cz,dr,dphi
      logical firsttime
      save firsttime
      data firsttime /.true./
c we store the potential at nrad intervals out to ng3.      


      if (firsttime) then
          firsttime=.false.
          write(*,*) 'Using spherically symmetric gravity solver'
      endif

      cx=ng1/2+1
      cy=ng2/2+1
      cz=ng3/2+1
c$omp parallel default(private) 
c$omp& shared(rho,u)
      do k=1,ng3
c$dir prefer_parallel_ext   
c$omp do      
         do j=1,ng2
            do i=1,ng1
               rho(i,j,k)=u(1,i,j,k)
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      call fakezerosum(rho,def)
      do i=0,nrad
         massradial(i)=0
         phiradial(i)=0
      enddo
c$omp parallel default(private) firstprivate(nrad)
c$omp& shared(def,massradial,rho)
      do k=1,ng3
c$dir prefer_parallel_ext  
c$omp do       
         do j=1,ng2
            do i=1,ng1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               x=i+(def(ip,j,k)-def(im,j,k))/2
               y=j+(def(i,jp,k)-def(i,jm,k))/2
               z=k+(def(i,j,kp)-def(i,j,km))/2
               r=sqrt((x-cx)**2+(y-cy)**2+(z-cz)**2)
               ir=2.*r*nrad/ng3+0.5
               if (ir .le. nrad) then
                  massradial(ir)=massradial(ir)+rho(i,j,k)
               endif
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      cumass=0
      do i=0,nrad/200
        cumass=cumass+massradial(i)
      enddo
      do i=0,nrad/200
        massradial(i)=cumass/(nrad/200.+1)
      enddo
      
      cumass=massradial(0)
      newtong=fpig/4/3.14159265
      phiradial(0)=0
      do i=1,nrad
         r=i*ng3*0.5/nrad
         dr=0.5*ng3/nrad
         dphi=newtong*cumass/(r-0.5*dr)**2
         phiradial(i)=phiradial(i-1)+dphi*dr
         cumass=cumass+massradial(i)
         if (mod(i,nrad/50) .eq. 0) write(15,*) r,cumass
      enddo
c      close(15)
      do i=0,nrad
         phiradial(i)=phiradial(i)-phiradial(nrad)
         r=i*ng3*0.5/nrad
c         write(15,*) r,phiradial(i)
      enddo
c$omp parallel default(private) firstprivate(nrad)
c$omp& shared(def,phi,phiradial)
      do k=1,ng3
c$dir prefer_parallel_ext 
c$omp do        
         do j=1,ng2
            do i=1,ng1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               x=i+(def(ip,j,k)-def(im,j,k))/2
               y=j+(def(i,jp,k)-def(i,jm,k))/2
               z=k+(def(i,j,kp)-def(i,j,km))/2
               r=sqrt((x-cx)**2+(y-cy)**2+(z-cz)**2)
               ir=2.*r*nrad/ng3+0.5
               if (ir .le. nrad) then
                  phi(i,j,k)=phiradial(ir)
               else
                  phi(i,j,k)=phiradial(nrad)
               endif
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      do i=2,ng1-1
c      write(17,*) i+(def(i+1,17,17)-def(i-1,17,17))/2,phi(i,17,17)
      enddo
            
c$omp parallel default(private) firstprivate(nrad)
c$omp& shared(def,phiradial,u,dts,deltae)
      do k=1,ng3
c$dir prefer_parallel_ext         
c$doacross local(i,j)
c 7.2
cc*$* assert no recurrence(u, deltae)         
c$omp do
         do j=1,ng2
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               if (.true.) then
                  x=i+(def(ip,j,k)-def(im,j,k))/2
                  y=j+(def(i,jp,k)-def(i,jm,k))/2
                  z=k+(def(i,j,kp)-def(i,j,km))/2
                  r=sqrt((x-cx)**2+(y-cy)**2+(z-cz)**2)
                  ir=2.*r*nrad/ng3+0.5
                  if (ir .lt. nrad .and. ir .ge. 1) then
         	     dr=0.5*ng3/nrad
                     force=-(phiradial(ir+1)-phiradial(ir-1))/dr/2
                     du=u(1,i,j,k)*force*dts*(x-cx)/r
                     de=(2*u(2,i,j,k)*du+du**2)/(2*u(1,i,j,k))
                     u(2,i,j,k)=u(2,i,j,k)+du
                     du=u(1,i,j,k)*force*dts*(y-cy)/r
                     de=de+(2*u(3,i,j,k)*du+du**2)/(2*u(1,i,j,k))
                     u(3,i,j,k)=u(3,i,j,k)+du
                     du=u(1,i,j,k)*force*dts*(z-cz)/r
                     de=de+(2*u(4,i,j,k)*du+du**2)/(2*u(1,i,j,k))
                     u(4,i,j,k)=u(4,i,j,k)+du
                     u(5,i,j,k)=u(5,i,j,k)+de
                     deltae(i,j,k)=de
                  endif
	       else
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
     &                 -a12**2*a33
               b11=(a22*a33-a23**2)/det
               b12=(a13*a23-a12*a33)/det
               b13=(a12*a23-a13*a22)/det
               b22=(a11*a33-a13**2)/det
               b23=(a12*a13-a11*a23)/det
               b33=(a11*a22-a12**2)/det
               phix=(phi(ip,j,k)-phi(im,j,k))/2
               phiy=(phi(i,jp,k)-phi(i,jm,k))/2
               phiz=(phi(i,j,kp)-phi(i,j,km))/2
               phi1=b11*phix+b12*phiy+b13*phiz
               phi2=b12*phix+b22*phiy+b23*phiz
               phi3=b13*phix+b23*phiy+b33*phiz
               du=-u(1,i,j,k)*phi1*dts
               de=(2*u(2,i,j,k)*du+du**2)/(2*u(1,i,j,k))
               u(2,i,j,k)=u(2,i,j,k)+du
               du=-u(1,i,j,k)*phi2*dts
               de=de+(2*u(3,i,j,k)*du+du**2)/(2*u(1,i,j,k))
               u(3,i,j,k)=u(3,i,j,k)+du
               du=-u(1,i,j,k)*phi3*dts
               de=de+(2*u(4,i,j,k)*du+du**2)/(2*u(1,i,j,k))
               u(4,i,j,k)=u(4,i,j,k)+du
               u(5,i,j,k)=u(5,i,j,k)+de
               deltae(i,j,k)=de
               endif
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end

c*********************************************************************
      subroutine gravity(phi,deltae,rho,u,def,nfluid,dt)
c dts is the effective 2*dt*a^2 in units of \tau      
      implicit none
      include 'relaxgl.fi'
      integer nfluid
      real phi(ng1,ng2,ng3),rho(ng1,ng2,ng3),u(nfluid,ng1,ng2,ng3)
     &     ,def(ng1,ng2,ng3),deltae(ng1,ng2,ng3)
      real dt
c locals
      integer i,j,k,ip,im,jp,jm,kp,km,kg,kgp,kgm
      real gbaj(ng1,ng2,7),fx1(ng1,ng2),fx2(ng1,ng2),fx3(ng1,ng2)
     &     ,fy1(ng1,ng2),fy2(ng1,ng2),fy3(ng1,ng2),fz1(ng1,ng2)
     &     ,fz2(ng1,ng2),fz3(ng1,ng2)
      real gb11mhx,gb12mhx,gb13mhx,gb22mhx,gb23mhx,gb33mhx,vxmhx,vzmhx
     &     ,vymhx,vymhy,vyphz,dts,v1,v2,v3,detmhx,detmhy,detphz,det
     &     ,v2mhx,fx1mhx,fx2mhx,fx3mhx,gb11mhy,gb12mhy,gb13mhy,gb22mhy
     &     ,gb23mhy,gb33mhy,vxmhy,vzmhy,v2mhy,fy1mhy,fy2mhy,fy3mhy
     &     ,phixx,phiyy,phizz,phixy,phixz
     &     ,phiyz,a11,a12,a13,a23,a33,gb11,gb12,gb13,gb23,gb33,gb22,a22
     &     ,gb11phz,gb12phz,gb13phz,gb22phz,gb23phz,gb33phz,vxphz,vzphz
     &     ,v2phz,fx1phz,fx2phz,fx3phz,dux,duy,duz,de,du,fz1n,fz2n,fz3n
      real fpig
      parameter(fpig=2./3.)
      


      dts=dt/fpig
c$omp parallel default(private)
c$omp& shared(rho,u)
ccc$doacross local(i,j,k)
cc*$* assert do (concurrent)         
      do k=1,ng3
c$dir prefer_parallel_ext 
c$omp do        
         do j=1,ng2
            do i=1,ng1
               rho(i,j,k)=fpig*u(1,i,j,k)
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel

      call multigrid(phi,rho,def,u,1.,ng1,ng2,ng3,nfluid,1,1)
      do i=1,0
      call multigrid(phi,rho,def,u,1.,ng1,ng2,ng3,nfluid,1,1)
      enddo


c$dir prefer_parallel_ext
c this loop crashes the SGI compiler...      
ccc$doacross local(i,j)
cc*$* assert do (concurrent)
c$omp parallel default(private)
c$omp& shared(def,gbaj,phi,fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3,u,dts,deltae) 
c$omp do        
      do j=1,ng2
         do i=1,ng1
            kg=1
            ip=mod(i,ng1)+1
            im=mod(i+ng1-2,ng1)+1
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            kgp=mod(kg,ng3)+1
            kgm=mod(kg+ng3-2,ng3)+1
               
            phixx=(def(ip,j,kg)-2*def(i,j,kg)+def(im,j,kg))
            phiyy=(def(i,jp,kg)-2*def(i,j,kg)+def(i,jm,kg))
            phizz=(def(i,j,kgp)-2*def(i,j,kg)+def(i,j,kgm))
            phixy=(def(ip,jp,kg)-def(im,jp,kg)
     &           -def(ip,jm,kg)+def(im,jm,kg))/4
            phiyz=(def(i,jp,kgp)-def(i,jp,kgm)
     &           -def(i,jm,kgp)+def(i,jm,kgm))/4
            phixz=(def(ip,j,kgp)-def(im,j,kgp)
     &           -def(ip,j,kgm)+def(im,j,kgm))/4
            a11=1+phixx
            a12=phixy
            a13=phixz
            a22=1+phiyy
            a23=phiyz
            a33=1+phizz
            det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &                 -a12**2*a33
            gb11=(a22*a33-a23**2)
            gb12=(a13*a23-a12*a33)
            gb13=(a12*a23-a13*a22)
            gb22=(a11*a33-a13**2)
            gb23=(a12*a13-a11*a23)
            gb33=(a11*a22-a12**2)
            gbaj(i,j,1)=gb11
            gbaj(i,j,2)=gb12
            gbaj(i,j,3)=gb13
            gbaj(i,j,4)=gb22
            gbaj(i,j,5)=gb23
            gbaj(i,j,6)=gb33
            gbaj(i,j,7)=det
         enddo
      enddo

c the equation we are trying to write is:
c \dot{\sqrt{g} \rho u^i} + \partial_\alpha F^{i\alpha}=0
c where F^{i\alpha}=\sqrt{g} e^\alpha_j G^{ij}
c and G^{ij}=(V,i V,j - V,k V,k \delta^{ij}/2)/(4\pi G)
c           =( e^{i\alpha}a^{j\beta}V_{,\alpha} V_{,\beta}
c                  - g^{\alpha\beta}V_{,\alpha\beta}/2  )/(4\pi G) .
c      
c Note that G is symmetric, but F is not since the product of
c two symmetric matrices is not necessarily symmetric.      
c The approach is to calculate F^{1\alpha}, etc
c$dir prefer_parallel_ext      
ccc$doacross local(i,j)
cc*$* assert do (concurrent)
c$omp do         
      do j=1,ng2
c$dir prefer_vector         
         do i=1,ng1
            kg=ng3
            k=1
            kp=mod(k,ng3)+1
            km=mod(k+ng3-2,ng3)+1
            ip=mod(i,ng1)+1
            im=mod(i+ng1-2,ng1)+1
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            kgp=mod(kg,ng3)+1
            kgm=mod(kg+ng3-2,ng3)+1
               
            phixx=(def(ip,j,kg)-2*def(i,j,kg)+def(im,j,kg))
            phiyy=(def(i,jp,kg)-2*def(i,j,kg)+def(i,jm,kg))
            phizz=(def(i,j,kgp)-2*def(i,j,kg)+def(i,j,kgm))
            phixy=(def(ip,jp,kg)-def(im,jp,kg)
     &           -def(ip,jm,kg)+def(im,jm,kg))/4
            phiyz=(def(i,jp,kgp)-def(i,jp,kgm)
     &           -def(i,jm,kgp)+def(i,jm,kgm))/4
            phixz=(def(ip,j,kgp)-def(im,j,kgp)
     &           -def(ip,j,kgm)+def(im,j,kgm))/4
            a11=1+phixx
            a12=phixy
            a13=phixz
            a22=1+phiyy
            a23=phiyz
            a33=1+phizz
            det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &                 -a12**2*a33
            gb11=(a22*a33-a23**2)
            gb12=(a13*a23-a12*a33)
            gb13=(a12*a23-a13*a22)
            gb22=(a11*a33-a13**2)
            gb23=(a12*a13-a11*a23)
            gb33=(a11*a22-a12**2)
c the label p is a misnomer in this loop, and should be an m
c because it refers to k-1.            
            gb11phz=(gb11+gbaj(i,j,1))/2
            gb12phz=(gb12+gbaj(i,j,2))/2
            gb13phz=(gb13+gbaj(i,j,3))/2
            gb22phz=(gb22+gbaj(i,j,4))/2
            gb23phz=(gb23+gbaj(i,j,5))/2
            gb33phz=(gb33+gbaj(i,j,6))/2
            detphz=(det+gbaj(i,j,7))/2
            vzphz=phi(i,j,k)-phi(i,j,km)
            vyphz=(phi(i,jp,k)-phi(i,jm,k)+phi(i,jp,km)-phi(i,jm,km)
     &                  )/4
            vxphz=(phi(ip,j,k)-phi(im,j,k)+phi(ip,j,km)-phi(im,j,km)
     &                  )/4

            v1=(gb11phz*vxphz+gb12phz*vyphz+gb13phz*vzphz)/detphz
            v2=(gb12phz*vxphz+gb22phz*vyphz+gb23phz*vzphz)/detphz
            v3=(gb13phz*vxphz+gb23phz*vyphz+gb33phz*vzphz)/detphz
            v2phz=(v1**2+v2**2+v3**2)/2
            fx1phz=v1*v3
            fx2phz=v2*v3
            fx3phz=v3**2-v2phz+fpig*(phi(i,j,km)+phi(i,j,k))/2
            fz1n=gb11phz*fx1phz+gb12phz*fx2phz+gb13phz*fx3phz
            fz2n=gb12phz*fx1phz+gb22phz*fx2phz+gb23phz*fx3phz
            fz3n=gb13phz*fx1phz+gb23phz*fx2phz+gb33phz*fx3phz

            fz1(i,j)=fz1n
            fz2(i,j)=fz2n
            fz3(i,j)=fz3n

         enddo
      enddo
      
      do k=1,ng3
c$dir prefer_parallel_ext         
ccc$doacross local(i,j)
cc*$* assert do (concurrent)
c$omp do         
         do j=1,ng2
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               gb11mhx=(gbaj(im,j,1)+gbaj(i,j,1))/2
               gb12mhx=(gbaj(im,j,2)+gbaj(i,j,2))/2
               gb13mhx=(gbaj(im,j,3)+gbaj(i,j,3))/2
               gb22mhx=(gbaj(im,j,4)+gbaj(i,j,4))/2
               gb23mhx=(gbaj(im,j,5)+gbaj(i,j,5))/2
               gb33mhx=(gbaj(im,j,6)+gbaj(i,j,6))/2
               detmhx=(gbaj(im,j,7)+gbaj(i,j,7))/2
               vxmhx=phi(i,j,k)-phi(im,j,k)
               vymhx=(phi(i,jp,k)-phi(i,jm,k)+phi(im,jp,k)-phi(im,jm,k)
     &                  )/4
               vzmhx=(phi(i,j,kp)-phi(i,j,km)+phi(im,j,kp)-phi(im,j,km)
     &                  )/4
               v1=(gb11mhx*vxmhx+gb12mhx*vymhx+gb13mhx*vzmhx)/detmhx
               v2=(gb12mhx*vxmhx+gb22mhx*vymhx+gb23mhx*vzmhx)/detmhx
               v3=(gb13mhx*vxmhx+gb23mhx*vymhx+gb33mhx*vzmhx)/detmhx
               v2mhx=(v1**2+v2**2+v3**2)/2
               fx1mhx=v1**2-v2mhx+fpig*(phi(im,j,k)+phi(i,j,k))/2
               fx2mhx=v1*v2
               fx3mhx=v1*v3
c fx1mhx=G^{xx}, fx2mhx=G^{xy}
               fx1(i,j)=gb11mhx*fx1mhx+gb12mhx*fx2mhx+gb13mhx*fx3mhx
               fx2(i,j)=gb12mhx*fx1mhx+gb22mhx*fx2mhx+gb23mhx*fx3mhx
               fx3(i,j)=gb13mhx*fx1mhx+gb23mhx*fx2mhx+gb33mhx*fx3mhx
c fx1=F^{xx}, fx2=F^{xy}      
               gb11mhy=(gbaj(i,jm,1)+gbaj(i,j,1))/2
               gb12mhy=(gbaj(i,jm,2)+gbaj(i,j,2))/2
               gb13mhy=(gbaj(i,jm,3)+gbaj(i,j,3))/2
               gb22mhy=(gbaj(i,jm,4)+gbaj(i,j,4))/2
               gb23mhy=(gbaj(i,jm,5)+gbaj(i,j,5))/2
               gb33mhy=(gbaj(i,jm,6)+gbaj(i,j,6))/2
               detmhy=(gbaj(i,jm,7)+gbaj(i,j,7))/2
               vxmhy=(phi(ip,j,k)-phi(im,j,k)+phi(ip,jm,k)-phi(im,jm,k)
     &                  )/4
               vymhy=phi(i,j,k)-phi(i,jm,k)
               vzmhy=(phi(i,j,kp)-phi(i,j,km)+phi(i,jm,kp)-phi(i,jm,km)
     &                  )/4
               v1=(gb11mhy*vxmhy+gb12mhy*vymhy+gb13mhy*vzmhy)/detmhy
               v2=(gb12mhy*vxmhy+gb22mhy*vymhy+gb23mhy*vzmhy)/detmhy
               v3=(gb13mhy*vxmhy+gb23mhy*vymhy+gb33mhy*vzmhy)/detmhy
               v2mhy=(v1**2+v2**2+v3**2)/2
               fy1mhy=v1*v2
               fy2mhy=v2**2-v2mhy+fpig*(phi(i,jm,k)+phi(i,j,k))/2
               fy3mhy=v2*v3
               fy1(i,j)=gb11mhy*fy1mhy+gb12mhy*fy2mhy+gb13mhy*fy3mhy
               fy2(i,j)=gb12mhy*fy1mhy+gb22mhy*fy2mhy+gb23mhy*fy3mhy
               fy3(i,j)=gb13mhy*fy1mhy+gb23mhy*fy2mhy+gb33mhy*fy3mhy
            enddo
         enddo
c$dir prefer_parallel_ext         
ccc$doacross local(i,j)
cc*$* assert do (concurrent) 
c$omp do        
         do j=1,ng2
            do i=1,ng1
               kg=mod(k,ng3)+1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kgp=mod(kg,ng3)+1
               kgm=mod(kg+ng3-2,ng3)+1
               
               phixx=(def(ip,j,kg)-2*def(i,j,kg)+def(im,j,kg))
               phiyy=(def(i,jp,kg)-2*def(i,j,kg)+def(i,jm,kg))
               phizz=(def(i,j,kgp)-2*def(i,j,kg)+def(i,j,kgm))
               phixy=(def(ip,jp,kg)-def(im,jp,kg)
     &              -def(ip,jm,kg)+def(im,jm,kg))/4
               phiyz=(def(i,jp,kgp)-def(i,jp,kgm)
     &              -def(i,jm,kgp)+def(i,jm,kgm))/4
               phixz=(def(ip,j,kgp)-def(im,j,kgp)
     &              -def(ip,j,kgm)+def(im,j,kgm))/4
               a11=1+phixx
               a12=phixy
               a13=phixz
               a22=1+phiyy
               a23=phiyz
               a33=1+phizz
               det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &                 -a12**2*a33
               gb11=(a22*a33-a23**2)
               gb12=(a13*a23-a12*a33)
               gb13=(a12*a23-a13*a22)
               gb22=(a11*a33-a13**2)
               gb23=(a12*a13-a11*a23)
               gb33=(a11*a22-a12**2)

c finished building the next level metric.
               gb11phz=(gb11+gbaj(i,j,1))/2
               gb12phz=(gb12+gbaj(i,j,2))/2
               gb13phz=(gb13+gbaj(i,j,3))/2
               gb22phz=(gb22+gbaj(i,j,4))/2
               gb23phz=(gb23+gbaj(i,j,5))/2
               gb33phz=(gb33+gbaj(i,j,6))/2
               detphz=(det+gbaj(i,j,7))/2
               vzphz=phi(i,j,kp)-phi(i,j,k)
               vyphz=(phi(i,jp,k)-phi(i,jm,k)+phi(i,jp,kp)-phi(i,jm,kp)
     &                  )/4
               vxphz=(phi(ip,j,k)-phi(im,j,k)+phi(ip,j,kp)-phi(im,j,kp)
     &                  )/4

               v1=(gb11phz*vxphz+gb12phz*vyphz+gb13phz*vzphz)/detphz
               v2=(gb12phz*vxphz+gb22phz*vyphz+gb23phz*vzphz)/detphz
               v3=(gb13phz*vxphz+gb23phz*vyphz+gb33phz*vzphz)/detphz
               v2phz=(v1**2+v2**2+v3**2)/2
               fx1phz=v1*v3
               fx2phz=v2*v3
               fx3phz=v3**2-v2phz+fpig*(phi(i,j,kp)+phi(i,j,k))/2

c there is a coefficient of \bar{\rho}=1 in front of phi               
               fz1n=gb11phz*fx1phz+gb12phz*fx2phz+gb13phz*fx3phz
               fz2n=gb12phz*fx1phz+gb22phz*fx2phz+gb23phz*fx3phz
               fz3n=gb13phz*fx1phz+gb23phz*fx2phz+gb33phz*fx3phz

               dux=fx1(ip,j)-fx1(i,j)
               duy=fy1(i,jp)-fy1(i,j)
               duz=fz1n-fz1(i,j)
               du=-(dux+duy+duz)*dts
               de=(2*u(2,i,j,k)*du+du**2)/(2*u(1,i,j,k))
               u(2,i,j,k)=u(2,i,j,k)+du
               dux=fx2(ip,j)-fx2(i,j)
               duy=fy2(i,jp)-fy2(i,j)
               duz=fz2n-fz2(i,j)
               du=-(dux+duy+duz)*dts
               de=de+(2*u(3,i,j,k)*du+du**2)/(2*u(1,i,j,k))
               u(3,i,j,k)=u(3,i,j,k)+du
               dux=fx3(ip,j)-fx3(i,j)
               duy=fy3(i,jp)-fy3(i,j)
               duz=fz3n-fz3(i,j)
               du=-(dux+duy+duz)*dts
               de=de+(2*u(4,i,j,k)*du+du**2)/(2*u(1,i,j,k))
               u(4,i,j,k)=u(4,i,j,k)+du
               u(5,i,j,k)=u(5,i,j,k)+de
               deltae(i,j,k)=de

               gbaj(i,j,1)=gb11
               gbaj(i,j,2)=gb12
               gbaj(i,j,3)=gb13
               gbaj(i,j,4)=gb22
               gbaj(i,j,5)=gb23
               gbaj(i,j,6)=gb33
               gbaj(i,j,7)=det
               fz1(i,j)=fz1n
               fz2(i,j)=fz2n
               fz3(i,j)=fz3n
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end

#endif

c*********************************************************************

      subroutine mgravity(phi,deltae,rho,u,def,defp,nfluid,dts,dtgrav,a)
c dts is the effective 2*dt*a^2 in units of \tau      
      implicit none
      include 'relaxgl.fi'
      integer nfluid
      real phi(ng1,ng2,ng3),rho(ng1,ng2,ng3),u(nfluid,ng1,ng2,ng3)
     &     ,def(ng1,ng2,ng3),defp(ng1,ng2,ng3)
#ifdef EXACTENERGY
      real deltae(ng1,ng2,ng3)
#else
      real deltae
#endif
      real dts,dtgrav,a
cmydist phi(*,block,*),rho(*,block,*),u(*,*,block,*)
cmydist def(*,block,*),defp(*,block,*)
c locals
      integer i,j,k,ip,im,jp,jm,kp,km
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &      , a23, a33, det, b11,b12,b13,b22,b23,b33,phix,phiy,phiz
     &      ,phi1,phi2,phi3,du,de,phixx1,phiyy2,phizz3,corr1,corr2,corr3
     &     ,rhox,rhoy,rhoz
      real fpig, dadtau, omegat
      external dadtau
      parameter(fpig=2./3.)
      logical firsttime, secondmom
      parameter(secondmom=.true.)
      save firsttime
      data firsttime /.true./


      if (firsttime) then
          firsttime=.false.
          write(*,*) 'Using non momentum conserving gravity solver'
          if (secondmom) write(*,*) ' with second order momentum'
      endif


#ifdef NBODY
c 7.2 f77
c*$* assert do(serial)
c$omp parallel default(private) shared(rho,u)
      do k=1,ng3
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               rho(i,j,k)=u(1,i,j,k)
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      call pcalcrho(rho,def)
#else
 
c$omp parallel private(i,j,k) shared(rho,u)           
      do k=1,ng3
c$dir prefer_parallel_ext
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               rho(i,j,k)=fpig*u(1,i,j,k)
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      
#endif
      call multigrid(phi,rho,def,u,1.,ng1,ng2,ng3,nfluid,1,1+32)
      do i=1,1
         call multigrid(phi,rho,def,u,1.,ng1,ng2,ng3,nfluid,1,1+32)
      enddo
         
c$omp parallel default(private)
c$omp& shared(def,phi,u,dts,deltae)           
      phixx1=0
      phiyy2=0
      phizz3=0
c 7.2 f77
c*$* assert do(serial)
      do k=1,ng3
         kp=mod(k,ng3)+1
         km=mod(k+ng3-2,ng3)+1
c$dir prefer_parallel_ext         
c*$* assert no recurrence(u, deltae)         
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            do i=1,ng1
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
     &                 -a12**2*a33
               b11=(a22*a33-a23**2)/det
               b12=(a13*a23-a12*a33)/det
               b13=(a12*a23-a13*a22)/det
               b22=(a11*a33-a13**2)/det
               b23=(a12*a13-a11*a23)/det
               b33=(a11*a22-a12**2)/det
               phix=(phi(ip,j,k)-phi(im,j,k))/2
               phiy=(phi(i,jp,k)-phi(i,jm,k))/2
               phiz=(phi(i,j,kp)-phi(i,j,km))/2
               phi1=b11*phix+b12*phiy+b13*phiz
               phi2=b12*phix+b22*phiy+b23*phiz
               phi3=b13*phix+b23*phiy+b33*phiz


               if (secondmom) then
c the second order correction due to diffusion near shocks:               
                  phixx=phi(ip,j,k)-2*phi(i,j,k)+phi(im,j,k)
                  phiyy=phi(i,jp,k)-2*phi(i,j,k)+phi(i,jm,k)
                  phizz=phi(i,j,kp)-2*phi(i,j,k)+phi(i,j,km)
                  phixy=(phi(ip,jp,k)-phi(im,jp,k)
     &                 -phi(ip,jm,k)+phi(im,jm,k))/4
                  phiyz=(phi(i,jp,kp)-phi(i,jp,km)
     &                 -phi(i,jm,kp)+phi(i,jm,km))/4
                  phixz=(phi(ip,j,kp)-phi(im,j,kp)
     &                 -phi(ip,j,km)+phi(im,j,km))/4
c we make the approximation that e^i_j does not change with space
                  rhox=(u(1,ip,j,k)-u(1,im,j,k))/2
                  rhoy=(u(1,i,jp,k)-u(1,i,jm,k))/2
                  rhoz=(u(1,i,j,kp)-u(1,i,j,km))/2
                  corr1=phixx*rhox+phixy*rhoy+phixz*rhoz
                  corr2=phixy*rhox+phiyy*rhoy+phiyz*rhoz
                  corr3=phixz*rhox+phiyz*rhoy+phizz*rhoz
                  phixx1=b11*corr1+b12*corr2+b13*corr3
                  phiyy2=b12*corr1+b22*corr2+b23*corr3
                  phizz3=b13*corr1+b23*corr2+b33*corr3
               endif
               du=-(u(1,i,j,k)*phi1+phixx1)*dts
               de=(2*u(2,i,j,k)*du+du**2)/(2*u(1,i,j,k))
               u(2,i,j,k)=u(2,i,j,k)+du
               du=-(u(1,i,j,k)*phi2+phiyy2)*dts
               de=de+(2*u(3,i,j,k)*du+du**2)/(2*u(1,i,j,k))
               u(3,i,j,k)=u(3,i,j,k)+du
               du=-(u(1,i,j,k)*phi3+phizz3)*dts
               de=de+(2*u(4,i,j,k)*du+du**2)/(2*u(1,i,j,k))
               u(4,i,j,k)=u(4,i,j,k)+du
               u(5,i,j,k)=u(5,i,j,k)+de
#ifdef EXACTENERGY               
               deltae(i,j,k)=de
#endif               
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel

      return
      end


      
c*********************************************************************

      subroutine stepgh(u,phi,def,defp,tmp,tmp2,isopt,a,tdt, dtfluid)
c parameter isopt: 0: regular time step
c                  1: initialize and step
c                  2: reentry after a time sync step
c      
c takes a second order split step, advances t by 2dt
c To call for halfstepping the final Nbodies, call with tdt<0,
c and a=a(t_f), defp unchanged.      
      implicit none
      include 'relaxgl.fi'
      integer nfluid,isopt
      parameter(nfluid=5)
      real u(5,ng1,ng2,ng3),def(ng1,ng2,ng3),defp(ng1,ng2,ng3)
     &     ,phi(ng1,ng2,ng3),tmp(ng1,ng2,ng3),tmp2(ng1+2,ng2,ng3)
cmydist u(*,*,block,*),def(*,block,*),defp(*,block,*)
cmydist phi(*,block,*),tmp(*,block,*),tmp2(*,block,*)
c tmp2 is for use in FFTs, and may be equivalenced with tmp()

      real tdt,a,rho2, dtfluid

c apply operator splitting on hydro and gravity.  In this framework,
c we use units where a=t**2, dx=1

c the prescription is simple:  L_H L_GG L_H, we take two
c gravity steps at once.



c locals:
      real engy1,engy2,engyk,engyp
     &     ,adot,tau, phiold(ng1,ng2,ng3),t
     &     ,engythermal, ethermal1, ethermal2, dtold,ptmp, dt, c
cmydist phiold(*,block,*)
c c is not really used except as redundant dummy argument      
c tmp() is used to store \delta \rho
#ifdef EXACTENERGY
      real rhoflux(ng1,ng2,ng3,3),defpm(ng1,ng2,ng3),defm(ng1,ng2,ng3)
     &     , deltae(ng1,ng2,ng3)
#else
      real deltae
#endif
#ifdef NBODY
      real aold
c the common block is only shared with checkpoint and restart      
      common /defold/ aold
c we will use tmp(,,) to store def(,,) from the last time step.
c Otherwise, it is also used to store rho(,,)
#endif      
      integer i,j,k
c deltae is only the non-conserved part of the energy.
c      
      real dadtau, dascale
      external  dadtau, dascale
      real adotk, adotw, engy0, engy0k, gdt, gdtold
c common block is only shared with checkpoint/restart      
      common /cstepgh/ phiold, adotk, adotw, engy0, engy0k
     &     , gdtold, dtold
      real fpig
      parameter(fpig=2./3.)

#ifdef _ALPHANNN
cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(phiold)
#endif


      dt=tdt
c we need to be a bit careful if the last time step was synched.
c       
      if (gdtold .le. 0 .and. isopt .eq. 0 ) isopt=2
      if (isopt .ge. 1 ) then
         dtold=0
         gdtold=1.e-20
         if (isopt .eq. 1 ) then
            adotk=0
            adotw=0
c we divide by gdtold later, so we cant make it too small
#ifdef NBODY
            call matchhydro(phiold,tmp,u,def)

c$omp parallel default(private) shared(u,tmp)
c*$* assert do(serial)
            do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do

               do j=1,ng2
                  do i=1,ng1
                     tmp(i,j,k)=u(1,i,j,k)
                  enddo
               enddo
c$omp enddo nowait
            enddo
c$omp end parallel
            call pcalcrho(tmp,def)
            call multigrid(phi,tmp,def,u,1.,ng1,ng2,ng3,nfluid,1,1+32)
            do i=1,1
            call multigrid(phi,tmp,def,u,1.,ng1,ng2,ng3,nfluid,1,1+32)
            enddo
            aold=a-dascale(t,a,dt)
c gdt is the time between the center of two steps, which is the
c gravity time step.
         endif
c some careful juggling with the tmp(,,) array:         
c 7.2 f77
c*$* assert do(serial)
c$omp parallel default(private) shared(tmp,def)
            do k=1,ng3
c*$* assert do(concurrent)
               do j=1,ng2
                  do i=1,ng1
                     tmp(i,j,k)=def(i,j,k)
                  enddo
               enddo
            enddo
c$omp end parallel
#else        
         endif
#endif
      endif

      if (tdt .lt. 0) then
         write(*,*) 'final exiting call to stepgh'
         dt=0
      endif



      
c let us recall the unit conventions:
c dx=1, rhomean=1
c 4 \pi G = 2/3
c      
c the first calcdefp call was in the calling routine.
c      call calcdefp(defp,tmp,tmp2,def,u,dt,nfluid)

      gdt=dt+dtold
#ifdef NBODY
      
      call matpp(def,defp,dt)

c     step to current center time step from last time step.
c tmp(,,) contains def(,,) at the last time step 

      call stepxv(phi,tmp,def,gdt,gdtold,aold)
      
c def(,,) contains now the future value of def(,,)

      call matpp(def,defp,-dt)
#endif
      
c we need to store \rho in order to find \dot{\rho} in the middle of the step. 
c*$* assert do(serial)
#ifdef EXACTENERGY
c$omp parallel default(private) firstprivate(isopt,gdt,gdtold)
c$omp& shared(deltae,phi,phiold,rhoflux)
#else
c$omp parallel default(private) firstprivate(isopt,gdt,gdtold)
c$omp& shared(deltae,phi,phiold)
#endif
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
#ifdef EXACTENERGY
               rhoflux(i,j,k,1)=0
               rhoflux(i,j,k,2)=0
               rhoflux(i,j,k,3)=0
               deltae(i,j,k)=0
#endif
               ptmp=phi(i,j,k)
               if (isopt .eq. 0) then
                  phi(i,j,k)=phi(i,j,k)*(1+gdt/gdtold)
     &                 -phiold(i,j,k)*gdt/gdtold
               endif
               phiold(i,j,k) = ptmp
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
#ifdef NBODY      
      
      aold=a
#endif

      if (tdt .lt. 0) then
c we need to assume that in this case, defp is the leftover from
c the last iteration
c*$* assert do(serial)
c$omp parallel private(i,j,k) shared(tmp,u)
         do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
            do j=1,ng2
               do i=1,ng1
#ifdef NBODY                  
                  tmp(i,j,k)=u(1,i,j,k)
#else
                  tmp(i,j,k)=fpig*u(1,i,j,k)
#endif                  
               enddo
            enddo
c$omp enddo nowait
         enddo
c$omp end parallel
#ifdef NBODY         
         call pcalcrho(tmp,def)
#endif         
         call multigrid(phi,tmp,def,u,1.,ng1,ng2,ng3,nfluid,1,1+32)
         do i=1,1
         call multigrid(phi,tmp,def,u,1.,ng1,ng2,ng3,nfluid,1,1+32)
         enddo
#ifdef NBODY         
         gdt=0
         call stepxv(phi,def,def,gdt,dtold,a)
c just half step the velocities and nothing else.         
#endif

c         return
c lets calculate the LI energy at the synchronized time step
c$omp parallel private(i,j,k) shared(defp)         
         do k=1,ng3
c$omp do
            do j=1,ng2
               do i=1,ng1
c                  defp(i,j,k)=0
               enddo
            enddo
c$omp enddo nowait
         enddo
c$omp end parallel
      endif
       
      call matpp(def,defp,dt/2)
#ifdef GMETRIC
      call gcalcbaj(def)
#endif

      gdtold=gdt


c tau should actually be the negative root, and adot=-18/tau**3      

      adot=dadtau(a)

      if (tdt .gt. 0) then
#ifdef EXACTENERGY      
         call rhorelaxing(u,rhoflux,def,defp,c,dt,nfluid)
#else      
         call relaxing(u,def,defp,c,dt,nfluid)
#endif      
         call matpp(def,defp,dt/2)
#ifdef GMETRIC
         call gcalcbaj(def)
#endif
      endif
      call sumenergy(engy1,ethermal1,u)
      if (tdt .gt. 0 ) then
      call mgravity(phi,deltae,tmp,u,def,defp,nfluid,2*a*dt,dt+dtold,a)
      endif
c      call sphericalgravity(phi,deltae,tmp,u,def,nfluid,2*a*dt)
c      call gravity(phi,deltae,tmp,u,def,nfluid,2*a*dt)

c calculate layzer-irvine energy      
      call sumenergy(engy2,ethermal2,u)
      engyk=(engy1+engy2)/2
      engythermal=(ethermal1+ethermal2)/2
      call sumpotential(engyp,phi,tmp,a)
c      call sumgradphipot(engyp,phi,tmp)
      adotk=adotk-adot*engyk*dt*2/a**2
      adotw=adotw+adot*engyp*dt*2
      if (isopt .eq. 1) then
         engy0 = engyk+engyp*a-adotw
         engy0k =  engyk/a+engyp-adotk
      endif
c      
c the conservation law follows from our equations and units:
c the energy source term \dot{e}=-a\rho v \grad V
c and we perform a couple of integration by parts to isolate the
c potential energy a\rho\phi/2 with a total source \dot{a}\rho\phi/2.
c This casts eqn (21.37) of peebles (1994) into
c       d/dt (a^2 K + a^2 W)=a\dot{a}W
c       where a^2 K = engyk,   a W= engyp
c There is also the alternative form (21.34)
c       d/dt (a K + a W ) = -\dot{a} K
c      
      write(*,99) engyk,100*engythermal/(engyk+1.e-10), engyp*a
     &     , engyk+engyp*a-engy0, adotw
c      
 99   format(' KE=',g10.3,' TE=',f6.1,'% PE=',g10.3,' sum=',g10.3
     &     ,'  L-I w=',g10.3)
c
 18   format( '  alternative sum=' , g10.3
     &      , '  L-I adotk=' , g10.3 , '  R=' , g11.5 )
      write(*,18) engyk/a+engyp-engy0k, adotk,
     &		-(engyk/a-engy0k-adotk)/engyp

      

c
#ifdef EXACTENERGY 
c*$* assert do(serial)
c$omp parallel default(private) shared(defpm,defp,def)
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               defpm(i,j,k)=defp(i,j,k)
               defm(i,j,k)=def(i,j,k)
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
#endif
      if (tdt .gt. 0) then
           call calcdefp(defp,tmp,tmp2,def,u,dtfluid,dt,nfluid)
c      	   write(29,'(7e11.3)') 
c     &	a,engyk, engythermal, engyp*a, tmp(1,1,1),tmp(2,1,1),tmp(3,1,1)
      endif
#ifndef EXACTENERGY
c$omp parallel private(i,j,k) shared(defp,tmp,def)
#else 
c$omp parallel private(i,j,k) shared(defp,tmp,def,defpm)
#endif
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
#ifdef EXACTENERGY
               defpm(i,j,k)=(defpm(i,j,k)+defp(i,j,k))/2
#endif
#ifdef NBODY
               tmp(i,j,k)=def(i,j,k)
#endif
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      if (tdt .lt. 0) then
         return
      endif
      call matpp(def,defp,dt/2)
#ifdef GMETRIC
         call gcalcbaj(def)
#endif
#ifdef EXACTENERGY      
      call rhorelaxing(u,rhoflux,def,defp,c,dt,nfluid)
      call rfengymoving(rhoflux,defm,defpm,u,2.)
      call engycomp(u,rhoflux,deltae,phi,nfluid,2*a*dt)
#else      
      call relaxing(u,def,defp,c,dt,nfluid)
#endif
      call matpp(def,defp,dt/2)
#ifdef GMETRIC
         call gcalcbaj(def)
#endif
      dtold=dt
      return
      end


#ifdef EXACTENERGY
      subroutine engycomp(u,rhoflux,deltae,phi,nfluid,dts)
      implicit none
      include 'relaxgl.fi'
      integer nfluid
      real u(nfluid,ng1,ng2,ng3),deltae(ng1,ng2,ng3),phi(ng1,ng2,ng3)
     &            ,rhoflux(ng1,ng2,ng3,3),dts

c locals
      integer i,j,k,ip,im,jp,jm,kp,km
      real de, ecorr, ecorrrms

      if (nfluid .ne. 5) then
         write(*,*)'engycomp: invalid argument nfluid=',nfluid
         stop
      endif
      ecorr=0
      ecorrrms=0
c*$* assert do(serial)
c$omp parallel default(private)
c$omp& shared(deltae,dts,rhoflux,phi,u) reduction(+:ecorr,ecorrrms)
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1

               de=-deltae(i,j,k)
     &           -dts*(rhoflux(i,j,k,1)*(phi(ip,j,k)+phi(i,j,k))/2
     &                 -rhoflux(im,j,k,1)*(phi(i,j,k)+phi(im,j,k))/2
     &                 +rhoflux(i,j,k,2)*(phi(i,jp,k)+phi(i,j,k))/2
     &                  -rhoflux(i,jm,k,2)*(phi(i,j,k)+phi(i,jm,k))/2
     &                 +rhoflux(i,j,k,3)*(phi(i,j,kp)+phi(i,j,k))/2
     &                  -rhoflux(i,j,km,3)*(phi(i,j,k)+phi(i,j,km))/2 )
     &              /2
     &            +dts*phi(i,j,k)*(rhoflux(i,j,k,1)-rhoflux(im,j,k,1)
     &                  +rhoflux(i,j,k,2)-rhoflux(i,jm,k,2)
     &                  +rhoflux(i,j,k,3)-rhoflux(i,j,km,3)  )
     &              /2

c               de=0
c relaxing is called twice, so we divide by 2              
               u(5,i,j,k)=u(5,i,j,k)+de               
               ecorr=ecorr+de
               ecorrrms=ecorrrms+de**2

            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel

c in the absence of flux limiters, these corrections should be small.      
      ecorrrms=ng1*ng2*ng3*sqrt(ecorrrms/(ng1*ng2*ng3))
      write(*,93) ecorr,ecorrrms
 93   format(' energy correction=',g10.3,'  rms=',g10.3)

      return
      end



c*********************************************************************
      subroutine rfengymoving(rhoflux,def,defp,u,scale)
c we assume that the flux vector index k refers to k+1/2
c      implicit none
      include 'relaxgl.fi'
      real rhoflux(ng1,ng2,ng3,3),def(ng1,ng2,ng3),defp(ng1,ng2,ng3)
     &     ,u(5,ng1,ng2,ng3)
c the nfluid doesn't really matter as long as the first four are
c density and momenta      
      real scale

c locals
      integer i,j,k,ip,im,jp,jm,kp,km
      real gbaj(ng1,ng2,6)


c$omp parallel default(private)
c$omp& shared(def,gbaj,defp,rhoflux,u,scale)
c$dir prefer_parallel_ext
c$omp do
      do j=1,ng2
         do i=1,ng1
            kg=1
            ip=mod(i,ng1)+1
            im=mod(i+ng1-2,ng1)+1
            jp=mod(j,ng2)+1
            jm=mod(j+ng2-2,ng2)+1
            kgp=mod(kg,ng3)+1
            kgm=mod(kg+ng3-2,ng3)+1
               
            phixx=(def(ip,j,kg)-2*def(i,j,kg)+def(im,j,kg))
            phiyy=(def(i,jp,kg)-2*def(i,j,kg)+def(i,jm,kg))
            phizz=(def(i,j,kgp)-2*def(i,j,kg)+def(i,j,kgm))
            phixy=(def(ip,jp,kg)-def(im,jp,kg)
     &           -def(ip,jm,kg)+def(im,jm,kg))/4
            phiyz=(def(i,jp,kgp)-def(i,jp,kgm)
     &           -def(i,jm,kgp)+def(i,jm,kgm))/4
            phixz=(def(ip,j,kgp)-def(im,j,kgp)
     &           -def(ip,j,kgm)+def(im,j,kgm))/4
            a11=1+phixx
            a12=phixy
            a13=phixz
            a22=1+phiyy
            a23=phiyz
            a33=1+phizz
            det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &                 -a12**2*a33
            gb11=(a22*a33-a23**2)
            gb12=(a13*a23-a12*a33)
            gb13=(a12*a23-a13*a22)
            gb22=(a11*a33-a13**2)
            gb23=(a12*a13-a11*a23)
            gb33=(a11*a22-a12**2)
            gbaj(i,j,1)=gb11
            gbaj(i,j,2)=gb12
            gbaj(i,j,3)=gb13
            gbaj(i,j,4)=gb22
            gbaj(i,j,5)=gb23
            gbaj(i,j,6)=gb33
         enddo
      enddo

c*$* assert do(serial)
      do k=1,ng3
         kg=mod(k,ng3)+1
c$dir no_recurrence
cc*$* assert no recurrence(gbaj)
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               ip=mod(i,ng1)+1
               im=mod(i+ng1-2,ng1)+1
               jp=mod(j,ng2)+1
               jm=mod(j+ng2-2,ng2)+1
               kp=mod(k,ng3)+1
               km=mod(k+ng3-2,ng3)+1
               
               kgp=mod(kg,ng3)+1
               kgm=mod(kg+ng3-2,ng3)+1


               gb11phx=(gbaj(ip,j,1)+gbaj(i,j,1))/2
               gb12phx=(gbaj(ip,j,2)+gbaj(i,j,2))/2
               gb13phx=(gbaj(ip,j,3)+gbaj(i,j,3))/2
               gb22phx=(gbaj(ip,j,4)+gbaj(i,j,4))/2
               gb23phx=(gbaj(ip,j,5)+gbaj(i,j,5))/2
               gb33phx=(gbaj(ip,j,6)+gbaj(i,j,6))/2
               pvx=defp(ip,j,k)-defp(i,j,k)
               pvy=(defp(i,jp,k)-defp(i,j,k)+defp(ip,jp,k)-defp(ip,j,k)
     &                              )/2
               pvz=(defp(i,j,kp)-defp(i,j,k)+defp(ip,j,kp)-defp(ip,j,k)
     &                              )/2
               p1=gb11phx*pvx+gb12phx*pvy+gb13phx*pvz
               
               rhoflux(i,j,k,1)=rhoflux(i,j,k,1)
     &              +p1*(u(1,i,j,k)+u(1,ip,j,k))*scale/2


               gb11phy=(gbaj(i,jp,1)+gbaj(i,j,1))/2
               gb12phy=(gbaj(i,jp,2)+gbaj(i,j,2))/2
               gb13phy=(gbaj(i,jp,3)+gbaj(i,j,3))/2
               gb22phy=(gbaj(i,jp,4)+gbaj(i,j,4))/2
               gb23phy=(gbaj(i,jp,5)+gbaj(i,j,5))/2
               gb33phy=(gbaj(i,jp,6)+gbaj(i,j,6))/2
               pvx=(defp(ip,j,k)-defp(i,j,k)+defp(ip,jp,k)-defp(i,jp,k)
     &                              )/2
               pvy=defp(i,jp,k)-defp(i,j,k)
               pvz=(defp(i,j,kp)-defp(i,j,k)+defp(i,jp,kp)-defp(i,jp,k)
     &                              )/2
               p2=gb12phy*pvx+gb22phy*pvy+gb23phy*pvz
               
               rhoflux(i,j,k,2)=rhoflux(i,j,k,2)
     &              +p2*(u(1,i,j,k)+u(1,i,jp,k))*scale/2

               phixx=(def(ip,j,kg)-2*def(i,j,kg)+def(im,j,kg))
               phiyy=(def(i,jp,kg)-2*def(i,j,kg)+def(i,jm,kg))
               phizz=(def(i,j,kgp)-2*def(i,j,kg)+def(i,j,kgm))
               phixy=(def(ip,jp,kg)-def(im,jp,kg)
     &              -def(ip,jm,kg)+def(im,jm,kg))/4
               phiyz=(def(i,jp,kgp)-def(i,jp,kgm)
     &              -def(i,jm,kgp)+def(i,jm,kgm))/4
               phixz=(def(ip,j,kgp)-def(im,j,kgp)
     &              -def(ip,j,kgm)+def(im,j,kgm))/4
               a11=1+phixx
               a12=phixy
               a13=phixz
               a22=1+phiyy
               a23=phiyz
               a33=1+phizz
               det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &                 -a12**2*a33
               b11=(a22*a33-a23**2)/det
               b12=(a13*a23-a12*a33)/det
               b13=(a12*a23-a13*a22)/det
               b22=(a11*a33-a13**2)/det
               b23=(a12*a13-a11*a23)/det
               b33=(a11*a22-a12**2)/det

               gb11phz=(b11+gbaj(i,j,1))/2
               gb12phz=(b12+gbaj(i,j,2))/2
               gb13phz=(b13+gbaj(i,j,3))/2
               gb22phz=(b22+gbaj(i,j,4))/2
               gb23phz=(b23+gbaj(i,j,5))/2
               gb33phz=(b33+gbaj(i,j,6))/2
               pvx=(defp(ip,j,k)-defp(i,j,k)+defp(ip,j,kg)-defp(i,j,kg)
     &                              )/2
               pvy=(defp(i,jp,k)-defp(i,j,k)+defp(i,jp,kg)-defp(i,j,kg)
     &                              )/2
               pvz=defp(i,j,kg)-defp(i,j,k)
               p3=gb13phz*pvx+gb23phz*pvy+gb33phz*pvz
               rhoflux(i,j,k,3)=rhoflux(i,j,k,3)
     &              +p3*(u(1,i,j,kg)+u(1,i,j,k))*scale/2
               gbaj(i,j,1)=b11
               gbaj(i,j,2)=b12
               gbaj(i,j,3)=b13
               gbaj(i,j,4)=b22
               gbaj(i,j,5)=b23
               gbaj(i,j,6)=b33
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end



#endif

c*********************************************************************

      subroutine cfldefp(c,def,defp)
      implicit none
      include 'relaxgl.fi'
      real c,def(ng1,ng2,ng3),defp(ng1,ng2,ng3)
cmydist def(*,block,*),defp(*,block,*)
c locals
      integer i,j,k,ip,im,jp,jm,kp,km,ndetn
      real vmax,dfp,v11,v12,v13,v22,v23,v33,vxx,vxy,vxz,vyx,vyy,vyz
     &     ,vzx,vzy,vzz,tmax
      real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22
     &      , a23, a33, det, b11,b12,b13,b22,b23,b33, a, b, d11, d22, d33
     &	    , theta, am, pi,trace
      complex r1, r2, r3, desc, s1, s2
      
      parameter(pi=3.14159265358979323846264338328)
      
      ndetn=0
      vmax=0
cc$doacross local(i,j,k), reduction(vmax)
c*$* assert do(serial)
c$omp parallel default(private)
c$omp& shared(def,defp) reduction(+:ndetn) reduction(max:vmax)
      do k=1,ng3
c$dir force_parallel_ext         
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
c$dir force_vector            
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
     &                 -a12**2*a33
               if (det .lt. 0) then
                  ndetn=ndetn+1
c                  write(*,*)'det<0!'
               endif
               b11=(a22*a33-a23**2)/det
               b12=(a13*a23-a12*a33)/det
               b13=(a12*a23-a13*a22)/det
               b22=(a11*a33-a13**2)/det
               b23=(a12*a13-a11*a23)/det
               b33=(a11*a22-a12**2)/det
               dfp=defp(i,j,k)
               v11=(defp(ip,j,k)-2*dfp+defp(im,j,k))/2
               v22=(defp(i,jp,k)-2*dfp+defp(i,jm,k))/2
               v33=(defp(i,j,kp)-2*dfp+defp(i,j,km))/2
               v12=(defp(ip,jp,k)-defp(im,jp,k)
     &              -defp(ip,jm,k)+defp(im,jm,k))/4
               v23=(defp(i,jp,kp)-defp(i,jp,km)
     &              -defp(i,jm,kp)+defp(i,jm,km))/4
               v13=(defp(ip,j,kp)-defp(im,j,kp)
     &              -defp(ip,j,km)+defp(im,j,km))/4
               vxx=v11*b11+v12*b12+v13*b13
               vxy=v11*b12+v12*b22+v13*b23
               vxz=v11*b13+v12*b23+v13*b33
               vyx=v12*b11+v22*b12+v23*b13
               vyy=v12*b12+v22*b22+v23*b23
               vyz=v12*b13+v22*b23+v23*b33
               vzx=v13*b11+v23*b12+v33*b13
               vzy=v13*b12+v23*b22+v33*b23
               vzz=v13*b13+v23*b23+v33*b33
c     
c now find the eigenvalues using the cubic equation.
c
               trace=vxx+vyy+vzz
               d11=vxx-trace/3
               d22=vyy-trace/3
c               d33=vzz-trace/3
c
c after we take out the trace, the characteristic equation does not
c have a quadratic term.
c
c the coefficient of \lambda:  \lambda^3+a\lambda+b=0             
c  a=d11*d22+d11*d33+d22*d33-(vxy*vyx+vxz*vzx)
               a=-(d11**2+d11*d22+d22**2+vxy*vyx+vxz*vzx+vyz*vzy)
c       b = - DET(V) (the traceless part)
               b= d11*vyz*vzy + d22*vxz*vzx - d11*vxy*vyx - vxy*vyx*d22
     &    + d11**2*d22 +  d11*d22**2 - vxy*vyz*vzx - vyx*vzy*vxz
c
c if the eigenvalues are real, then a < 0.               
c V is the product of symmetric matrices.  I postulate that its eigenvalues
c are in fact real.               
c
c               if (abs(a) .lt. 1.e-20) a=1.e-20
	       desc=(a/3)**3+(b/2)**2
	       s1=(-b/2+sqrt(desc))**(1./3.)
	       s2=(-b/2-sqrt(desc))**(1./3.)
	       r1=s1+s2
 	       r2=-(s1+s2)/2+cmplx(0.,sqrt(3.))/2*(s1-s2)
 	       r3=-(s1+s2)/2-cmplx(0.,sqrt(3.))/2*(s1-s2)
               tmax=max(abs(r1+trace/3),abs(r2+trace/3),abs(r3+trace/3))
c     
c     this temporary is necessary because the CRAY compiler does not
c recognize multiple argument critical loop MAX reductions.
c
               vmax=max(vmax,tmax)
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel

c remember that we take two steps with each dt.  The following constraint
c says that we want to limit the time step to 1/10 the zeldovich time.      
c      c=10*vmax
      c=4*vmax
      if (ndetn .gt. 0) then
         write(*,*)'det<0 ',ndetn,' times'
      endif
      return
      end


c*********************************************************************
      subroutine sumenergy(engy,engythermal,u)
      implicit none
c everything is local
      include 'relaxgl.fi'
      integer nfluidcmp
      parameter(nfluidcmp=5)
      real u(nfluidcmp,ng1,ng2,ng3),engy,engythermal
cmydist u(*,*,block,*)
#ifdef NBODY
      include 'nbody.fi'
      real p1,epart
#endif
c locals
      include 'globalpa.fi'
      integer i,j,k
      real engyk

      engy=0
      engyk=0
c$omp parallel private(i,j,k)
c$omp& shared(u) reduction(+:engy,engyk)
c*$* assert do(serial)
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
               engy=engy+u(5,i,j,k)
               engyk=engyk+(u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2)
     &              /(2*u(1,i,j,k))
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      engythermal=engy-engyk
#ifdef NBODY
      p1=1./npartcell
      epart=0
      do i=1,npmassive
         p1=pmass(i)/npartcell
         epart=epart+(xv(4,i)**2+xv(5,i)**2+xv(6,i)**2)*p1/2
      enddo
      write(*,*)'sumenergy:E(particle)/E(gas)=',epart/max(1.e-10,engy)
      engy=engy*omegab/omega0+(1-omegab/omega0)*epart
#endif      
      return
      end


c*********************************************************************
      subroutine sumerho2(erho2,u,def)
      implicit none
      include 'relaxgl.fi'
      integer nfluidcmp
      parameter(nfluidcmp=5)
      real u(nfluidcmp,ng1,ng2,ng3),erho2
     &		,def(ng1,ng2,ng3)
cmydist u(*,*,block,*),def(*,block,*)
c locals
      real a11,a12,a13,a22,a23,a33,phixx,phiyy,phizz,phixy,phixz,phiyz
     &		,det
      integer i,j,k,ip,im,jp,jm,kp,km

      erho2=0
c$omp parallel default(private)
c$omp& shared(def,u) reduction(+:erho2)
c*$* assert do(serial)
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
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
     &                 -a12**2*a33
               erho2=erho2+u(1,i,j,k)/det* ( u(5,i,j,k) 
     &			-  (u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2)
     &              			/(2*u(1,i,j,k))  )
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end

c*********************************************************************
      subroutine sumpotential(engy,phi,rho,a)
      implicit none
c everything is local
      include 'relaxgl.fi'
      real phi(ng1,ng2,ng3),rho(ng1,ng2,ng3),engy
cmydist phi(*,block,*),rho(*,block,*)
c locals
      real fpig, omegat, dadtau, a
      external dadtau
      parameter(fpig=2./3.)
      integer i,j,k

      engy=0
c$omp parallel private(i,j,k)
c$omp& shared(phi,rho) reduction(+:engy)
c*$* assert do(serial)
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
c rho actually has one 4\pi G too many.               
               engy=engy+phi(i,j,k)*rho(i,j,k)/2/fpig
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end


c*********************************************************************
      subroutine sumgradphipot(engy,phi,rho)
      implicit none
c everything is local
      include 'relaxgl.fi'
      real phi(ng1,ng2,ng3),rho(ng1,ng2,ng3),engy
cmydist phi(*,block,*),rho(*,block,*)
c locals
      real fpig
      parameter(fpig=2./3.)
      integer i,j,k,ip,jp,kp

      engy=0
c$omp parallel default(private) shared(phi) reduction(-:engy)
c*$* assert do(serial)
      do k=1,ng3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,ng2
            do i=1,ng1
            ip=mod(i,ng1)+1
            jp=mod(j,ng2)+1
            kp=mod(k,ng3)+1
               engy=engy-(phi(ip,j,k)-phi(i,j,k))**2/2/fpig
               engy=engy-(phi(i,jp,k)-phi(i,j,k))**2/2/fpig
               engy=engy-(phi(i,j,kp)-phi(i,j,k))**2/2/fpig
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
      return
      end

      
c
c


c*********************************************************************
      subroutine initialize
      include 'relaxgl.fi'
      integer omp_get_max_threads
c$    external omp_get_max_threads
#ifdef COLD
      include 'cold.fi'
      integer i,j,k

c$omp parallel private(i,j,k) shared(cold)
      do k=1,ng3
c$omp do
         do j=1,ng2
            do i=1,ng1
               cold(i,j,k)=.true.
            enddo
         enddo
c$omp enddo nowait
      enddo
c$omp end parallel
#endif      

c      call wmemory
#ifdef EXACTENERGY
      write(*,*) 'enforcing spatially exact energy conservation'
#else
      write(*,*)' WARNING: no exact spatial energy conservation!'
#endif
c$    write(*,*) ' OMP_max_threads=',omp_get_max_threads()
#ifdef _SGI_SOURCE
       if (.false.) then
       write(*,*) 'running on SGI_SOURCE, with floating point traps'
       write(*,*) 'this may cause unexpected crashes at high opt level'
       call set_fpe()
       else
       write(*,*) 'not using float traps on SGI to allow aggressive opt'
       endif
#endif      

      if (ng3 .lt. 16) then
         write(*,*) 'ERROR: temporary buffers require 16 <= ng3=',ng3
         stop
      endif

      return
      end










