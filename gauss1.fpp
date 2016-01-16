c -*- Fortran -*-
c File Seidel.f
c written May 1994 by Ue-Li Pen, upen@astro.princeton.edu
c modified Nov 17 1995 to use static memory      
c We now use the global header files to allocate the biggest
c necessary amount of memory for the bottom level grid, and
c simply only use a smaller fraction of it at higher levels.
c    
#ifdef _T3D
#define BAJ6 8
#else
#define BAJ6 6
#endif  
#ifdef _SX5
c the extra padding reduces bank conflicts.
#define BAJI2 9+BAJ6
#else
#define BAJI2 8+BAJ6      
#endif
c for test runs uncomment next two lines
#ifdef TESTRELAX    
      call testrelax
      end
#endif
c
c Compiler options note:  on SGI, compile with
c      -WK,-ndr,-directives=-AK,-roundoff=2
c     so that it will ignore the convex and cray directives.
c     
c
c general styles and paradigms:
c
c temporary arrays:
c     are all allocated using the index of a common block stub.
c The interface wrapper then calls the actual work routine using
c the appropriate sectioning of the temporary arrays.
c This will result in an index violation at run time.  If compiling with
c array check, one needs to turn index checking off the the wrappers.
c
c Arguments:  I always try to put them in the order:
c  0.  everything < Temporary arrays      
c  1.  Arrays < Scalars
c  2.  intent(out) < intent(inout) < intent(in)
c  3.  floats < integers < logicals
c  4.  for arrays: bigger rank < smaller rank
c  5.  Arrays:     bigger size < smaller size
c  6.  For integers: dimensions < other
c  7.  Dimensions: same order as the arrays they declare
c
c Limitations:
c  the number of rows (na3) must be ab even multiple of nrelax.
c     
c  Machine specifics:
c    I have inserted vectorization and parallelization directives
c     for convex fc-8.0.  The algorithm is fully data parallel, so
c     it should run at full speed on a parallel/vector machine.
c    The convex compilers have some peculiar bugs, which parallelization
c     of the second layer loops due to some apparent varying inner trip
c     count.  This is fixed by vectorizing the second and parallelizing
c     the inner one.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               
      subroutine relax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu
     &     ,nrelax,iopt)
      implicit none
      integer iopt,na1,na2,na3,nrelax,nu
      real arr(na1,na2,na3), resid((na1+1)/2,(na2+1)/2,(na3+1)/2)
     &     , rhs(na1,na2,na3), rinf, u(nu,na1,na2,na3)
     &     ,deformpot(na1,na2,na3),dx
cdir$ shared *arr(:block,:block,:),*deformpot(:block,:block,:)
cdir$ shared *rhs(:block,:block,:),*resid(:block,:block,:)
cdir$ shared *u(:,:block,:block,:)
c
c
c perform a four color Gauss-Seidel relaxation given only the
c deformation potential at some reduced grid scale.  Optionally compute
c residual.
c
c
c the idea is to use a 2-D work array, so we maintain parallel/vector
c efficiency.
c Parameters:
c       arr(na1,na2,na3)                        the unknown field
c       deformpot(na1,na2,na3)   deformation potential
c       nrelax                           number of Gauss-Seidel iterations.
c                                        zero is legal.
c       resid(na1,na2,na3)              optional residual array
c       dx                              delta x, the grid spacing 
c       integer iopt                    the operation to be performed:
c   mod(iopt,8): 0 -- dont calculate residual
c                1 -- calculate residual on half grid size
c                2 -- return infinity norm of residual in rinf
c                3 -- return infinity norm in rinf, and residual
c                            in resid(,,)
c     The residual has the same sign as RHS, i.e. rhs-L[u].
c      
c   Let j=mod(iopt/8,8).  Then the procedure does one of the following for j:
c       1: regular Poisson iteration using metric g_{ab}
c       2: potential flow solver using triad e^a_b
c       3: implicit hydro solver (not implemented).
c
c  Return codes (passed in iopt):
c       iopt>=0:        successfully completed relaxation
c       iopt = -1:      invalid parameter
c
c it is assumed that all indeces integer powers of 2.
c
c ANSI Fortran-77 issues: I use ENDDO, long variable names, INCLUDE
c and IMPLICIT NONE.
c I have tried to make the variable names unique, even if truncated
c to six characters.  I do not know of compilers where it fails.
c The IMPLICIT NONE statement can simply be removed if not supported.
c The ENDDOs would need to be replaced by continue statements.
c
c parallel/porting issues: temporary arrays are allocated on a common
c block.  To be more portable, one could simply eliminate the common lines.
c

c up equal or less space than reals.
c      
c levels are counted increasing down the pyramid:
c  1                  X
c  2                  X   X
c  3                  X X X X
c etc.      
c

c locals:
      integer i,nrmax,iiopt, nrmleft,itmax


c there are problems when nrelax is either too small, there will be
c problems with index wraparound.
c      
      nrmax = min(4,max(1,na3/4))
c somehow, it doesnt like nrelax=2 for na3=8, so simply surpress it.      
      nrmleft=mod(nrelax,nrmax)


c      write(*,*)'RELAX:def=',deformpot
c      write(*,*)'RELAX:rhs=',rhs
c      write(*,*)'RELAX:arr=',arr


      if (nrelax .gt. nrmax) then
c dont calculate residuals until the last sweep         
         iiopt=iopt-mod(iopt,8)
         itmax=nrelax/nrmax
         do i=1,itmax
            if (i .eq. itmax .and. nrmleft .eq. 0) iiopt=iopt
            call r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu
     &            ,nrmax,iiopt)
         enddo
         if ( nrmleft .gt. 0) then
            call r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu
     &            ,nrmleft,iopt)
         endif
      else
         call r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu
     &         ,nrelax,iopt)
      endif
c      write(*,*)'RELAX(exit):arr=',arr
      return
      end


c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
#ifdef _ALPHA
      subroutine r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3
     &      ,nu,nrelax,iopt)
#include "dimen.fh"
      real flat(BAJ6*(BAJI2)*NG*NG+8192)
      real flat2(BAJ6*(BAJI2)*NG*NG+8192)
c locals
      integer iwordsize,ioff,ichk,ioff2
      logical firsttime,firsthalf
      save firsttime,ioff,ioff2,firsthalf
      data firsttime /.true./
      data firsthalf /.true./

cdec$ alias get_page_size, "getpagesize"
      integer get_page_size,ipagesize
      external get_page_size

      if (firsttime) then
         firsttime=.false.
      ipagesize=get_page_size()
      iwordsize=%loc(flat(2))-%loc(flat(1))
      ioff=(ipagesize-mod(%loc(flat(1)),ipagesize))/iwordsize+1
      ioff2=(ipagesize-mod(%loc(flat2(1)),ipagesize))/iwordsize+1
      write(*,*) 'r1lax: ioff=',ioff,ioff2
cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(flat)
      endif
      if (firsthalf .and. na1 .eq. NG/2) then
         firsthalf=.false.
         write(*,*) 'aligning halfsize baj',na1,na2,na3
cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(flat2)
      endif

      
      if (na1 .eq. NG) then
      call r2lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3
     &      ,nu,nrelax,iopt,flat(ioff))
      else
      call r2lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3
     &      ,nu,nrelax,iopt,flat2(ioff2))
      endif
      return
      end

      subroutine r2lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3
     &      ,nu,nrelax,iopt,baj)
#else
      subroutine r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3
     &      ,nu,nrelax,iopt)
#endif
      implicit none
      integer iopt,na1,na2,na3,nrelax,iclvls,it,nu
      real arr(na1,na2,na3), resid((na1+1)/2,(na2+1)/2,(na3+1)/2)
     & , rhs(na1,na2,na3),rinf,u(nu,na1,na2,na3)
     &     ,deformpot(na1,na2,na3),dx
cdir$ shared *arr(:block,:block,:),*deformpot(:block,:block,:)
cdir$ shared *rhs(:block,:block,:),*resid(:block,:block,:)
cdir$ shared *u(:,:block,:block,:)
c
c Locals
c
#include "dimen.fh"
      integer idxf(NG)
#define MAXRELAX 4
c I think I actually only use baj(1:2*nrelax+1)
#      
#ifdef DYNAMIC_MEMORY
      real baj(BAJ6,BAJI2,na1,na2)
cdir$ shared baj(:,:,:block,:block)
#else      
c      real baj(BAJ6,2*(MAXRELAX+1/(MAXRELAX+1))+BAJ6,NG,NG)
      real baj(BAJ6,BAJI2,NG,NG)
#endif      
      logical lsquare
      integer k,nfudge,ki,ko,idateline,kp,iresidopt

#ifdef _ALPHA
ccdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(baj)
#endif

c always a special case for nrelax=0
      nfudge=2*(nrelax+1/(nrelax+1))+6
      nfudge=min(nfudge,na3)
      if (iopt .lt. 0 .or. iopt .gt. 256 ) then
         write(*,*) 'relax: iopt exceeds legal range:',iopt
         pause
      endif


      if (mod(na3,2) .ne. 0) then
         write(*,*)'relax: dimensions must be even'
c         tmp=0
c         tmp=1/tmp
         pause
      endif

      if (mod(iopt/8,8) .eq. 1) then
         lsquare=.true.
      else
         lsquare=.false.
      endif
      rinf=0
      iresidopt=mod(iopt,8)
      if (iresidopt .eq. 4) rinf=1.0e20
      
c computation strategy:
c the basic idea is red-black Gauss-Seidel.
c  There is a significant operation count to construct the triad or
c  metric, so we need to retain it for each of the sweeps.
c  Temporaries are kept in planes
c
c this looks rather complicated.  To illustrate an example
c with nrelax=2:
c
c    18 20  6 10  9 14 13 17 16 19
c     1  3  2  5  4  8  7 12 11 15
c----------------------------------
c     1  2  3  4  5  6  7  8  9 10     
c
c or with nrelax=3  (which we will use in the examples)
c
c      
c    44 47 46 48 15 21 20 27 26 33 32 38 37 42 41 45
c    39 43  6 10  9 14 13 19 18 25 24 31 30 36 35 40
c     1  3  2  5  4  8  7 12 11 17 16 23 22 29 28 34
c---------------------------------------------------
c     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
c
c     which unraveled is the sequence:
c 1 3 2 5 4 3 7 6 5 4 9 8 7 6 5 11 10 9 8 7 6 13 12 11 10 9 8 15 14 13 12 11 10
c 16 15 14 13 12 1 16 15 14 2 1 16 3 2 4
c      
      if (iresidopt .eq. 1 .or. iresidopt .eq. 3) then
	 call szero(resid,(na1+1)/2,(na2+1)/2,(na3+1)/2)
      endif


c Now relax first planes

      idxf(na3-1)=nfudge-1
      idxf(na3)=nfudge
      
c first set up metric matrix for 3 tiers:
c  do the bottom plane.      

      idateline=1
      idxf(1)=idateline
      idateline=mod(idateline,nfudge)+1
      idxf(2)=idateline
c note that calcbajxy must be called with an odd plane
c plow na3,1      
      call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge
     &    ,1,lsquare)

      if (nrelax .eq. 0) then
         if (iresidopt .eq. 0) then
            write(*,*) 'relax: nothing to do!'
         endif
         do k=2,na3,2
            idateline=mod(idateline,nfudge)+1
            kp=mod(k,na3)+1
            idxf(kp)=idateline
            idateline=mod(idateline,nfudge)+1
            kp=mod(kp,na3)+1
            idxf(kp)=idateline
            kp=mod(k,na3)+1
            call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu
     &           ,nfudge,kp,lsquare)
            call calcresid(resid,baj,arr,rhs,idxf,rinf,dx
     &           ,na1,na2,na3,nfudge,k,iresidopt)
         enddo
         return
      endif


      idateline=mod(idateline,nfudge)+1
      idxf(3)=idateline
      idateline=mod(idateline,nfudge)+1
      idxf(4)=idateline
c and 2:3      
      call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge
     &    ,3,lsquare)



c the initial pyramid:
c      
c                15
c           6 10  9 14 13
c     1  3  2  5  4  8  7 12 11
c---------------------------------------------------
c     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
c
c so this is #1: (1)
      k=1
      call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)

      do ko=2,(nrelax-1)*4+1,2
         iclvls=(ko-2)/4
         idateline=mod(idateline,nfudge)+1
         kp=mod(ko+2,na3)+1
         idxf(kp)=idateline
         idateline=mod(idateline,nfudge)+1
         kp=mod(ko+3,na3)+1
         idxf(kp)=idateline
         kp=mod(ko+2,na3)+1
         call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge
     &        ,kp,lsquare)
c to visualize: when ko=2, iclvls=0, and we do #2, #3 ( 3 2 )
c second time around, ko=4, iclvls=0, do       #4, #5 ( 5 4 )
c and the third time, ko=6, iclvls=1, do       #7 #8 #9 #10 (7 6 5 4)
c and so on.         
         do ki=0,iclvls
            k=mod(ko-2*ki+4*na3,na3)+1
            call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
            k=mod(ko-2*ki+4*na3-1,na3)+1
            call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
         enddo
c on the second time, fill #6 (3), etc 
         if (mod(ko,4) .eq. 0) then
            k=mod(ko/2,na3)+1
            call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
         endif
      enddo

      
c do intermediate plane
c     
c                   21 20 27 26 33 32
c                         19 18 25 24 31 30
c                               17 16 23 22 29 28 
c---------------------------------------------------
c     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
c      
      do ko=(nrelax-1)*4+2,na3-1,2
         k=ko
         kp=mod(ko+2,na3)+1
         idateline=mod(idateline,nfudge)+1
         idxf(kp)=idateline
         kp=mod(kp,na3)+1
         idateline=mod(idateline,nfudge)+1
         idxf(kp)=idateline
         k=mod(ko+2,na3)+1
         call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge
     &    ,k,lsquare)
         do ki=1,nrelax*2,2
            k=ko-ki+2
            call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
            k=ko-ki+1
            call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
         enddo
         if (iresidopt .gt. 0 .and. ko .gt. (nrelax-1)*4+2 ) then
c we want tier 6 for nrelax=3
c and tier 2 for nrelax=1            
            k=mod(4*na3+ko-nrelax*2+1,na3)+1
            call calcresid(resid,baj,arr,rhs,idxf,rinf,dx
     &           ,na1,na2,na3,nfudge,k,iresidopt)
         endif
      enddo
      
c now clean up the leftovers: the Finale
c     
c    44 47 46 48                      38 37 42 41 45
c    39 43                                  36 35 40
c                                                 34
c---------------------------------------------------
c     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
c
c
      if (.false.) then
      idateline=mod(idateline,nfudge)+1
      idxf(1)=idateline
      idateline=mod(idateline,nfudge)+1
      idxf(2)=idateline
c note that calcbajxy must be called with an odd plane
c plow na3,1      
      call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge
     &    ,1,lsquare)
      endif
      k=na3
c do #34 (16)
      call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
c loop through #35 #36 #37 #38 (15 14 13 12)
      do ki=2,nrelax
         k=mod(4*na3-2*ki+2,na3)+1
         call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
         k=mod(4*na3-2*ki+1,na3)+1
         call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
      enddo
      if (iresidopt .gt. 0) then
         k=mod(4*na3-nrelax*2+1,na3)+1
         call calcresid(resid,baj,arr,rhs,idxf,rinf,dx
     &        ,na1,na2,na3,nfudge,k,iresidopt)
      endif
      k=3
      idateline=mod(idateline,nfudge)+1
      idxf(k)=idateline
      k=4
      idateline=mod(idateline,nfudge)+1
      idxf(k)=idateline
      k=3
      call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge
     &    ,k,lsquare)
      do ko=1,nrelax-1
c the first time we enter the ki loop, it=1, ko=1, and we do a full
c step sweep: #39 #40 #41 #42 (1 16 15 14)
c the second time round, it=2, ko=1, and the first step only has one
c brick left.
         k=mod(4*na3+2*ko+2,na3)+1
         idateline=mod(idateline,nfudge)+1
         idxf(k)=idateline
         k=mod(4*na3+2*ko+3,na3)+1
         idateline=mod(idateline,nfudge)+1
         idxf(k)=idateline
         k=mod(4*na3+2*ko+2,na3)+1
         call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu
     &              ,nfudge,k,lsquare)
         do it=1,2
            do ki=1,nrelax-ko
c we update three new rows each iteration, so we need to alternate
c one and two double calcbajxys.               
               if ( it .ne. 2 .or. ki .ne. 1 ) then
                  k=mod(4*na3+2*ko-2*ki+2*it-2,na3)+1
                  call relaxplane(arr,baj,rhs,idxf,dx,na1,na2
     &                    ,na3,nfudge,k)
               endif
               k=mod(4*na3+2*ko-2*ki+2*it-3,na3)+1
               call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
            enddo
         enddo
c end it=1,2         
         if (iresidopt .gt. 0) then
c update the remaining residuals at the top and bottom:
c normalize by: if nrelax=2, ko=1, we want k=na3
            k=mod(4*na3+4*ko-2*nrelax-1,na3)+1
            call calcresid(resid,baj,arr,rhs,idxf,rinf,dx
     &           ,na1,na2,na3,nfudge,k,iresidopt)
            k=mod(4*na3+4*ko-2*nrelax+1,na3)+1
            call calcresid(resid,baj,arr,rhs,idxf,rinf,dx
     &           ,na1,na2,na3,nfudge,k,iresidopt)
         endif
      enddo
c end ko=1,nrelax
      if (iresidopt .gt. 0) then
          k=mod(k+1,na3)+1
          if (nrelax .eq. 1) then
            k=2
          endif
          call calcresid(resid,baj,arr,rhs,idxf,rinf,dx
     &           ,na1,na2,na3,nfudge,k,iresidopt)
       endif
c
c do you really believe that worked?
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c compute the residual on a coarsened grid
c
c when called on an even grid, it calculates the residual at ki/2,
c and side effects resid(ki/2+1).  It uses accesses arrays at levels
c ki-2, ki-1, ki, ki+1      
c
      subroutine calcresid(resid,fudge,phi,rhs,idxf,rinf,dx
     &    ,na1,na2,na3,nfudge,ki,iopt)

c WORKER routine.  This is a very expensive routine.
c
c operation count: 46.375 A  and 26.875 M for nred=2
c typically, iter=4, this is called every eight relax iterations,
c and the effective cost should be 9.0312 FLOP/element to be
c added to the 67 for planerelax.      
c In practice, most machines spend much more time than this
c requires.  An exception is the cray, which does this very efficiently.
c      
      implicit none
      integer na1,na2,na3,nfudge,idxf(na3),ki,iopt
      real resid((na1+1)/2,(na2+1)/2,(na3+1)/2),phi(na1,na2,na3)
     &    ,fudge(BAJ6,BAJI2,na1,na2),rhs(na1,na2,na3),dx,rinf
cdir$ shared *rhs(:block,:block,:),*resid(:block,:block,:)
cdir$ shared *fudge(:,:,:block,:block),*phi(:block,:block,:)
c locals
      integer i,j,km,ko,kp,k,jp,jm,ip,im,ib,ibp,ibm,jb,jbp,jbm,kb,kbp
#ifdef DYNAMIC_MEMORY
      real rt(na1,na2,2)
cdir$ shared rt(:block,:block,:)
#else      
      real rt(NG,NG,2)
#endif
      real b11p,b11m,b22p,b22m,b33p,b33m,tdiag,t1,t2
     &     ,di2,tresid,rinfmax,rinfmin
      logical lnorminf,lresid
#ifdef T3D
      real rinf1
cdir$ shared rinf1
#endif

      if (mod(ki,2) .eq. 1) then
          write(*,*) 'calcresid called with odd k'
c          di2=0
c          di2=1/di2
          pause
      endif


      lresid=.false.
      lnorminf=.false.
      if (mod(iopt,2) .eq. 1) lresid=.true.
      if (mod(iopt/2,2) .eq. 1) then
         lnorminf=.true.
      endif
      di2=1/dx**2
#if   defined(_SGI_SOURCE) || defined(_SX5)
c the SGI does not know how to promote tests!
      rinfmax=0
      rinfmin=1.e10
c 7.2 f77
c*$* assert do(serial)
c$omp parallel default(private) firstprivate(ki,na1,na2,na3,di2,lresid)
c$omp& shared(fudge,idxf,phi,rhs,rt) reduction(max:rinfmax)
c$omp& reduction(min:rinfmin)  
      do ko=1,2
         k=ki-2+ko
         kp=mod(k,na3)+1
         km=mod(k+na3-2,na3)+1
c help the convex:            
c$dir prefer_vector
c help the SGI:         
cc*$* assert do(concurrent)         
c$doacross local(i,j,TRESID,t1,t2,tdiag,b33m,b33p,b22m,b22p,b11m,b11p
c$&         ,im,ip,j,jp,jm), reduction(rinfmax,rinfmin)
c$omp do
         do j=1,na2
            jp=mod(j,na2)+1
            jm=mod(j+na2-2,na2)+1
c$dir prefer_parallel_ext            
c*$* assert do(serial)
            do i=1,na1
               ip=mod(i,na1)+1
               im=mod(i+na1-2,na1)+1
c calculate a symmetric Laplacian
c \partial_a G^{ab} \partial_b \phi
c where G^{ab} is a 3x3 symmetric matrix given by fudge.
c We build a symmetric discretization by splitting the tensor
c sum into the diagonal and offdiagonal pieces, t1 and t2.
c The diagonal stencil is as follows:
c               
c                      o
c                     b22p 
c               o b11m o b11p o
c                     b22m
c                      o
c
c we get (\phi_{i+1}-\phi_i)*b_{11}^p-(\phi_i-\phi_{i-1})*b_{11}^m
c and so on.
c     
               b11p=(fudge(1,idxf(k),ip,j)+fudge(1,idxf(k),i,j))/2
               b11m=(fudge(1,idxf(k),im,j)+fudge(1,idxf(k),i,j))/2
               b22p=(fudge(4,idxf(k),i,jp)+fudge(4,idxf(k),i,j))/2
               b22m=(fudge(4,idxf(k),i,jm)+fudge(4,idxf(k),i,j))/2
               b33p=(fudge(6,idxf(kp),i,j)+fudge(6,idxf(k),i,j))/2
               b33m=(fudge(6,idxf(km),i,j)+fudge(6,idxf(k),i,j))/2
c 6 A,  6 M
               tdiag=-(b11p+b11m+b22p+b22m+b33p+b33m)
               t1=(b11p*phi(ip,j,k)+b11m*phi(im,j,k)+b22p*phi(i,jp,k)
     &              +b22m*phi(i,jm,k)+b33p*phi(i,j,kp)+b33m*phi(i,j,km))
c 10 A, 6 M              
c t2 contains the diagonal terms
               t2=(  (phi(ip,jp,k)-phi(im,jp,k))*fudge(2,idxf(k),i,jp)
     &              -(phi(ip,jm,k)-phi(im,jm,k))*fudge(2,idxf(k),i,jm)
     &              +(phi(ip,jp,k)-phi(ip,jm,k))*fudge(2,idxf(k),ip,j)
     &              -(phi(im,jp,k)-phi(im,jm,k))*fudge(2,idxf(k),im,j)
     &              
     &              +(phi(i,jp,kp)-phi(i,jm,kp))*fudge(5,idxf(kp),i,j)
     &              -(phi(i,jp,km)-phi(i,jm,km))*fudge(5,idxf(km),i,j)
     &              +(phi(i,jp,kp)-phi(i,jp,km))*fudge(5,idxf(k),i,jp)
     &              -(phi(i,jm,kp)-phi(i,jm,km))*fudge(5,idxf(k),i,jm)
     &              
     &              +(phi(ip,j,kp)-phi(im,j,kp))*fudge(3,idxf(kp),i,j)
     &              -(phi(ip,j,km)-phi(im,j,km))*fudge(3,idxf(km),i,j)
     &              +(phi(ip,j,kp)-phi(ip,j,km))*fudge(3,idxf(k),ip,j)
     &              -(phi(im,j,kp)-phi(im,j,km))*fudge(3,idxf(k),im,j)
     &              )/4
c 3*8-1=23 A, 13 M
               tresid=rhs(i,j,k)-(t1+t2+tdiag*phi(i,j,k))*di2
c 3 A, 2 M
               rinfmax=max(rinfmax,abs(tresid))
               if (lresid) rt(i,j,ko)=tresid
               rinfmin=min(rinfmin,abs(tresid))
            enddo
         enddo
      enddo
c$omp end parallel
      if (lnorminf) rinf=max(rinf,rinfmax)
      if (iopt .eq. 4) rinf=min(rinf,rinfmin)
#else
c zpj: I parallelized the following section by divide the condition
c zpj:  sentences if (lnormin) and if (iopt .eq.4 ) into two separate
c zpj:  section. (1)  if (lnormin) and (2) if (iopt .eq.4 ) because 
c zpj: (1) -> not (2) and (2) -> not (1)

      if (lnorminf) then
c 7.2 f77
c$omp parallel default(private) shared(rt,fudge,phi,rhs)
c$omp& firstprivate(idxf,na1,na2,na3,di2,ki,lresid) reduction(max:rinf)
c*$* assert do(serial)
      do ko=1,2
         k=ki-2+ko
         kp=mod(k,na3)+1
         km=mod(k+na3-2,na3)+1
c help the convex:         
c$dir prefer_vector
cdir$ doshared(j,i) on rhs(i,j,1)      
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,na2
            do i=1,na1
c$dir prefer_parallel_ext            
               jp=mod(j,na2)+1
               jm=mod(j+na2-2,na2)+1
               ip=mod(i,na1)+1
               im=mod(i+na1-2,na1)+1
c calculate a symmetric Laplacian
c \partial_a G^{ab} \partial_b \phi
c where G^{ab} is a 3x3 symmetric matrix given by fudge.
c We build a symmetric discretization by splitting the tensor
c sum into the diagonal and offdiagonal pieces, t1 and t2.
c The diagonal stencil is as follows:
c               
c                      o
c                     b22p 
c               o b11m o b11p o
c                     b22m
c                      o
c
c we get (\phi_{i+1}-\phi_i)*b_{11}^p-(\phi_i-\phi_{i-1})*b_{11}^m
c and so on.
c     
               b11p=(fudge(1,idxf(k),ip,j)+fudge(1,idxf(k),i,j))/2
               b11m=(fudge(1,idxf(k),im,j)+fudge(1,idxf(k),i,j))/2
               b22p=(fudge(4,idxf(k),i,jp)+fudge(4,idxf(k),i,j))/2
               b22m=(fudge(4,idxf(k),i,jm)+fudge(4,idxf(k),i,j))/2
               b33p=(fudge(6,idxf(kp),i,j)+fudge(6,idxf(k),i,j))/2
               b33m=(fudge(6,idxf(km),i,j)+fudge(6,idxf(k),i,j))/2
c 6 A,  6 M
               tdiag=-(b11p+b11m+b22p+b22m+b33p+b33m)
               t1=(b11p*phi(ip,j,k)+b11m*phi(im,j,k)+b22p*phi(i,jp,k)
     &              +b22m*phi(i,jm,k)+b33p*phi(i,j,kp)+b33m*phi(i,j,km))
c 10 A, 6 M              
c t2 contains the diagonal terms
               t2=(  (phi(ip,jp,k)-phi(im,jp,k))*fudge(2,idxf(k),i,jp)
     &              -(phi(ip,jm,k)-phi(im,jm,k))*fudge(2,idxf(k),i,jm)
     &              +(phi(ip,jp,k)-phi(ip,jm,k))*fudge(2,idxf(k),ip,j)
     &              -(phi(im,jp,k)-phi(im,jm,k))*fudge(2,idxf(k),im,j)
     &              
     &              +(phi(i,jp,kp)-phi(i,jm,kp))*fudge(5,idxf(kp),i,j)
     &              -(phi(i,jp,km)-phi(i,jm,km))*fudge(5,idxf(km),i,j)
     &              +(phi(i,jp,kp)-phi(i,jp,km))*fudge(5,idxf(k),i,jp)
     &              -(phi(i,jm,kp)-phi(i,jm,km))*fudge(5,idxf(k),i,jm)
     &              
     &              +(phi(ip,j,kp)-phi(im,j,kp))*fudge(3,idxf(kp),i,j)
     &              -(phi(ip,j,km)-phi(im,j,km))*fudge(3,idxf(km),i,j)
     &              +(phi(ip,j,kp)-phi(ip,j,km))*fudge(3,idxf(k),ip,j)
     &              -(phi(im,j,kp)-phi(im,j,km))*fudge(3,idxf(k),im,j)
     &              )/4
c 3*8-1=23 A, 12 M
               tresid=rhs(i,j,k)-(t1+t2+tdiag*phi(i,j,k))*di2
c 3 A, 2 M
c zpj: I send these two if ( if (lnorminf) before max and if (iopt
c               .eq. 4) before min ahead of the parallel unit) 
c zpj: the original command is:
c zpj:         rinf=max(rinf,abs(tresid)) 
c zpj:         if (lresid) rt(i,j,ko)=tresid
c zpj:         rinf=max(rinf,abs(tresid)

               rinf=max(rinf,abs(tresid)) 
               if (lresid) rt(i,j,ko)=tresid
            enddo
         enddo
      enddo
c$omp end parallel
      endif
      if (iopt .eq. 4) then
c 7.2 f77
c$omp parallel default(private) shared(rt,fudge,phi,rhs)
c$omp& firstprivate(idxf,na1,na2,na3,di2,ki,lresid) reduction(min:rinf)
c*$* assert do(serial)
      do ko=1,2
         k=ki-2+ko
         kp=mod(k,na3)+1
         km=mod(k+na3-2,na3)+1
c help the convex:         
c$dir prefer_vector
cdir$ doshared(j,i) on rhs(i,j,1)      
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,na2
            do i=1,na1
c$dir prefer_parallel_ext            
               jp=mod(j,na2)+1
               jm=mod(j+na2-2,na2)+1
               ip=mod(i,na1)+1
               im=mod(i+na1-2,na1)+1
c calculate a symmetric Laplacian
c \partial_a G^{ab} \partial_b \phi
c where G^{ab} is a 3x3 symmetric matrix given by fudge.
c We build a symmetric discretization by splitting the tensor
c sum into the diagonal and offdiagonal pieces, t1 and t2.
c The diagonal stencil is as follows:
c               
c                      o
c                     b22p 
c               o b11m o b11p o
c                     b22m
c                      o
c
c we get (\phi_{i+1}-\phi_i)*b_{11}^p-(\phi_i-\phi_{i-1})*b_{11}^m
c and so on.
c     
               b11p=(fudge(1,idxf(k),ip,j)+fudge(1,idxf(k),i,j))/2
               b11m=(fudge(1,idxf(k),im,j)+fudge(1,idxf(k),i,j))/2
               b22p=(fudge(4,idxf(k),i,jp)+fudge(4,idxf(k),i,j))/2
               b22m=(fudge(4,idxf(k),i,jm)+fudge(4,idxf(k),i,j))/2
               b33p=(fudge(6,idxf(kp),i,j)+fudge(6,idxf(k),i,j))/2
               b33m=(fudge(6,idxf(km),i,j)+fudge(6,idxf(k),i,j))/2
c 6 A,  6 M
               tdiag=-(b11p+b11m+b22p+b22m+b33p+b33m)
               t1=(b11p*phi(ip,j,k)+b11m*phi(im,j,k)+b22p*phi(i,jp,k)
     &              +b22m*phi(i,jm,k)+b33p*phi(i,j,kp)+b33m*phi(i,j,km))
c 10 A, 6 M              
c t2 contains the diagonal terms
               t2=(  (phi(ip,jp,k)-phi(im,jp,k))*fudge(2,idxf(k),i,jp)
     &              -(phi(ip,jm,k)-phi(im,jm,k))*fudge(2,idxf(k),i,jm)
     &              +(phi(ip,jp,k)-phi(ip,jm,k))*fudge(2,idxf(k),ip,j)
     &              -(phi(im,jp,k)-phi(im,jm,k))*fudge(2,idxf(k),im,j)
     &              
     &              +(phi(i,jp,kp)-phi(i,jm,kp))*fudge(5,idxf(kp),i,j)
     &              -(phi(i,jp,km)-phi(i,jm,km))*fudge(5,idxf(km),i,j)
     &              +(phi(i,jp,kp)-phi(i,jp,km))*fudge(5,idxf(k),i,jp)
     &              -(phi(i,jm,kp)-phi(i,jm,km))*fudge(5,idxf(k),i,jm)
     &              
     &              +(phi(ip,j,kp)-phi(im,j,kp))*fudge(3,idxf(kp),i,j)
     &              -(phi(ip,j,km)-phi(im,j,km))*fudge(3,idxf(km),i,j)
     &              +(phi(ip,j,kp)-phi(ip,j,km))*fudge(3,idxf(k),ip,j)
     &              -(phi(im,j,kp)-phi(im,j,km))*fudge(3,idxf(k),im,j)
     &              )/4
c 3*8-1=23 A, 12 M
               tresid=rhs(i,j,k)-(t1+t2+tdiag*phi(i,j,k))*di2
c 3 A, 2 M
               if (lresid) rt(i,j,ko)=tresid
               rinf=min(rinf,abs(tresid))
            enddo
         enddo
      enddo
c$omp end parallel
      endif
#  ifdef T3D
c this is the standard multiple processor reduction:
cdir$ master
      rinf1=rinf
cdir$ end master
      if (lnorminf) then
cdir$ atomic update
         rinf1=max(rinf1,rinf)
cdir$ barrier
cdir$ suppress
         rinf=rinf1
      else if (iopt .eq. 4) then
cdir$ atomic update
         rinf1=min(rinf1,rinf)
cdir$ barrier
cdir$ suppress
         rinf=rinf1
      endif
#  endif
#endif
c the next section uses (27A+7M)/2/nred**2 operations per element
      if (lresid) then
c$omp parallel default(private) firstprivate(na1,na2,na3,ki)
c$omp& shared(resid,rt)
         k=ki/2
         kp=mod(k,na3/2)+1
         kb=1
         kbp=2
c$dir prefer_vector      
cdir$ doshared(j,i) on rt(i,j,1)      
c 7.2 f77 
c*$* assert do(concurrent)
c$omp do
         do j=1,(na2+1)/2
c*$* assert do(serial)
            do i=1,(na1+1)/2
               ib=2*i-1
               ibp=mod(ib,na1)+1
               ibm=mod(ib+na1-2,na1)+1
               jb=2*j-1
               jbp=mod(jb,na2)+1
               jbm=mod(jb+na2-2,na2)+1
               resid(i,j,k)=resid(i,j,k)+rt(ib,jb,kb)/8
     &              +(rt(ibp,jb,kb)+rt(ibm,jb,kb)+rt(ib,jbp,kb)
     &             +rt(ib,jbm,kb) +rt(ib,jb,kbp))/16
     &              +(rt(ibp,jbp,kb)+rt(ibp,jbm,kb)+rt(ibm,jbp,kb)
     &              +rt(ibm,jbm,kb)+rt(ibp,jb,kbp)+rt(ibm,jb,kbp)
     &              +rt(ib,jbp,kbp)+rt(ib,jbm,kbp))/32
     &              +(rt(ibp,jbp,kbp)+rt(ibp,jbm,kbp)+rt(ibm,jbp,kbp)
     &              +rt(ibm,jbm,kbp))/64
c 18 A 4 M,  but we divide by eight for nred=2
c the kbm part was taken care of in the previous iteration.
               resid(i,j,kp)=resid(i,j,kp)+rt(ib,jb,kbp)/16
     &              +(rt(ibp,jb,kbp)+rt(ibm,jb,kbp)+rt(ib,jbp,kbp)
     &              +rt(ib,jbm,kbp))/32
     &              +(rt(ibp,jbp,kbp)+rt(ibp,jbm,kbp)+rt(ibm,jbp,kbp)
     &              +rt(ibm,jbm,kbp))/64
c 9 A 3 M, ditto         
            enddo
         enddo
c$omp end parallel
      endif
      return
      end
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine relaxplane(phi,fudge,rhs,idxf,dx,na1,na2,na3,nfudge,k)
c
c WORKER routine.  This is the single most expensive routine, typically
c accounting for about 80% of the total execution time.
c
      implicit none
      integer na1,na2,na3,k,nfudge,idxf(na3)
      real phi(na1,na2,na3),rhs(na1,na2,na3)
      real fudge(BAJ6,BAJI2,na1,na2),dx
cdir$ shared *rhs(:block,:block,:),*phi(:block,:block,:)
cdir$ shared *fudge(:,:,:block,:block)
c
c Operation count: Floating point - 41 Additions, 26 Multiplicatins, 1 Division
c   per element (=> 67 FLOP),    plus various integer stuff
c      
c This should be compared to 6 Additions and 1 Multiplication
c      for a Euclidian Poisson iteration.
c      
c locals
      real dx2,t1,t2,tdiag,b11p,b11m,b22p,b22m,b33p,b33m
      integer i,j,io,jo,ip,im,jp,jm,kp,km
#ifdef _SX5
       integer ij
#endif
#ifdef _ALPHA
cdec$ alias addr_to_rad, "_OtsAddrToRad"
      integer  addr_to_rad,irad,irad1,irad2
      external addr_to_rad
      integer omp_get_thread_num,cpuid,cpu_get_rad
      external omp_get_thread_num,cpuid,cpu_get_rad
#endif

c executable section
cc$omp parallel default(private) 
cc$omp& shared(fudge,phi,rhs,dx,na1,na2,na3,idxf,k)
      dx2=dx**2
      kp=mod(k,na3)+1
      km=mod(k+na3-2,na3)+1

c tell the cray preprocessor that there is no data dependency in this routine.
cfpp$ nodepchk r
c Integer division means that io,jo loop to 1 unless na=1.
cfpp$ skip
c$dir scalar      
      do io=0,1-1/na1
cfpp$ skip
c$dir scalar         
         do jo=0,1-1/na2
c on the SGI Mipspro 7.2 compiler, the no recurrence directive
c no longer exists.  There the doacross will still force parallelization.
cfpp$ select concur
c$dir no_recurrence, force_parallel_ext
cdir$ doshared(j,i) on phi(i,j,1)      
cc*$* assert no recurrence(phi)
c 7.2 f77
c*$* assert do(concurrent)
#ifdef _SX5
!cdir nodep
             do ij=0,na2*na1/4-1
                i=mod(ij,na1/2)*2+1+io
                j=ij*2/na1
               j=j*2+1+jo
#else  
cc$omp do       
C$OMP PARALLEL DO
C$OMP+PRIVATE (i,j,ip,im,jp,jm,b11p,b11m,b22p,b22m,
C$OMP+         b33p,b33m,tdiag,t1,t2)
C$OMP+SHARED (k,idxf,phi,fudge,rhs,dx2,na1,na2,io,jo)
            do j=1+jo,na2,2
#ifdef _ALPHAXXXXX
            if (k .eq. 1 .and. na2 .ge. 128) then
               if(mod(j,na2/32) .eq. 1) then
               irad1=addr_to_rad(rhs(1,j,1))
               irad2=addr_to_rad(fudge(1,1,1,j))
               irad=cpu_get_rad()
               if ((irad1 .ne. irad) .or. (irad2 .ne. irad)) then
                       write(*,*) j,irad1,irad2,irad,na2
               endif
               endif
               endif
#endif
cfpp$ select vector
c$dir no_recurrence, force_vector
cdir$ ivdep
               do i=1+io,na1,2
#endif
               jp=mod(j,na2)+1
               jm=mod(j+na2-2,na2)+1
                  ip=mod(i,na1)+1
                  im=mod(i+na1-2,na1)+1
                  b11p=(fudge(1,idxf(k),ip,j)+fudge(1,idxf(k),i,j))/2
                  b11m=(fudge(1,idxf(k),im,j)+fudge(1,idxf(k),i,j))/2
                  b22p=(fudge(4,idxf(k),i,jp)+fudge(4,idxf(k),i,j))/2
                  b22m=(fudge(4,idxf(k),i,jm)+fudge(4,idxf(k),i,j))/2
                  b33p=(fudge(6,idxf(kp),i,j)+fudge(6,idxf(k),i,j))/2
                  b33m=(fudge(6,idxf(km),i,j)+fudge(6,idxf(k),i,j))/2
c 6 A,  6 M
                  tdiag=-(b11p+b11m+b22p+b22m+b33p+b33m)
              t1=(b11p*phi(ip,j,k)+b11m*phi(im,j,k)+b22p*phi(i,jp,k)
     &            +b22m*phi(i,jm,k)+b33p*phi(i,j,kp)+b33m*phi(i,j,km))
c 10 A, 6 M              
c t2 contains the diagonal terms
              t2=(   (phi(ip,jp,k)-phi(im,jp,k))*fudge(2,idxf(k),i,jp)
     &              -(phi(ip,jm,k)-phi(im,jm,k))*fudge(2,idxf(k),i,jm)
     &              +(phi(ip,jp,k)-phi(ip,jm,k))*fudge(2,idxf(k),ip,j)
     &              -(phi(im,jp,k)-phi(im,jm,k))*fudge(2,idxf(k),im,j)
     &
     &              +(phi(i,jp,kp)-phi(i,jm,kp))*fudge(5,idxf(kp),i,j)
     &              -(phi(i,jp,km)-phi(i,jm,km))*fudge(5,idxf(km),i,j)
     &              +(phi(i,jp,kp)-phi(i,jp,km))*fudge(5,idxf(k),i,jp)
     &              -(phi(i,jm,kp)-phi(i,jm,km))*fudge(5,idxf(k),i,jm)
     &          
     &              +(phi(ip,j,kp)-phi(im,j,kp))*fudge(3,idxf(kp),i,j)
     &              -(phi(ip,j,km)-phi(im,j,km))*fudge(3,idxf(km),i,j)
     &              +(phi(ip,j,kp)-phi(ip,j,km))*fudge(3,idxf(k),ip,j)
     &              -(phi(im,j,kp)-phi(im,j,km))*fudge(3,idxf(k),im,j)
     &              )/4
c 3*8-1=23 A, 13 M
                     phi(i,j,k)=(rhs(i,j,k)*dx2-t1-t2)/tdiag
c 2 A, 1 M, 1 D  
#ifndef _SX5                
               enddo
#endif
            enddo
C$OMP END PARALLEL DO
         enddo
      enddo
cc$omp end parallel
      return
      end

c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c calculate the metric at level k-1 and k, interpolating if necessary.
c



      subroutine calcbajxy(baj,deformpot,u,
     &  idxf,dx,na1,na2,na3,nu,nbaj,k,lsquare)
      implicit none
      integer nu,na1,na2,na3,k,idxf(na3),nbaj
      real u(nu,na1,na2,na3)
      real  baj(BAJ6,BAJI2,na1,na2),deformpot(na1,na2,na3),dx
cdir$ shared *deformpot(:block,:block,:),*u(:,:block,:block,:)
cdir$ shared *baj(:,:,:block,:block)
      logical lsquare
c      
c locals

      integer km

      if (mod(k,2) .ne. 1) then
         write(*,*)' calcbajxy must be called on odd level'
         pause
      endif
      km=mod(k+na3-2,na3)+1
      call cbajraw(baj,deformpot,u,dx
     &        ,na1,na2,na3,nu,nbaj,idxf(k),k,lsquare)
      call cbajraw(baj,deformpot,u,dx
     &        ,na1,na2,na3,nu,nbaj,idxf(km),km,lsquare)
      return
      end
      
      
      subroutine cbajraw(baj,def,u,dx,na1,na2,na3,nu,nbaj,ik,k,lsquare)
c WORKER routine.  The second most expensive call in the whole program.
      implicit none
      logical lsquare
      integer nu,na1,na2,na3,k,ik,nbaj
      real baj(BAJ6,BAJI2,na1,na2),def(na1,na2,na3),dx
     &   ,u(nu,na1,na2,na3)
cdir$ shared *def(:block,:block,:),*u(:,:block,:block,:)
cdir$ shared *baj(:,:,:block,:block)
c calculate the metric in the plane k
c
c Operation count:  28+12 A, 38+18 M, 1 D  =>  97 FLOP/element
c      
c locals
      integer kp,km,j,jp,jm,i,ip,im
      real phixx,phiyy,phizz,phixy,phixz,phiyz,a11,a12,a13,a22,a23,a33
     &     ,b11,b12,b13,b22,b23,b33,di2,det,urg
      logical singular
c      real dtmin
c      save dtmin
c      data dtmin /1./
#ifdef _SX5
       integer ij
#endif
#ifdef GMETRIC
      include 'relaxgl.fi'
      include 'gmetric.fi'

      if (ng1.eq.na1 .and. ng2.eq.na2 .and. ng3.eq.na3) then
         urg=1
c$omp parallel do default(private) firstprivate(na1,na2,na3,k,ik,isquare)
c$omp& shared(gbaj,baj,u)
c$omp do
         do j=1,na2
            do i=1,na1
               det=gbaj(7,i,j,k)
               b11=gbaj(1,i,j,k)
               b12=gbaj(2,i,j,k)
               b13=gbaj(3,i,j,k)
               b22=gbaj(4,i,j,k)
               b23=gbaj(5,i,j,k)
               b33=gbaj(6,i,j,k)
               if (lsquare) then
               baj(1,ik,i,j)=(b11**2+b12**2+b13**2)*det
               baj(2,ik,i,j)=(b11*b12+b12*b22+b13*b23)*det
               baj(3,ik,i,j)=(b11*b13+b12*b23+b13*b33)*det
               baj(4,ik,i,j)=(b12**2+b22**2+b23**2)*det
               baj(5,ik,i,j)=(b12*b13+b22*b23+b23*b33)*det
               baj(6,ik,i,j)=(b13**2+b23**2+b33**2)*det
               else
#ifndef NORHORG
                kp=mod(k,na3)+1
                km=mod(k+na3-2,na3)+1
                jp=mod(j,na2)+1
                jm=mod(j+na2-2,na2)+1
                ip=mod(i,na1)+1
                im=mod(i+na1-2,na1)+1
                urg=
     &      u(1,im,jm,km)+u(1,im,jm,k)+u(1,im,jm,kp)
     &     +u(1,im,j,km)+u(1,im,j,k)+u(1,im,j,kp)
     &     +u(1,im,jp,km)+u(1,im,jp,k)+u(1,im,jp,kp)
     &     +u(1,i,jm,km)+u(1,i,jm,k)+u(1,i,jm,kp)
     &     +u(1,i,j,km)+u(1,i,j,k)+u(1,i,j,kp)
     &     +u(1,i,jp,km)+u(1,i,jp,k)+u(1,i,jp,kp)
     &     +u(1,ip,jm,km)+u(1,ip,jm,k)+u(1,ip,jm,kp)
     &     +u(1,ip,j,km)+u(1,ip,j,k)+u(1,ip,j,kp)
     &     +u(1,ip,jp,km)+u(1,ip,jp,k)+u(1,ip,jp,kp)
               urg=urg/27.
#endif
               baj(1,ik,i,j)=b11*urg
               baj(2,ik,i,j)=b12*urg
               baj(3,ik,i,j)=b13*urg
               baj(4,ik,i,j)=b22*urg
               baj(5,ik,i,j)=b23*urg
               baj(6,ik,i,j)=b33*urg
               endif
            enddo
         enddo
c$omp end parallel
      else
#endif

      di2=1/dx**2
      kp=mod(k,na3)+1
      km=mod(k+na3-2,na3)+1
      singular=.false.
c the 7.2 f77 gets confused if the urg is left dangling in a conditional
c assignment.
      urg=1
c$dir prefer_vector
cdir$ doshared(j,i) on def(i,j,1)      
c 7.2 f77
c*$* assert do(concurrent)
#ifdef _SX5
      do ij=0,na1*na2-1
         i=mod(ij,na1)+1
         j=ij/na1+1
#else
c$omp parallel do default(private)
c$omp&     shared(na1,na2,na3,dx,k,ik,lsquare,def,baj,u,di2,kp,km)
      do j=1,na2
         do i=1,na1
#endif
         jp=mod(j,na2)+1
         jm=mod(j+na2-2,na2)+1
            ip=mod(i,na1)+1
            im=mod(i+na1-2,na1)+1
            phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))*di2
            phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))*di2
            phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))*di2
c 3*( 2 A, 2 M )            
            phixy=(def(ip,jp,k)-def(im,jp,k)
     &           -def(ip,jm,k)+def(im,jm,k))*di2/4
            phiyz=(def(i,jp,kp)-def(i,jp,km)
     &           -def(i,jm,kp)+def(i,jm,km))*di2/4
            phixz=(def(ip,j,kp)-def(im,j,kp)
     &           -def(ip,j,km)+def(im,j,km))*di2/4
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
c 7 A, 11 M, 1 D
#ifndef _SGI_SOURCE
c the SGI compiler has problems            
c            singular = singular .or. det .le. 0
#endif            
c these bij are actually \sqrt{g} e^i_\alpha            
            b11=(a22*a33-a23**2)
            b12=(a13*a23-a12*a33)
            b13=(a12*a23-a13*a22)
            b22=(a11*a33-a13**2)
            b23=(a12*a13-a11*a23)
            b33=(a11*a22-a12**2)
c 6*( 1 A, 2 M )
#ifndef _SGI_SOURCE
            if (det .le. 1.e-5) then
c if a coarse mesh point gets tangled, skip this point.
               b11=1
               b12=0
               b13=0
               b22=1
               b23=0
               b33=1
               det=1
            endif
#endif
            if (lsquare) then
c we want \sqrt{g} g^{\alpha\beta}               
               baj(1,ik,i,j)=(b11**2+b12**2+b13**2)/det
               baj(2,ik,i,j)=(b11*b12+b12*b22+b13*b23)/det
               baj(3,ik,i,j)=(b11*b13+b12*b23+b13*b33)/det
               baj(4,ik,i,j)=(b12**2+b22**2+b23**2)/det
               baj(5,ik,i,j)=(b12*b13+b22*b23+b23*b33)/det
               baj(6,ik,i,j)=(b13**2+b23**2+b33**2)/det
c 6*( 2 A, 4 M ) plus 1 D
            else
c this branch is not taken very often.
#ifndef NORHORG
               urg=
     &      u(1,im,jm,km)+u(1,im,jm,k)+u(1,im,jm,kp)
     &     +u(1,im,j,km)+u(1,im,j,k)+u(1,im,j,kp)
     &     +u(1,im,jp,km)+u(1,im,jp,k)+u(1,im,jp,kp)
     &     +u(1,i,jm,km)+u(1,i,jm,k)+u(1,i,jm,kp)
     &     +u(1,i,j,km)+u(1,i,j,k)+u(1,i,j,kp)
     &     +u(1,i,jp,km)+u(1,i,jp,k)+u(1,i,jp,kp)
     &     +u(1,ip,jm,km)+u(1,ip,jm,k)+u(1,ip,jm,kp)
     &     +u(1,ip,j,km)+u(1,ip,j,k)+u(1,ip,j,kp)
     &     +u(1,ip,jp,km)+u(1,ip,jp,k)+u(1,ip,jp,kp)
               urg=urg/27.
#endif
c               urg=u(1,i,j,k)
c               urg=1
c
               baj(1,ik,i,j)=urg*b11/det
               baj(2,ik,i,j)=urg*b12/det
               baj(3,ik,i,j)=urg*b13/det
               baj(4,ik,i,j)=urg*b22/det
               baj(5,ik,i,j)=urg*b23/det
               baj(6,ik,i,j)=urg*b33/det
            endif
#ifndef _SX5
         enddo
#endif
      enddo
c$omp end parallel      do
#ifdef GMETRIC
      endif
#endif
      return
      end



      subroutine rgzerosum(arr,def,dx,na1,na2,na3)
c WORKER routine.  expensive.
c set the mean of arr() to zero, by subtracting a constant multiple
c of \sqrt{g}.
c
      implicit none
      integer na1,na2,na3,k
      real arr(na1,na2,na3),dx ,def(na1,na2,na3)
cdir$ shared *def(:block,:block,:), *arr(:block,:block,:)
c locals
      integer kp,km,j,jp,jm,i,ip,im,ns1,ns2,ns3,ndim,ib,jb,kb
     &     ,io,ii,jo,ji,ko,ki,kbm,jbm,ibm
      real phixx,phiyy,phizz,phixy,phixz,phiyz,a11,a12,a13,a22,a23,a33
     &    ,di2,det,dsum,asum,dfact,csum, cfact
      include 'globalpa.fi'

#ifdef T3D
      real asum1, dsum1
cdir$ shared asum1,dsum1
      intrinsic sum
#endif
#ifdef GMETRIC
      include 'relaxgl.fi'
      include 'gmetric.fi'

      dsum=0
      if (ng1.eq.na1 .and. ng2.eq.na2 .and. ng3.eq.na3) then
c$omp parallel shared(gbaj,na1,na2,na3) private(i,j,k) reduction(+:dsum)
         do k=1,na3
c$omp do
            do j=1,na2
               do i=1,na1
                  dsum=dsum+gbaj(7,i,j,k)
               enddo
            enddo
         enddo
c$end parallel
      else
#endif    

      di2=1/dx**2
      dsum=0
#ifdef T3D
      asum=sum(arr)
#else
      csum=0
c 7.2 f77
c*$* assert do(serial)
c$omp parallel default(private) firstprivate(na1,na2,na3,di2)
c$omp& shared(def) reduction(+:dsum)      
      do k=1,na3
         kp=mod(k,na3)+1
         km=mod(k+na3-2,na3)+1
c$dir prefer_vector      
cdir$ doshared(j,i) on arr(i,j,1)      
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,na2
            do i=1,na1
            jp=mod(j,na2)+1
            jm=mod(j+na2-2,na2)+1
               ip=mod(i,na1)+1
               im=mod(i+na1-2,na1)+1
               phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))*di2
               phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))*di2
               phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))*di2
               phixy=(def(ip,jp,k)-def(im,jp,k)
     &              -def(ip,jm,k)+def(im,jm,k))*di2/4
               phiyz=(def(i,jp,kp)-def(i,jp,km)
     &              -def(i,jm,kp)+def(i,jm,km))*di2/4
               phixz=(def(ip,j,kp)-def(im,j,kp)
     &              -def(ip,j,km)+def(im,j,km))*di2/4
               a11=1+phixx
               a12=phixy
               a13=phixz
               a22=1+phiyy
               a23=phiyz
               a33=1+phizz
               det=(a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &              -a12**2*a33)
c This csum business only works if we start off with a regular
c initial grid.  The implicit assumption is that a distorted grid
c implies a non-linear region, where the exact form of the
c zero offset does not matter.  This is not true in the non-linear
c region.
c               if (det .gt. 1) csum=csum+det-1
c we use the same correction trick as in routine recut
               dsum=dsum+det
            enddo
         enddo
      enddo
c$omp end parallel
#endif
#ifdef GMETRIC
      endif
#endif
      asum=0
      di2=1/dx**2

#if 0
cdir$ master
      asum1=0
cdir$ end master
cdir$ atomic update      
      asum1=asum1+asum
c now make sure  the sum reduction finished:
cdir$ barrier
cdir$ suppress
      asum=asum1
#endif      
c      dsum=na1*na2*na3
c      if (abs(dsum-na1*na2*na3) .lt. csum) then
c         cfact=(dsum-na1*na2*na3)/csum
c      else
c         cfact=0
c      endif
      

c*$* assert do(serial)
c$omp parallel default(private) firstprivate(na1,na2,na3)
c$omp& shared(arr) reduction(+:asum)
      do k=1,na3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,na2
            do i=1,na1
               asum=asum+arr(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel
      dfact=asum/(dsum)
      csum=0
c      write(*,*)'rgzerosum: dfact=',dfact

#ifdef GMETRIC
      if (ng1.eq.na1 .and. ng2.eq.na2 .and. ng3.eq.na3) then
c*$* assert do(serial)
c$omp parallel default(private) shared(gbaj,arr)
c$omp& reduction(+:csum) firstprivate(na1,na2,na3,dfact)
         do k=1,na3
 7.2 f77
c*$* assert do(concurrent)
c$omp do
            do j=1,na2
               do i=1,na1
                  det=gbaj(7,i,j,k)
                  arr(i,j,k)=arr(i,j,k)-dfact*det
                  csum=csum+det
               enddo
            enddo
         enddo
c$omp end parallel
      else
#endif
c*$* assert do(serial)
c$omp parallel default(private) firstprivate(na1,na2,na3,di2,dfact)
c$omp& shared(arr,def) reduction(+:csum)
      do k=1,na3
         kp=mod(k,na3)+1
         km=mod(k+na3-2,na3)+1
         
c$dir prefer_vector      
cdir$ doshared(j,i) on def(i,j,1)      
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,na2
            do i=1,na1
            jp=mod(j,na2)+1
            jm=mod(j+na2-2,na2)+1
               ip=mod(i,na1)+1
               im=mod(i+na1-2,na1)+1
               phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))*di2
               phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))*di2
               phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))*di2
               phixy=(def(ip,jp,k)-def(im,jp,k)
     &              -def(ip,jm,k)+def(im,jm,k))*di2/4
               phiyz=(def(i,jp,kp)-def(i,jp,km)
     &              -def(i,jm,kp)+def(i,jm,km))*di2/4
               phixz=(def(ip,j,kp)-def(im,j,kp)
     &              -def(ip,j,km)+def(im,j,km))*di2/4
               a11=1+phixx
               a12=phixy
               a13=phixz
               a22=1+phiyy
               a23=phiyz
               a33=1+phizz
               det=(a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2
     &              -a12**2*a33)
C               if (det .gt. 1) det=det-cfact*(det-1)
               arr(i,j,k)=arr(i,j,k)-dfact*det
               csum=csum+det
            enddo
         enddo
      enddo
c$omp end parallel
#ifdef GMETRIC
      endif
#endif
      if (abs(csum-na1*na2*na3) .gt. 0.5) then
c         write(*,*) 'rgzerosum: determinant fudge failed, csum=',csum
      endif
      return
      end

      
      subroutine zerosum(arr,n1,n2,n3)
c WORKER routine
      implicit none
      integer n1,n2,n3
      real arr(n1,n2,n3)
cdir$ shared  *arr(:block,:block,:)
c locals
      integer i,j,k
      real sum1
#ifdef T3D
      intrinsic sum

      sum1=sum(arr)
#else
      sum1=0
c$omp parallel default(private) firstprivate(n1,n2,n3) shared(arr)
c$omp&    reduction(+:sum1)
c*$* assert do(serial)
      do k=1,n3
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,n2
            do i=1,n1
               sum1=sum1+arr(i,j,k)
            enddo
         enddo
      enddo
c$omp end parallel
#endif
      sum1=sum1/(n1*n2*n3)
c*$* assert do(serial)
c$omp parallel default(private) shared(arr) firstprivate(n1,n2,n3,sum1)
      do k=1,n3
cdir$ doshared(j,i) on arr(i,j,1)      
c 7.2 f77
c*$* assert do(concurrent)
c$omp do
         do j=1,n2
            do i=1,n1
               arr(i,j,k)=arr(i,j,k)-sum1
            enddo
         enddo
      enddo
c$omp end parallel
      return
      end


      subroutine szero(a,nx,ny,nz)
c WORKER routine
      implicit none
      integer nx,ny,nz
      real a(nx,ny,nz)
cdir$ shared *a(:block,:block,:)
c locals
      integer i,j,k


#ifdef T3D
      a=0
#else      
cdir$ doshared(j,i) on a(i,j,1)      
c$omp parallel private(i,j,k) firstprivate(nx,ny,nz) shared(a)
      do i=1,nx
c$omp do
         do j=1,ny
            do k=1,nz
               a(i,j,k)=0
            enddo
         enddo
      enddo
c$omp end parallel     
#endif      
      return
      end



      
c
c----------------------------------------------------------------------
c Now some self-test routines to test various aspects of the iteration.
c

      subroutine testrelax
      implicit none
c      call tst1drelax
      call test3d
      return
      end

      subroutine tst1drelax
      implicit none
      integer n1,n2,i,j,k,iopt,nrelax,it,ip,im
      parameter(n1=16,n2=1)
      real phi(n2,n2,n1), rhs(n2,n2,n1), dpot(n2,n2,n1), dx,
     &     resid(n1), rinf
     &     ,exact(n1)
c$omp parallel default(private) shared(phi,rhs,dpot,exact)
      do k=1,n1
c$omp do
         do j=1,n2
            do i=1,n2
               phi(i,j,k)=0
               rhs(i,j,k)=0
               dpot(i,j,k)=0
               exact(k)=0
            enddo
         enddo
      enddo
c$omp end parallel

      rhs(n2,n2,n1/2-1)=n1
      call zerosum(rhs,n2,n2,n1)
      iopt=0
      dx=1
      nrelax=2
      do k=1,nrelax*2
      do it=1,2
         do i=it,n1,2
            ip=mod(i,n1)+1
            im=mod(i+n1-2,n1)+1
            exact(i)=(exact(ip)+exact(im)+4*exact(i)-rhs(1,1,i))/6
         enddo
      enddo
      call zerosum(exact,n2,n2,n1)
      enddo

      call relax(phi,resid,rhs,dpot,rinf,dx,n2,n2,n1,nrelax
     &     ,iopt)
      call relax(phi,resid,rhs,dpot,rinf,dx,n2,n2,n1,nrelax
     &     ,iopt)
      call zerosum(phi,n2,n2,n1)
      write(*,'(8f10.6)')(phi(1,1,j),j=1,n1)
      write(*,*)'exact:'
      write(*,'(8f10.6)')(exact(j),j=1,n1)
      write(*,*)'real run:'
      do i=1,100
         iopt=3
         call relax(phi,resid,rhs,dpot,rinf,dx,n2,n2,n1,nrelax
     &        ,iopt)
         if (mod(i,10).eq. 0) then
            write(*,*) 'rinf=',rinf
            call zerosum(phi,n2,n2,n1)
            write(*,'(8f10.6)')(resid(j),j=1,n1)
            write(*,'(8f10.6)')(phi(1,1,j),j=1,n1)
         endif
      enddo

      iopt=2
      nrelax=0
      call relax(phi,resid,rhs,dpot,rinf,dx,n2,n2,n1,nrelax
     &        ,iopt)
      write(*,*) 'rinf=',rinf
      return
      end


      subroutine test3d
      implicit none
c Solve a 3-D periodic Poisson Problem, and compare with exact solution.
c Then measure the convergence rate, as function of nrelax, etc.
      integer n
#include "dimen.fh"
      parameter(n=NG)
      real phi(n,n,n),rhs(n,n,n),dpot(n,n,n),exact(n,n,n),dx
     &     ,rinf, resid(n,n,n), u(1,n,n,n)
      integer n1,n2,i,j,k,iopt,nrelax

c$omp parallel default(private) shared(phi,rhs,dpot,resid,u)
      do k=1,n
c$omp do
         do j=1,n
            do i=1,n
               phi(i,j,k)=0
               rhs(i,j,k)=0
               dpot(i,j,k)=0
               resid(i,j,k)=0
               u(1,i,j,k)=1
            enddo
         enddo
      enddo
c$omp end parallel
      rhs(1,1,1)=n**3
      call zerosum(rhs,n,n,n)

      iopt=3
      nrelax=4
      dx=1.
      write(*,*) 'starting relax'
      call relax(phi,resid,rhs,dpot,u,rinf,dx,n,n,n,1,nrelax
     &        ,iopt)   
      write(*,*) 'rinf=',rinf
      
      return
      end
