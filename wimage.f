      subroutine wimage(a,n1,n2,prefix)
c assume that each time step it is first called with out, then with gas
      real*4 a(n1,n2)
      character*3 prefix
!locals
      parameter(ibar=70)
      character*1 b(n1+ibar,n2)
      integer bint(ibar,n2)
      character*9 fn
      real amin,amax,zero,one, amin2
      integer icount
      data icount /0/
      logical firsttime,flog
      data firsttime /.true./
      save amin,amax,firsttime

      n11=n2+ibar
      zero=0
      one=1
      if (firsttime) then
         firsttime=.false.
         do i=1,90
            if (i.lt.10) then
               write(fn,'(a3,1h0,i1,4h.gif)')prefix,i
            else
               write(fn,'(a3,i2,4h.gif)')prefix,i
            endif
            open(19,file=fn,err=20,status='old')
            close(19)
         enddo
         i=1
 20      continue
         icount=i-1
      endif
      if (prefix .eq. 'gas') then
         icount=icount+1
      endif
      amax=maxval(a)
      if (amax .lt. 0) then
          write(*,*) 'wimage: amax<0, using abs()'
          a=abs(a)
          amax=abs(amax)
      endif
      amin=minval(a)
      if (amin>=0) then      ! use log scale if positive
c     we want to take the range of 'a' not counting zeros.
         if (amin .le. 0) then
            amin2=2*amax
            where (a .le. 0) a=amin2
            amin=minval(a)
            where (a .eq. amin2) a=amin
         endif
         amin=max(amax/10000,amin)
         a=log10(max(amin,a))
         amax=log10(amax)
         amin=log10(amin)
      endif
      b(:n1,:)=char(int(255*min(one,max(zero,
     &              (a-amin)/(amax-amin)))))
      if (icount.lt.10) then
         write(fn,'(a3,1h0,i1,4h.gif)')prefix,icount
      else
         write(fn,'(a3,i2,4h.gif)')prefix,icount
      endif
      flog=amin>=0
      call label(bint,ibar,n2,amin,amax,flog)
      b(n1+1:,:)=char(bint)
      write(*,*) 'wimage: file=',fn,' aminmax=',amin,amax
      call wgif(b,n11,n2,fn)
      return
      end

      subroutine label(bint,n1,n2,x,x2,flog)
      integer bint(n1,n2)
      parameter(nx=5,ny=9,ioff1=10,nlabel=7)
      integer template(7,nx,ny), digits(12,nx,ny)
      logical flog
      integer segments(7,12)
      data segments /1,1,1,0,1,1,1,   ! 0
     &               0,0,1,0,0,1,0,
     &               1,0,1,1,1,0,1,
     &               1,0,1,1,0,1,1,
     &               0,1,1,1,0,1,0,   ! 4
     &               1,1,0,1,0,1,1,
     &               1,1,0,1,1,1,1,
     &               1,0,1,0,0,1,0,   ! 7
     &               1,1,1,1,1,1,1,
     &               1,1,1,1,0,1,1,
     &               0,0,0,1,0,0,0,   ! -
     &               1,1,0,1,1,0,1/   ! E

      joffset=1
      template=0
      digits=0
      do i=1,nx
         template(1,i,1)=1
         template(4,i,(ny+1)/2)=1
         template(7,i,ny)=1
      enddo
      do i=1,(ny+1)/2
         template(2,1,i)=1
         template(3,nx,i)=1
      enddo
      do i=(ny+1)/2,ny
         template(5,1,i)=1
         template(6,nx,i)=1
      enddo
      
      bint=0
c$omp parallel private(i,j) shared(n2,bint)
      do i=1,ioff1
c$omp do
         do j=1,n2
            bint(i,j)=(j-1)*255./n2
         enddo
      enddo
c$omp end parallel
      ioffset=ioff1

      
      do ii=1,nlabel
         w2=(ii-1.)/(nlabel-1.)
         w1=1-w2
         joffset=w1+w2*(n2-ny-1)
         if (flog) then
            xt=10**(w1*x+w2*x2)
         else
            xt=(w1*x+w2*x2)
         endif
!         write(*,*) xt,joffset,w1,w2
         call plabel(xt,ioffset)
      enddo
      contains


      subroutine plabel(x,ioff1)
      integer ioffset
      ioffset=ioff1
      xsign=sign(1.,x)
      x=abs(x)
      if (x.ne.0.0) then 
       exponent=floor(log10(x))+1
      else
       exponent=0     
      endif
      mantissa=x/10**exponent*1000+0.5
      if (xsign<0) then
c         write(*,*) 'adding -'
         do j=1,nx
            bint(j+ioffset,(ny+1)/2+joffset)=(255)
         enddo
         ioffset=ioffset+nx+1
      endif
      bint(ioffset+nx/2,ny+joffset)=255
      bint(ioffset+nx/2,ny+1+joffset)=255
      bint(ioffset+nx/2+1,ny+joffset)=255
      bint(ioffset+nx/2+1,ny+1+joffset)=255
      ioffset=ioffset+4
      do i=2,0,-1
         idigit=mantissa/10**i
!      write(*,*) 'mantissa=',mantissa,idigit
      mantissa=mantissa-idigit*10**i
      idigit=idigit+1
      call pchar(idigit,ioffset)
      enddo
      ioffset=ioffset+nx/2
      call pchar(12,ioffset)
      ioffset=ioffset+nx/2
      xsign=sign(1.,exponent)
      exponent=abs(exponent)
      if (exponent .gt. 0.5) then
         iexp=log10(exponent)
      else
         iexp=0
      endif
      if (xsign<0) then
c$omp parallel private(j) shared(bint,ioffset,joffset)
c$omp do
         do j=1,nx
            bint(j+ioffset,(ny+1)/2+joffset)=(255)
         enddo
c$omp end parallel
         ioffset=ioffset+nx+1
      endif
      mantissa=exponent
      do i=iexp,0,-1
         idigit=mantissa/10**i
         mantissa=mantissa-idigit*10**i
         idigit=idigit+1
         call pchar(idigit,ioffset)
      enddo
      end subroutine
      
      subroutine pchar(idigit,ioffset)
c$omp parallel private(i,i1,j1)
c$omp&      shared(segments,idigit,ioffset,joffset,bints,template)
c$omp do 
      do j=1,7
         if (segments(j,idigit) .gt. 0) then
            do i1=1,nx
               do j1=1,ny
                  bint(i1+ioffset,j1+joffset)=max(bint(i1+ioffset
     &                 ,j1+joffset),(255*template(j,i1,j1)))
               enddo
            enddo
         endif
      enddo
c$omp end parallel
      ioffset=ioffset+nx+1
      end subroutine
      end
      



