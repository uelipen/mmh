       program project1
       
c 8/10/99 power spectrum code added back in.
c
c Contains fix for planes with reverse a values (no negative tau's).
c
c
c 7/22/99 power spectrum code removed for Origin compiling
c both in program, and subroutine/functions.
c
c projects the simulation outputs using the weights according 
c to the field, eg SZ, kinetic SZ, mass density... Can be used as
c a weak lensing map without displacement corrections

        implicit real(a-h,o-z)
        real tol
        parameter (tol=1.0d-6,nx=1000,twopi=6.2831853,nmax=4096)
        parameter (nm=1024,nrm=2*nm,nxyz=1,nplane=500,nreamax=5)
        real xtau(nx),xa(nx),xpr(nx)
        real func1,a0,at,a1,a2,a3,xaold,xaromb,apmtemp
        real xaplo,xaphi,xtauold
        real omegam,omegav,h0,curv,tau0,omegabh
        real dx,omegampm,omegavpm,h0pm
        real kappapm(nm,nm,nxyz),a(nplane)
        real apm,size(nplane),kappa(nrm+2,nrm,nreamax)
        real kappaold(nrm+2,nrm,nreamax)
        real dchipl(nplane),chipl(nplane),atau(nplane)
        real*8 deltag2(nmax)
        common /cosmos/ omegam,omegav,h0,curv,tau0
        external func1
        character*80 filename
        character*3 prefix
        character*9 gifname
        prefix='img'
c  Read projected fields
        write(*,*) 'Enter input filename'
        read(*,'(a)') filename
        write(*,*)'Enter size rescale:'
        read(*,*)istretch
        open(10,file=filename,status='old',form='unformatted')
	write(*,*) 'Opened',filename
        rewind 10
	write(*,*) 'Rewind',filename
        read(10)n1,n2,n3,nimage
        write(*,*) 'Read n1,n2,n3,nimage from ', filename
	read(10)dx,omegampm,omegavpm,h0pm,sigma8pm
	write(*,*)'omegam,omegav,h0,sigma8pm'
     1  ,omegampm,omegavpm,h0pm,sigma8pm
	write(*,*) 'Read parameters from',filename
	write(*,*) 'dx=',dx,'n1=',n1
        write(*,*)'sigma8,omega_bh'
        read(*,*)sigma8,omegabh
        deltag2(:nmax)=0.0d0
c Inserted from later in the code        
        write(*,*)'Insert integer seed, number of realizations'
        read(*,*)idum,nre
        write(*,*) 'Enter output power spectrum filename'
        read(*,'(a)') filename
        open(13,file=filename,status='unknown',form='formatted')
        write(13,*)nre,n1
        write(13,*)dx*n1,omegampm,omegavpm,h0pm
        write(13,*)sigma8pm,sigma8,omegabh
        write(*,*) 'Enter output map filename'
        read(*,'(a)') filename  
c choose between formatted or unformatted output
        open(22,file=filename,status='unknown',form='unformatted',
     1  convert='little_endian')
c        open(22,file=filename,status='unknown',form='formatted') 
        write(*,*) 'Enter output gif filename'
        read(*,'(a)') gifname
        
        n1r=n1*istretch
        if (n1.ne.n2.or.n1.ne.n3) pause 'Only cubic box allowed'
        q=sigma8/sigma8pm
        if (q.gt.1.0) pause 'sigma8 larger than in the simulation'
c rescale omega's
        curvpm=omegampm+omegavpm-1.0
        curv2=curvpm/omegampm*q
        omegav2=omegavpm/omegampm*q*q*q
        omegam=1.0/(1.0+omegav2-curv2)
        curv1=curv2*omegam
        omegav=omegav2*omegam
        h0=h0pm*100.0
        write(*,*)n1*dx,omegam,omegav,curv1,h0

        write(*,*)'Maximum redshift of galaxies'
c altered for multiple means
        read(*,*)z
        a1=1.0/(1.0+z)
        hc=2.998e5/h0
        a0=0.
        at=1.0
        tau0=hc*rombint(func1,a0,at,tol)
        tau1=hc*rombint(func1,a0,a1,tol)
        curv=curv1/hc/hc
        write(*,*)'galaxy distance/horizon ratio:',(1-tau1/tau0)
c set up tau,a table
        xaold=a1
        xtauold=tau1
        do 15 i=1,nx
         xa(i)=xaold+(at-a1)/nx
         xtau(i)=xtauold+hc*rombint(func1,xaold,xa(i),tol)
         xaold=xa(i)
         xtauold=xtau(i)
15      continue
        xaplo=1.0e32
        xaphi=1.0e32
	write(*,*)'starting spline'
        call spline(xtau,xa,nx,xaplo,xaphi,xpr)
	write(*,*)'spline finished'
        chi=0.0
        ipl=1

c seed must be negative
        if (idum.gt.0.0) idum=-idum
        write(*,*)'Which map [0-4]:'
        read(*,*)ifact
        write(*,*)'ifact=',ifact
c initialize
         kappa=0.0

40      continue
c transverse size
        size(ipl)=n1r*dx
c radial size
        dchipl(ipl)=dx*n1
c advance by one plane in time
        chi=chi+dx*n1
        chipl(ipl)=chi
c find corresponding a in the middle of the plane
        tau=tau0-(chi-dx*n1/2)
c	write(*,*)'ipl=',ipl
        call splint(xtau,xa,xpr,nx,tau,a(ipl))
c find corresponding a in the simulation
c make a somewhat smaller for inequality tests
       atau(ipl)=a(ipl)*q-1.0e-5
       ipl=ipl+1
       if (tau.gt.tau1) goto 40
       write(*,*)'passed goto 40'
       npl=ipl-1
       read(10,end=10)apm
       l=1
       apmtemp=apm
       apmax=0.0d0
       tau=hc*rombint(func1,a0,apmtemp,tol)
       theta0=size(npl)/r(chipl(npl))

c begin loop over planes
       write(*,*) 'Beginning loop over planes'
       jjj=1
       do ipl=npl,1,-1
c        ixyz=mod(ipl-1,3)+1
         ixyz=1
c read input
20      continue
        if ((atau(ipl).gt.apm).or.(apm.lt.apmax)) then
          read(10)((kappapm(k,j,1), k=1,n1), j=1,n1)
          amax=maxval(kappapm(:n1,:n1,1))
          amin=minval(kappapm(:n1,:n1,1))
	  s=sum(kappapm(:n1,:n1,1))
          s=s/n1/n1
        write(*,'(1I4,1E12.3,1I4,5E12.3)')l,apm,
     1  ipl,atau(ipl),tau,amax,amin,s
          l=l+1
          apmold=apm
          if (apm.gt.apmax) apmax=apm
          read(10,end=10)apm
          tau=hc*rombint(func1,apmold,apm,tol)
          goto 20
        else
         if (l.eq.1) pause 'initial z too high'


c got the plane, now rescale and add with appropriate weight 

c weights: 
c weak lensing
          if (ifact.eq.0) then
           weight=r(chipl(ipl)-dchipl(ipl)/2.0)/r(chipl(npl))*
     1     r(chipl(npl)+dchipl(npl)/2.0-chipl(ipl))
           fact=1.5/hc**2*omegam*dchipl(ipl)/a(ipl)*weight
c           xm=0.0
c           xm=sum(kappapm(:n1,:n1,1))
c           xm=xm/n1/n1
c           kappapm(:n1,:n1,1)=kappapm(:n1,:n1,1)-xm
          endif

c thermal SZ (velocity squared rescaling gives q)
           if (ifact.eq.1) fact=omegabh/a(ipl)**2*dchipl(ipl)*q
c kinetic SZ (velocity rescaling)
           if (ifact.eq.2) fact=omegabh/a(ipl)**2*dchipl(ipl)*sqrt(q)
c Rees-Sciama (dphi/da passed)
           if (ifact.eq.3) fact=dchipl(ipl)/func1(a(ipl))/hc*q*q
c bolometric X-ray, sqrt(T) gives sqrt(q)
      if (ifact.eq.4) fact=dchipl(ipl)*omegabh**2/a(ipl)**4*sqrt(q)

c rescale the array by scale
           scale=r(chipl(ipl)-dchipl(ipl)/2.0)/r(chipl(npl)-
     1     dchipl(npl)/2.0)*size(ipl)/size(npl)

c test
c           kappapm(:n1,:n1,ixyz)=kappapm(:n1,:n1,ixyz)+1.0
c loop over realizations
          do ir=1,nre
            kappaold=kappa
c randomly choose the center of the plane
            i1=ran1(idum)*n1
            j1=ran1(idum)*n1
c positions of photons
           do i=1,n1r
c unperturbed position
            x1=i*scale+i1
            do j=1,n1r
              y1=j*scale+j1
c find position in the grid
              ix=x1
              iy=y1
              wx=x1-ix
              wy=y1-iy
              ix = mod(ix+n1-1,n1)+1
              iy = mod(iy+n1-1,n1)+1
              w11=wx*wy
              w10=wx*(1-wy)
              w01=(1-wx)*wy
              w00=(1-wx)*(1-wy)
              ixp=mod(ix,n1)+1
              iyp=mod(iy,n1)+1
          kappa(j,i,ir)=kappa(j,i,ir)+w00*kappapm(iy,ix,ixyz)*fact
          kappa(j,i,ir)=kappa(j,i,ir)+w01*kappapm(iyp,ix,ixyz)*fact
          kappa(j,i,ir)=kappa(j,i,ir)+w10*kappapm(iy,ixp,ixyz)*fact
          kappa(j,i,ir)=kappa(j,i,ir)+w11*kappapm(iyp,ixp,ixyz)*fact
            enddo
          enddo
         enddo
        endif
       kappaold=kappa-kappaold
c       call wimage(kappaold(1,1,1),n1r,nrm,prefix,istretch)
	s=sum(kappa(:n1r,:n1r,1))
        amax=maxval(kappa(:n1r,:n1r,1))
        amin=minval(kappa(:n1r,:n1r,1))
        write(*,*)s/n1r/n1r,amax,amin
c        if (mod(ipl,5).eq.0) then
c        write(20+jjj,*)1.0/apm-1,theta0*360.0/twopi
c        call ps1(kappa(1,1,1),theta0,n1r,nrm,jjj,apm)
c        jjj=jjj+1
c        endif
       enddo
10     continue
       istretch=n1r/512
       write(*,*)'istretch=',istretch
       call wimage1(kappa(1,1,1),n1r,nrm,gifname,istretch)
       close(10)
       l=l-1
       write(*,*)'Number of input planes and final a is',l,apm
c compute power spectrum 
         write(*,*)'angle in degrees:',theta0*360.0/twopi
         do ir=1,nre
	  write(*,*) 'call ps'
          call ps(kappa(1,1,ir),theta0,deltag2,n1r,nrm)
c          write(*,*)'write kappa. n1r=',n1r,'ir=',ir 
c 	  write(22)kappa(:n1r,:n1r,ir)
          write(22)((kappa(k,j,1), k=1,n1), j=1,n1)
          write(*,*)(kappa(k,1,1),k=1,10)
          write(33,*)(kappa(k,1,1),k=1,10)
          s=s+sum(kappa(:n1r,:n1r,ir))
          s2=s2+sum(kappa(:n1r,:n1r,ir)**2)
         enddo
          s=s/n1r/n1r/2/nre
          s2=s2/n1r/n1r/4/nre
          rms=sqrt(s2-s*s)
          write(*,*)'mean y,rms(y)=',s,rms
	  write(*,*)'Initial z was z=',z
          write(13,*)theta0*360.0/twopi,s,rms,z
         dk=twopi/theta0
          do i=1,n1r
            ak=i*dk
	   if(deltag2(i).gt.0) write(13,'(3E13.5)') 
     1     ak,deltag2(i)/nre
          enddo
      stop
      end

       function func1(a)
        implicit real (a-h,o-z)
        common /cosmos/ omegam,omegav,h0,curv,tau0
        real omegam,omegav,h0,curv,tau0
c      integrates the a(tau) relation
        parameter (tcmb=2.726,nnur=3)
        parameter (grhog=1.4952d-13*tcmb**4)
        parameter (grhor=3.3957d-14*tcmb**4)
        grhom=3.3379d-11*h0*h0
        aeq=(grhog+grhor*nnur)/grhom
        func1=1./sqrt(a*omegam+aeq+a*a*(1.-omegam-omegav)+
     2  a*a*a*a*omegav)
       return
       end

        function rombint(f,a,b,tol)
c  Rombint returns the integral from a to b of using Romberg integration.
c  The method converges provided that f(x) is continuous in (a,b).
c  f must be double precision and must be declared external in the calling
c  routine.  tol indicates the desired relative accuracy in the integral.
c
        parameter (MAXITER=16,MAXJ=5)
        implicit real (a-h,o-z)
        dimension g(MAXJ+1)
        real f
        external f
c
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol))
     2      go to 40
c  Calculate next trapezoidal rule approximation to integral.
          g0=0.0d0
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1.0d0
            do 30 j=1,jmax
c  Use Richardson extrapolation.
            fourj=4.0d0*fourj
            g1=g0+(g0-g(j))/(fourj-1.0d0)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1.0d0-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)
     2    write(*,*) 'Rombint failed to converge; integral, error=',
     3    rombint,error
        return
        end

        function r(chi)
         implicit real(a-h,o-z)
         common /cosmos/ omegam,omegav,h0,curv,tau0
         real omegam,omegav,h0,curv,tau0
         r=chi
         if (curv.ne.0.0d0) then
         if (curv.gt.0.) then
          r=sin(sqrt(curv)*r)/sqrt(curv)
         else
          r=sinh(sqrt(-curv)*r)/sqrt(-curv)
         endif
         endif
         if (r.lt.0.0d0) r=0.0d0
        return
        end

      subroutine fft2(a,n,nm,isign)
       implicit real*4 (a-h,o-z)
c 2-d real-to-complex fft with normalization at the R to C transform
c isign=1 means forward transform
c isign=-1 means inverse transform
      real*4 a(nm+2,nm)
      integer sfft_2d, dfft_2d
      external sfft_2d, dfft_2d
      n2=n*n
      xn2inv=1.0/n2
      if (isign .eq. 1) then
         ierr=sfft_2d('R','C','F',a,a,n,n,nm+2,1,1)
      else
         ierr=sfft_2d('C','R','B',a,a,n,n,nm+2,1,1)
      endif
c reverse the normalization
      do j=1,n
       do i=1,n+2
        if (isign .eq. 1) then
         a(i,j)=a(i,j)*xn2inv
        else
         a(i,j)=a(i,j)*n2
        endif
       enddo
      enddo
      ierr=0
      if (ierr .ne. 0) then
                write(*,*) 'forward fft failed, ier=',ierr
                write(*,*) 'isign,n=',isign,n
                stop
      endif
      return
      end

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      real yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=1000)
      INTEGER i,k
      real p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
      END

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      real x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      real a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0
      return
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
      END

       subroutine ps1(kappag,size,n1,nm,ind,apm)
        implicit real*4 (a-h,o-z)
        parameter (twopi=6.2831853,nmax=4096)
        parameter (sqr2=1.41421356)
        real*4 kappag(nm+2,nm)
        real*4 powg(nmax)
        real*8 deltag2(nmax)
        integer np(nmax) 
        dk=twopi/size
        d2k=dk**2
        akmax=n1*dk
        n12=n1/2
c initialize power spectrum arrays
        do i=1,n1
         powg(i)=0.0
         np(i)=0
        enddo 
c fft gamma to k space: 
        call fft2(kappag,n1,nm,1)
        do k2=1,n1
          ak2=(k2-1)*dk
          if (k2.gt.n12) ak2=ak2-akmax
          do k1=1,n12+1
            ak1=(k1-1)*dk
            akk=ak1*ak1+ak2*ak2
              ak=sqrt(akk)
c correct for k_n12+1=+/-k_Nyquist
              ka=ak/dk+0.5
              dpowg=(kappag(2*k1-1,k2)**2+kappag(2*k1,k2)**2)
c  Count twice for conjugate harmonic.
             if ((k1.ne.1.and.k1.ne.(n12+1)).or.(k2.eq.1).or.
     1       (k2.eq.(n12+1))) then
               powg(ka)=powg(ka)+2.0d0*dpowg
               np(ka)=np(ka)+2
             else
               powg(ka)=powg(ka)+dpowg
               np(ka)=np(ka)+1
             endif
          enddo
         enddo
c  Normalize power per mode.
          do i=1,n1
           if (np(i).gt.0) fact=1.0d0/np(i)
            ak=i*dk
            powg(i)=powg(i)*fact
            write(20+ind,*)ak,powg(i)*i*i*twopi
          enddo
          call fft2(kappag,n1,nm,-1)
          return
         end

       subroutine ps(kappag,size,deltag2,n1,nm)
        implicit real*4 (a-h,o-z)
        parameter (twopi=6.2831853,nmax=4096)
        parameter (sqr2=1.41421356)
        real*4 kappag(nm+2,nm)
        real*4 powg(nmax)
        real*8 deltag2(nmax)
        integer np(nmax) 
        dk=twopi/size
        d2k=dk**2
        akmax=n1*dk
        n12=n1/2
c initialize power spectrum arrays
        do i=1,n1
         powg(i)=0.0
         np(i)=0
        enddo 
c fft gamma to k space: 
        call fft2(kappag,n1,nm,1)
        do k2=1,n1
          ak2=(k2-1)*dk
          if (k2.gt.n12) ak2=ak2-akmax
          do k1=1,n12+1
            ak1=(k1-1)*dk
            akk=ak1*ak1+ak2*ak2
              ak=sqrt(akk)
c correct for k_n12+1=+/-k_Nyquist
              ka=ak/dk+0.5
              dpowg=(kappag(2*k1-1,k2)**2+kappag(2*k1,k2)**2)
c  Count twice for conjugate harmonic.
             if ((k1.ne.1.and.k1.ne.(n12+1)).or.(k2.eq.1).or.
     1       (k2.eq.(n12+1))) then
               powg(ka)=powg(ka)+2.0d0*dpowg
               np(ka)=np(ka)+2
             else
               powg(ka)=powg(ka)+dpowg
               np(ka)=np(ka)+1
             endif
          enddo
         enddo
c  Normalize power per mode.
          do i=1,n1
           if (np(i).gt.0) fact=1.0d0/np(i)
            ak=i*dk
            powg(i)=powg(i)*fact
            deltag2(i)=deltag2(i)+powg(i)*i*i*twopi
          enddo
          call fft2(kappag,n1,nm,-1)
          return
         end

      FUNCTION ran1(idum)
       implicit real*4 (a-h,o-z)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*4 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      end

      subroutine wimage1(a,n1,n2,fn,istretch)
c assume that each time step it is first called with out, then with gas
      real*4 a(n2+2,n2)
!locals
      parameter(ibar=70)
      character*1 b(n1+istretch*ibar,n1)
      integer bint(ibar,n1/istretch)
      character*9 fn
      real amin,amax,zero,one
      integer icount
      data icount /0/
      logical firsttime
      data firsttime /.true./
      save amin,amax,firsttime

      n11=n1+ibar*istretch
      zero=0
      one=1
      if (firsttime) then
         firsttime=.false.
         i=1
      endif
      amax=maxval(a(:n1,:n1))
      if (amax .lt. 0) then
          write(*,*) 'wimage: amax<0, using abs()'
          a=abs(a)
          amax=abs(amax)
      endif
c      a=log10(max(amax/10000,a))
      amin=minval(a(:n1,:n1))
c      amax=log10(amax)
      b(:n1,:n1)=char(int(255*min(one,max(zero,
     &              (a(:n1,:n1)-amin)/(amax-amin)))))
      call label(bint,ibar,n1/istretch,amin,amax)
      do i=0,n1/istretch-1
      do j=0,ibar-1
      do k=1,istretch
      do l=1,istretch
      b(n1+j*istretch+k,i*istretch+l)=char(bint(j+1,i+1))
      enddo
      enddo
      enddo
      enddo
      write(*,*) 'wimage: file=',fn,' aminmax=',amin,amax
      call wgif(b,n11,n1,fn)
      return
      end

      subroutine wimage(a,n1,n2,prefix,Istretch)
c assume that each time step it is first called with out, then with gas
      real*4 a(n2+2,n2)
      character*3 prefix
!locals
      parameter(ibar=70)
      character*1 b(n1+istretch*ibar,n1)
      integer bint(ibar,n1/istretch)
      character*9 fn
      real amin,amax,zero,one
      integer icount
      data icount /0/
      logical firsttime
      data firsttime /.true./
      save amin,amax,firsttime
      n11=n1+ibar*istretch
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
            write(*,*)fn
         enddo
         i=1
 20      continue
         icount=i-1
      endif
      if (prefix .eq. 'img') then
         icount=icount+1
      endif
c      a(:n1,:n1)=log(4.0*a(:n1,:n1)+1.0)
      amax=maxval(a(:n1,:n1))
      if (amax .lt. 0) then
          write(*,*) 'wimage: amax<0, using abs()'
          a=abs(a)
          amax=abs(amax)
      endif
      amin=minval(a(:n1,:n1))
c increase the contrast
c      a=a-amin
c      a=log10(max(amax/10000,a))
c      amax=log10(amax)
c      amin=log10(amax/10000)
      b(:n1,:)=char(int(255*min(one,max(zero,
     &              (a(:n1,:n1)-amin)/(amax-amin)))))
      if (icount.lt.10) then
         write(fn,'(a3,1h0,i1,4h.gif)')prefix,icount
      else
         write(fn,'(a3,i2,4h.gif)')prefix,icount
      endif
      call label(bint,ibar,n1/istretch,amin,amax)
      do i=0,n1/istretch-1
      do j=0,ibar-1
      do k=1,istretch
      do l=1,istretch
      b(n1+j*istretch+k,i*istretch+l)=char(bint(j+1,i+1))
      enddo
      enddo
      enddo
      enddo
      write(*,*) 'wimage: file=',fn,' aminmax=',amin,amax
      call wgif(b,n11,n1,fn)
      return
      end

      subroutine label(bint,n1,n2,x,x2)
      integer bint(n1,n2)
      parameter(nx=5,ny=9,ioff1=10,nlabel=7)
      integer template(7,nx,ny), digits(12,nx,ny)
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
      do i=1,ioff1
         do j=1,n2
            bint(i,j)=(j-1)*255./n2
         enddo
      enddo
      ioffset=ioff1

      
      do ii=1,nlabel
         w2=(ii-1.)/(nlabel-1.)
         w1=1-w2
         joffset=w1+w2*(n2-ny-1)
         xt=w1*x+w2*x2
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
         do j=1,nx
            bint(j+ioffset,(ny+1)/2+joffset)=(255)
         enddo
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
      ioffset=ioffset+nx+1
      end subroutine
      end
      

