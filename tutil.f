
c*******************************************************************
c
c      
c     we start with the Friedman equation
c
c   [(da/dt)^2+k']/a^2 = 8\pi G \rho/3 + \Lambda'
c and substitute dt = a^2 d\tau, 8\pi G=4/3, \rho=1/a^3 (in our units)
c      
c    (da/d\tau)^2 = (4/9) a^3 [ 1 + a^3 \Lambda - a k ]
c      
c which we then normalize at a=1: let \Omega_0 be the total non-relativistic
c dark matter, then
c      
c    \Lambda=\Omega_L/\Omega_0,    k=(\Omega_0+\Omega_L-1)/\Omega_0      
c
c  this constraint is necessary since (?)      
c
c   When \Omega_L=0, the exact solution is \sqrt{1/a-k} = -\tau/3+c
c   For a closed universe, this coordinate system is symmetric around
c   \tau=0, and runs from -\infty < \tau < \infty.  For an open universe,
c   the system ends at some finite \tau.      
c      
c   Unfortunately, unlike in the proper time case, the \Omega_L case is
c   not solvable.      
c 

      function dadtau(a)
      implicit none
      include 'globalpa.fi'
      real dadtau, a,ak

      ak=(omega0+omegal-1)/omega0
      dadtau=2*sqrt(a**3)/3*sqrt(1 + a**3*(omegal)/omega0 - a*ak)
      return
      end


      function dascale(t,a,dtau)
c assume that a_f=1, and that omega0 is given at a_f.      
      implicit none
      real a, dtau, dascale, t
c t is passed back to the calling program
c assume universe with only lambda and matter      
c locals
      real adot, addot, atdot, anew, asinh, x, dadtau, rc, ak, oterm
     &     , otd
      external dadtau
      include 'globalpa.fi'
      logical firsttime
      data firsttime /.true./
      asinh(x)=log(x+sqrt(x**2+1))

      if (firsttime) then
         firsttime = .false.
         write(*,*) 'dascale: Omega_0,L=',omega0,omegal
      endif

      adot=dadtau(a)


      ak=(omega0+omegal-1)/omega0
      rc=(omegal)/omega0
      oterm=1+a**3*rc-a*ak
      otd=3*a**2*rc-ak
      
      
      addot=2./3.*a**2 + 4./3.*a**5*rc - 8./9.*a**3*ak
      atdot=adot*(4./3.*a + 20./3.*a**4*rc - 8./3.*a**2*ak)

      dascale=adot*dtau+addot*dtau**2/2+atdot*dtau**3/6

      anew=a+dascale

      rc=max(1.e-8,rc)
      t=asinh(sqrt(anew**3*rc))/sqrt(rc)

      return
      end
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc
        function tauconf(a,omegam,omegav)
c  Evaluate tau(a) (inverse of a(tau)) for FLRW cosmology.
c  Omegam := Omega today (a=1) in matter.
c  Omegav := Omega today (a=1) in vacuum energy.
c  dtau := H0*dt/a^2.
        real tau,a,omegam,omegav
        double precision om,ov,adp,rombint,dtaudaconf
        common /omegas/ om,ov
        external dtaudaconf
c
        om=omegam
        ov=omegav
        adp=a
        tauconf=rombint(dtaudaconf,1.0d0,adp,1.0d-8)
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc
        function dtaudaconf(a)
        implicit double precision (a-h,o-z)
        common /omegas/ omegam,omegav
c
        eta=sqrt(omegam/a+omegav*a*a+1.0d0-omegam-omegav)
        dtaudaconf=1.0d0/(a*eta)
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
        function rombint(f,a,b,tol)
c  Rombint returns the integral from a to b of f(x)dx using Romberg integration.
c  The method converges provided that f(x) is continuous in (a,b).  The function
c  f must be double precision and must be declared external in the calling
c  routine.  tol indicates the desired relative accuracy in the integral.
c
        parameter (MAXITER=16,MAXJ=5)
        implicit double precision (a-h,o-z)
        dimension g(MAXJ+1)
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
     2       write(*,*) 'Rombint failed to converge; integral, error=',
     3       rombint,error
        return
        end

