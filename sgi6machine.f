

      subroutine adjustthread
c adjust the number of threads used such that it does not exceed the
c system load.
c we need to give it enough time that the load actually reflects
c the current system load.  This is about one hour.
      implicit none
      include 'relaxgl.fi'
      real*8 sysload 
      integer lastcall, minthread, maxthread, shortcall, newthread
      integer ncpu
      parameter (minthread=2)
      logical firsttime
      save lastcall, firsttime, maxthread, shortcall, ncpu
      data firsttime /.true./
      integer numthreads, idealthread, mp_numthreads, time
      external time, mp_numthreads

      if (firsttime) then
         firsttime = .false.
         ncpu=14
         call cncpu(ncpu)
         write(*,*) 'adjustthread: we think machine has ', ncpu,' CPUs'
         maxthread=ncpu-1
         if (ng1 .le. 64) maxthread=8
         lastcall=time()
         shortcall=time()
         return
      endif
c
c check for manual thread adjustment by reading file 'adjthread.dat'.
c Set the number of threads to its contents, and reset the periodic time
c checking.  If the number of threads is negative, set to the absolute value,
c and don't check the load for the next 1000 hours.
c
c if the file does not exist, do not rethread, but keep checking for the file
      if (time() - shortcall .gt. 60 ) then
         shortcall=time()
         open(32,file='adjthread.dat',status='old',err=13)
         read(32,*,err=25) newthread
         rewind(32,err=35)
         write(32,*,err=45) 'adjthread: DONE READING'
         close(32,err=65)
         write(*,*) 'adjthread: using manual file, newthread=',newthread
         if (newthread .lt. 0) then
            newthread=-newthread
            lastcall=time()+3600*1000
         else
            lastcall=time()
         endif
         call mp_destroy
         call mp_set_numthreads(newthread)
         return
13       continue
         lastcall=time()+3600*10000
         return
25       continue
         close(32,err=65)
         goto 15
35       continue
         write(*,*) 'adjthread: rewind failed'
         goto 15
45       continue
         write(*,*) 'adjthread: write failed'
         goto 15
65       continue
         write(*,*) 'adjthread: close(32) failed!!'
15       continue
      endif
 
c         
c Once per hour, attempt to set the system load such that the total load leaves
c one CPU empty.
c
      if (time() - lastcall .lt. 3600 ) return
      lastcall=time()
      numthreads=mp_numthreads()
      call csysload(sysload)
      write(*,*) 'adjthread: load=',sysload
      if (sysload .le. 1 .or. sysload .gt. 100) return
c      idealthread=14.4-sysload+numthreads
      idealthread=ncpu-0.1-sysload+numthreads
      idealthread=min(maxthread,idealthread)
      idealthread=max(minthread,idealthread)
      if (idealthread .ne. numthreads) then
         write(*,55) sysload,idealthread
   55      format('adjthread: load=',f7.2,', readjusting to ',i2
     &               ,' threads')
         call mp_destroy
         call mp_set_numthreads(idealthread)
      endif


      return
      end
