      program main
      include 'real8.h'

      call tst02

      stop
      end
      subroutine tst02

      include 'real8.h'

      yld = 1.0e0
      hob = 0.0d0
      grm = 18.986d0 / 3.2808d0


      t   = 4.665d-4
      call shock(yld,hob,t,grm,opr,odr,vr)

      t   = t + 1.0d-6
      call shock(yld,hob,t,grm,opr,odr,vr)

      return
      end
      subroutine tst01

      include 'real8.h'

      yld   = 1.000d0
      rxft  = 18.896d0
      grm   = rxft   / 3.2808d0
      hob   = 0.0d0

c  timar1 return msec... convert to sec

      ta = timar1(rxft,hob,yld,xm) * 0.001

c  pad 

      ta = 2.00d0 * ta

      cp2psi = 0.000145038d0

      tm  = ta
      tm  = 5.00d-3

      dt  = tm / 1000.0d0
      t   = 0.0d0

c  move up past estimated time of arrival

      do while (t.le.tm)
         call shock(yld,hob,t,grm,opr,odr,vr)
         opr = opr*cp2psi
         write(10,10)t,opr,odr,vr
         t = t + dt
      enddo

      stop
 10   format(7(1x,1pe15.7))
      end
