      double precision function timar1(grft,hobft,yld,xm)

c  input
c     grft    ground range from burst (ft)
c     hobft   height of burst (ft)
c     yld     yield of detonation (Kt)

      implicit double precision (a-h,o-z)

      y13 = yld**0.3333333333d0

c  compute scaled gr,height

      x     = grft / y13
      y     = hobft / y13

c  handle 0's

      x = max(x,1.0d-9)
      y = max(y,1.0d-9)

      capr = sqrt(x * x + y * y)

      r  = capr / 1000.0d0
      z  = y / x

      if (z.gt.100.0d0)z = 100.0d0

      r2 = r * r
      r3 = r2 * r

c  compute time of shock arrival

      xm = 170.0d0 * y / (1.0d0 + 60.0d0 * y**0.25d0) + 
     *     2.89d0 * (y / 100.0d0)**2.5d0

      top = (0.543d0 - 21.8d0 * r + 386.0d0 * r2 + 
     *       2383.0d0 * r3) * r**8

      bot = 2.99d-14 - 1.91d-10 * r2 + 1.032d-6 * r**4 - 
     *      4.43d-6 * r**6 + (1.028d0 + 2.087d0 * r + 
     *      2.69d0 * r2) * r**8
      u   = top / bot

      timar1 = u

      if (x.lt.xm)return

c  if in mach region

      top = (1.086d0 - 34.605d0 * r + 486.30d0 * r2 + 
     *     2383.0d0 * r3) * r**8
  
      bot = 3.0137d-13 - 1.2128d-9 * r2 + 4.128d-6 * r**4 - 
     *     1.116d-5 * r**6 + (1.632d0 + 2.629d0 * r +
     *    2.69d0 * r2) * r**8

      w   = top / bot

      timar1 = u * xm / x + w * (1.0d0 - xm / x)

      return
      end
