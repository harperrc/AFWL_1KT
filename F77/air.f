      subroutine air(eee,rrr,gmone,p,temp)

      include 'real8.h'

      save

      e     = eee * 1.0d-10
      rho1n = log(773.39520495d0 * rrr)
      e1    = (8.5d-3) * 1.0256410256d0

      if (abs(e1).lt.5.0d0)goto 20
      if (e1.le.0.0d0)goto 10

      fo  = exp(-0.22421524664d0 * e)
      fon = 0.0d0
      ws  = 1.0d0
      goto 30

 10   fo  = 0.0d0
      fon = exp(-0.15082856259d0 * e)
      ws  = 0.0
      goto 30

 20   ws  = (8.5d0 + 0.15504314d0 * rho1n - e) * 
     *       exp(-0.05d0 * rho1n + 0.02531780798d0)
      ws  = max(-60.0d0,ws,min(60.0d0,ws))
      ws  = 1.0d0 / (exp(-ws) + 1.0d0)
      fo  = exp(-0.22421524664d0 * e) * ws
      fon = exp(-0.15082856259d0 * e) * (1.0 - ws)

 30   beta = 0.0d0
      if (e.gt.1.0d0)beta = (6.9487d-3 * ws + 1.38974d-2) * log(e)
      e2  = (e - 40.0d0) * 0.33333333333333333d0
      fn  = 0.0d0

      if (abs(e2).lt.5.0d0)goto 50
      if (e2.gt.0.0d0)goto 40

      fn  = 0.0d0
      ws  = 0.0d0
      goto 60

 40   fn  = exp(-0.039215686275d0 * e)
      ws  = 1.0d0
      goto 60

 50   ws = (e - exp(0.0157d0 * rho1n + 3.806662489d0)) * 
     *  exp(-0.085d0 * rho1n - 1.38629436d0)
      ws = 1.0d0 / (exp(-ws) + 1.0d0)
      fn = exp(-0.03921568275 * e) * ws

 60   beta = beta + ws * (0.045d0 - beta)
      ws   = (e - 160.0d0) / max(30.0d0 + 3.474356d0 * rho1n,6.0d0)
      fe   = 0.0d0
      if (ws.gt.-5.0d0)fe = 1.0d0 / (exp(-ws) + 1.0d0)

      gm   = (0.161d0 + 0.225d0 * fo + 0.28d0 * fon + 
     *       0.137d0 * fn + 0.05d0 * fe) * 
     *       exp(beta * rho1n)

      gmone = gm
      p     = gmone * eee * rrr
      rho1s = rho1n * rho1n
      if (e.gt.120.0d0)goto 111

      f  =  9.720d-1 - 2.714340d-3 * rho1n + 6.582549d-5 * rho1s
      g  =  2.645d-2 - 6.418873d-4 * rho1n - 2.338785d-5 * rho1s
      h  = -9.210d-5 + 5.971549d-6 * rho1n + 3.923123d-7 * rho1s

      con1 = 3480.0d0 * gm * e
      con2 = f + g * e + h * e * e
      temp = con1 / con2
      goto 222

 111  a1r  = log10(rrr)
      c1   = 1.69081d-5 + 2.99265d-7 * a1r
      c2   = 7.69813d-1 + 3.86180d-3 * a1r
      temp = c1 * eee**c2 + 102.6275735d0

 222  beta = (e - 3.0d0) * 1.5151515151d0
      if (beta.gt.10.0d0) return
      swit = 1.0d0 / (exp(beta) + 1.0d0)
      temp = con1 / (con2 * (1.0d0 - swit) + swit)
      if (temp.gt.0.0d0) return
      a1r  = log10(rrr)
      c1   = 1.69081d-5 + 2.9265d-7 * a1r
      c2   = 7.69813d-1 + 3.8618d-3 * a1r
      temp = c1 * eee**c2 + 102.6275735d0

      return
      end
