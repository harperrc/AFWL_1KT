      function wfpr(t)

      include 'real8.h'

      save

      r1 = 24210.0d0 * t**0.371d0 *
     *    (1.0d0 + (1.23d0 * t + 0.123d0) * 
     *    (1.0d0 - exp(-26.25d0 * t**0.79d0)))

      r2 = (1.0d0 - 0.03291d0 * t**(-1.086d0)) * 
     *    (33897.0d0 * t + 8490.0d0) + 
     *    8.36d3 + 2.5d3 *
     *    log(t) + 800.0d0 * t**(-0.21d0)

      if (t.lt.0.21d0) then
         wfpr = r1
      else if (t.gt.0.28d0) then
         wfpr = r2
      else
         wfpr = (r2 * (t - 0.21d0) + r1 * (0.28d0 - t)) / 0.07d0
      endif

      return
      end
