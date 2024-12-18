      function speed(t)

      include 'real8.h'

      save

      r1 = 24210.0d0 * t**0.371d0 * (1.0d0 + (1.23d0 * t + 0.123d0) * 
     *     (1.0d0 - exp(-26.25d0 * t**0.79d0)))

      r2 = (1.0d0 - 0.03291d0 * t**(-1.086d0)) *
     *     (33897.0d0 * t + 8490.0d0) + 8.36d3 + 2.5d3 *
     *     log(t) + 800.0d0 * t**(-0.21d0)

      v1 = 10086.685d0 * t**(-0.629d0) + 40826.049d0 * t**0.371d0 +
     *  exp(-26.25d0 * t**0.79d0) *
     * (617527.5d0 * t**1.161d0 - 40826.048d0 * t**0.371d0 -
     *  1104.7749d0 * t**(-0.629d0) + 
     *  61752.75d0 * t**0.161d0)

      v2 = 33897.0d0 + 95.937326d0 * t**(-1.086d0) + 
     *     303.43481d0 * t**(-2.086d0) + 2.5d3 / t - 
     *     168.0d0 * t**(-1.21d0)

      v3 = (v2 * (t - 0.21d0) + v1 * (0.28d0 - t)) / 0.07d0

      if (t.lt.0.21d0) then
         speed = v1
      else if (t.gt.0.28d0) then
         speed = v2
      else
         speed = v3
      endif

      return
      end
