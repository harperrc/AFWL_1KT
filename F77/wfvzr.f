      function wfvzr(t)

      include 'real8.h'

      save

      data b  /0.08459d0/,
     *     c  /-1.34d0/,
     *     cz /33897.0d0/,
     *     bz /8490.0d0/

      wfvzr = 0.0d0

      if (t.lt.0.0d0) then
         call terror(t,'wfvzr ')
      else if (t.gt.0.0d0.and.t.lt.0.263d0) then
         wfvzr = 2.5d4 * t**0.8d0
      else
         wfvzr = (1.0d0 - b * t**c) * (cz * t + bz)
      endif

      return
      end
