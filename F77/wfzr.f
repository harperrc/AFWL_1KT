      function wfzr(t)

      include 'real8.h'

      save

      data b  /0.03291d0/,
     *     c  /-1.086d0/,
     *     cz /33897.0d0/,
     *     bz /8490.0d0/

      wfzr = 0.0d0

      if (t.lt.0.0d0) then
         call terror(t,'wfzr  ')
      else
         wfzr = (1.0d0 - b*t**c) * (cz * t + bz)
      endif

      return
      end
