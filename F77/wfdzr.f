      function wfdzr(t)

      include 'real8.h'

      save

      data b  /0.03499d0/,
     *     c  /-1.068d0/,
     *     cz /33897.0d0/,
     *     bz /8490.0d0/

      wfdzr = 0.0d0

      if (t.lt.0.0) then
         call terror(t,'wfdzr ')
      else if(t.gt.0.0d0.and.t.lt.0.265d0) then
         wfdzr = 2.568d4 * t**0.395d0
      else
         wfdzr = (1.0d0 - b * t**c) * (cz * t + bz) + 500.0d0
      endif

      return
      end
