      subroutine well(t,r,depth)

      include 'real8.h'

      save

      include "wfrt.h"

      depth = odr

      rmsw  = min(0.7d0 * rzd,1.55d4)

      if (t.le.1.2d0.and.r.le.rmsw) then
         depth = max(-1.21d-3,-1.5d-3*exp(-8.0d-13 * r * r * r))

         if (t.gt.1.0d0)depth = (1.2d0 - t) * depth * 5.0d0

         trm = 0.90d0 * rmsw

         if (r.ge.trm) then
            depth = 10.0d0 * (depth * (rmsw - r) + 
     *              odr * (r - trm)) / rmsw
         endif
      endif

      return
      end
