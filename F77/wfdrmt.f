      subroutine wfdrmt(t,r)

      include 'real8.h'

      save

      include "wfrt.h"

      data ttold /0.0d0/

      rpk = prad

      if (t.ge.0.20) then
         rpk = rpk * 1.0d-5
         r   = r * 1.0d-5

         if (t.ne.ttold) then
            rz     = rzd
            rmn    = rz - 9.7163d3 * t**0.12115d0
            rz     = rz * 1.0d-5
            rmn    = rmn * 1.0d-5
            odmn   = -0.50d0 * odpk + 2.2d-5*t**(-1.6026d0)
            rneg   = rz - rmn
            rp1s   = rpk - rz
            rnp    = rpk - rmn
            a1n    = odpk / rp1s
            b1n    = odpk - a1n * rpk
            odmn1n = a1n * rmn + b1n
            fmlt   = abs(odmn / odpk)
            odmhy  = odmn + fmlt * (odmn1n - odmn)
            alpha  = (rnp / (odmhy - odpk) + rp1s / odpk) / rneg
            beta   = -rp1s * (alpha + 1.0d0 / odpk)
            fngz   = rnp / rp1s
            dnom   = alpha * rnp + beta
            bcrmn  = 1.0d0 - odmn / (rnp / dnom + odpk)
            crmn1b = log(bcrmn)
            cgz1   = (beta * (1.0 / bcrmn - 1.0d0) / dnom) /
     *               ((fngz * crmn1b * rnp**(fngz - 1.0d0)) * 
     *               (rnp + odpk * dnom))
            cgz    = exp(cgz1)
            bgz1   = crmn1b / cgz**(rnp**fngz)
            bgz    = exp(bgz1)
         endif
         
         rbr  = rpk - r
         gr   = 1.0d0 - bgz**(cgz**(rbr**fngz))
         if (r.gt.rz)gr = (rpk - r) / rp1s * gr +
     *               (r - rz) / rp1s
         hr  = rbr / (alpha * rbr + beta) + odpk
         odr = gr * hr
         rpk = rpk * 1.0d5
         r   = r * 1.0d5
      endif

      if (t.le.1.0d0) then
         if (t.ne.ttold) then
            a = -1.2d-3
            c = log(max(1.0d-20,-a/(odpk - a))) / (rzd - rpk)
            b = (odpk - a) * exp(-c * rpk)
         endif

         wflt = a + b * exp(c * r)

         if (t.le.0.20d0) then
            odr = wflt
         else
            odr = (wflt * (1.0d0 - t) + odr * (t - 0.20d0)) * 1.25d0
         endif
      endif

      ttold = t

      return
      end
